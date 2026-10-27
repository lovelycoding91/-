"""공정 물리/생물 모델.

- 멸균: Del factor(nabla) 와 F0
- 교반/통기: 비동력 P/V, 겉보기 기체유속 vs, kLa (비뉴턴 점도 보정)
- 발효: 로지스틱 성장 + Monod 기질제한 + 산소제한, Luedeking-Piret 생산
- 다운스트림: 수율 전파 + 불순물 제거 전파
"""
import math
import random
from dataclasses import dataclass, field
from typing import Dict, List, Optional


# ============================================================================
# 멸균
# ============================================================================
def del_factor(spore_load_per_L: float, volume_L: float, target_survivors: float) -> float:
    """nabla = ln(N0/N). 사멸시켜야 할 로그수."""
    N0 = spore_load_per_L * volume_L
    return math.log(N0 / target_survivors)


def kd_at(T_C: float, D121_min: float, z_C: float) -> float:
    """온도 T에서의 1차 사멸속도상수 [1/min]."""
    D_T = D121_min * 10 ** ((121.0 - T_C) / z_C)
    return 2.303 / D_T


def F0_equivalent(T_C: float, minutes: float, z_C: float = 10.0) -> float:
    """임의 온도 체류를 121 degC 등가시간(F0)으로 환산."""
    return minutes * 10 ** ((T_C - 121.0) / z_C)


def sterilization_profile(spec, volume_L: float) -> Dict[str, float]:
    """승온-유지-냉각 3구간의 사멸 기여를 적분해 유지시간을 역산."""
    nabla_req = del_factor(spec.spore_load_per_L, volume_L, spec.target_survivors)

    # 승온/냉각 구간: 100 -> 121 degC 선형 가정, 1분 간격 수치적분
    def ramp_contrib(T_start, T_end, minutes):
        n = max(int(minutes), 1)
        dt = minutes / n
        s = 0.0
        for i in range(n):
            T = T_start + (T_end - T_start) * (i + 0.5) / n
            if T < 95.0:        # 95 degC 미만은 기여 무시
                continue
            s += kd_at(T, spec.D121_min, spec.z_value_C) * dt
        return s

    nabla_up = ramp_contrib(100.0, spec.T_hold_C, spec.heatup_min)
    nabla_dn = ramp_contrib(spec.T_hold_C, 95.0, spec.cooldown_min)
    kd_hold = kd_at(spec.T_hold_C, spec.D121_min, spec.z_value_C)
    hold_min = max((nabla_req - nabla_up - nabla_dn) / kd_hold, 0.0)

    # 영양소 파괴 (티아민/글루코스-아미노산 Maillard). Ea ~ 100 kJ/mol 수준
    total_equiv = F0_equivalent(spec.T_hold_C, hold_min) + \
        F0_equivalent(110.0, spec.heatup_min + spec.cooldown_min) * 0.5
    nutrient_loss = 1.0 - math.exp(-0.0042 * total_equiv)

    return {
        "nabla_required": nabla_req,
        "nabla_heatup": nabla_up,
        "nabla_hold": nabla_req - nabla_up - nabla_dn,
        "nabla_cooldown": nabla_dn,
        "hold_min": hold_min,
        "F0_min": F0_equivalent(spec.T_hold_C, hold_min),
        "total_cycle_min": spec.heatup_min + hold_min + spec.cooldown_min,
        "nutrient_loss_frac": nutrient_loss,
        "nonsterility_risk": spec.target_survivors,
    }


# ============================================================================
# 교반 / 산소전달
# ============================================================================
def apparent_viscosity_cP(P_g_L: float, shear_rate_1s: float = 50.0,
                          K0: float = 0.0022, n_index: float = 0.38) -> float:
    """다당 농도에 따른 겉보기 점도 (power-law, 전단담화).

    mu_app = K * gamma^(n-1), K는 EPS 농도의 멱함수로 급상승.
    """
    K = K0 * (1.0 + P_g_L) ** 2.35          # Pa.s^n
    mu_Pa_s = K * shear_rate_1s ** (n_index - 1.0)
    return max(mu_Pa_s * 1000.0, 1.0)        # cP, 하한 = 물


def tank_diameter_m(working_volume_L: float, aspect_ratio: float = 2.0) -> float:
    """작업용적과 H/D 비로부터 조 내경 [m]. V = (pi/4) D^2 * (AR*D)."""
    V = working_volume_L / 1000.0
    return (4.0 * V / (math.pi * aspect_ratio)) ** (1.0 / 3.0)


def impeller_diameter_m(ferm) -> float:
    return ferm.D_T_ratio * tank_diameter_m(ferm.working_volume_L)


def tip_speed_m_s(ferm, rpm: float) -> float:
    return math.pi * impeller_diameter_m(ferm) * rpm / 60.0


def power_number_gassed(Po: float, vvm: float) -> float:
    """통기에 의한 동력 저하 Pg/P0 (간이 상관식)."""
    ratio = max(0.45, 1.0 - 0.45 * vvm ** 0.5)
    return Po * ratio


def power_per_volume(ferm, rpm: float, vvm: float, rho: float = 1020.0) -> float:
    """통기 상태 비동력 P/V [W/m3]."""
    D = impeller_diameter_m(ferm)
    N = rpm / 60.0
    Pog = power_number_gassed(ferm.Po, vvm)
    P = Pog * rho * N ** 3 * D ** 5 * ferm.n_impeller ** 0.8
    return P / (ferm.working_volume_L / 1000.0)


def rpm_for_PV(ferm, target_PV: float, vvm: float, rho: float = 1020.0) -> float:
    """목표 P/V를 내는 회전수 [rpm] (스케일업 기준 비교용)."""
    D = impeller_diameter_m(ferm)
    Pog = power_number_gassed(ferm.Po, vvm)
    V = ferm.working_volume_L / 1000.0
    N3 = target_PV * V / (Pog * rho * D ** 5 * ferm.n_impeller ** 0.8)
    return (N3 ** (1.0 / 3.0)) * 60.0


def superficial_gas_velocity(ferm, vvm: float) -> float:
    """겉보기 기체 선속도 vs [m/s]."""
    V_m3 = ferm.working_volume_L / 1000.0
    Q = vvm * V_m3 / 60.0
    d = tank_diameter_m(ferm.working_volume_L)
    return Q / (math.pi * d ** 2 / 4.0)


def kla(ferm, rpm: float, vvm: float, mu_app_cP: float) -> float:
    """kLa [1/h]. Van't Riet 형 상관식 + 비뉴턴 점도 보정."""
    PV = min(power_per_volume(ferm, rpm, vvm), ferm.max_P_V_W_m3)
    vs = superficial_gas_velocity(ferm, vvm)
    if vs <= 0:
        return 0.5
    base = ferm.kLa_coeff * (PV ** ferm.kLa_a) * (vs ** ferm.kLa_b) * 3600.0
    return max(base * (mu_app_cP ** -0.47), 0.5)


def c_star_O2(T_C: float, salinity_g_L: float = 8.0,
              head_pressure_barg: float = 0.0) -> float:
    """포화 용존산소 [mmol/L]. 가압 운전 시 Henry 법칙으로 선형 증가."""
    c_mg = 14.6 - 0.41 * T_C + 0.008 * T_C ** 2 - 0.000077 * T_C ** 3
    c_mg *= (1.0 - 0.0045 * salinity_g_L)
    c_mg *= (1.0 + head_pressure_barg / 1.013)
    return c_mg / 32.0


# ============================================================================
# 발효 (Stage 4,5 공용)
# ============================================================================
@dataclass
class FermentationResult:
    t: List[float] = field(default_factory=list)
    X: List[float] = field(default_factory=list)
    S: List[float] = field(default_factory=list)
    P: List[float] = field(default_factory=list)
    V: List[float] = field(default_factory=list)
    DO: List[float] = field(default_factory=list)
    mu: List[float] = field(default_factory=list)
    visc: List[float] = field(default_factory=list)
    kla: List[float] = field(default_factory=list)
    OTR: List[float] = field(default_factory=list)
    OUR: List[float] = field(default_factory=list)
    rpm: List[float] = field(default_factory=list)
    vvm: List[float] = field(default_factory=list)
    fO2: List[float] = field(default_factory=list)
    fed_kg: List[float] = field(default_factory=list)
    t_end: float = 0.0
    stop_reason: str = ""
    S_total_g_L: float = 0.0

    def at(self, i: int) -> Dict[str, float]:
        return {"t": self.t[i], "X": self.X[i], "S": self.S[i], "P": self.P[i],
                "V": self.V[i], "DO": self.DO[i], "mu": self.mu[i],
                "visc": self.visc[i], "kLa": self.kla[i], "rpm": self.rpm[i],
                "vvm": self.vvm[i], "fO2": self.fO2[i]}

    def final(self) -> Dict[str, float]:
        return self.at(len(self.t) - 1)


def simulate_fermentation(
    strain, ferm, S0: float, t_max_h: float, dt: float = 0.25,
    X0: float = 0.25, DO_setpoint: float = 0.30,
    T_C: Optional[float] = None,
    V0_L: Optional[float] = None,
    fed_batch: bool = False,
    feed_S_setpoint: float = 6.0,       # 유가식 기질 유지농도 [g/L]
    feed_kp: float = 0.45,              # 비례제어 게인 [1/h]
    feed_max_g_L_h: float = 2.2,
    feed_start_h: float = 48.0,
    feed_S_conc: float = 550.0,         # 공급액 당 농도 [g/L]
    feed_stop_h: Optional[float] = None,
    rng: Optional[random.Random] = None,
    head_pressure_barg: float = 0.0,
    strict_viscosity: bool = True,
    stop_on_substrate: bool = True,
    dP_stop: float = 0.0012,            # 생산속도 하한 [g/L/h]
) -> FermentationResult:
    """교반조 액침배양 수치적분 (explicit Euler, 용적변화 포함).

    상태변수: X(균체), S(기질), P(산물), V(배양액 용적)
    제어: DO 캐스케이드(rpm -> vvm), 유가식 기질 비례제어
    제한: 점도 상승 -> kLa 저하 -> 산소전달 율속 f_O2
    """
    rng = rng or random.Random(0)
    T_C = T_C if T_C is not None else strain.T_opt
    Cstar = c_star_O2(T_C, head_pressure_barg=head_pressure_barg)
    V = V0_L if V0_L is not None else ferm.working_volume_L
    V_max = ferm.total_volume_L * 0.82

    t, X, S, P = 0.0, X0, S0, 0.0
    rpm = max(ferm.max_rpm * 0.45, 1.0)
    vvm = 0.25
    fed_kg = 0.0
    R = FermentationResult()
    R.S_total_g_L = S0
    stop = "계획 배양시간 도달"

    while t <= t_max_h + 1e-9:
        mu_app = (apparent_viscosity_cP(
            P, K0=0.0022 * getattr(strain, "visc_scale", 1.0))
            if strict_viscosity else 1.0)
        kla_now = kla(ferm, rpm, vvm, mu_app)

        # --- 산소 ---
        monod = S / (strain.Ks + S) if S > 0 else 0.0
        qO2 = strain.qO2_max * (0.20 + 0.80 * monod)
        OUR_pot = qO2 * X
        OTR_max = kla_now * Cstar
        if OUR_pot > OTR_max:
            f_O2 = OTR_max / OUR_pot
            DO = 0.0
            OUR = OTR_max
        else:
            f_O2 = 1.0
            DO = max(Cstar - OUR_pot / max(kla_now, 1e-6), 0.0)
            OUR = OUR_pot
        OTR = kla_now * (Cstar - DO)
        do_pct = DO / Cstar

        # --- DO 제어 캐스케이드 ---
        if do_pct < DO_setpoint:
            if rpm < ferm.max_rpm * 0.999:
                rpm = min(rpm * 1.05, ferm.max_rpm)
            else:
                vvm = min(vvm * 1.04, ferm.max_vvm)
        elif do_pct > DO_setpoint + 0.15:
            if vvm > 0.26:
                vvm = max(vvm * 0.99, 0.25)
            else:
                rpm = max(rpm * 0.99, ferm.max_rpm * 0.30)
        if power_per_volume(ferm, rpm, vvm) > ferm.max_P_V_W_m3:
            rpm = min(rpm, rpm_for_PV(ferm, ferm.max_P_V_W_m3, vvm))

        # --- 유가식 공급 ---
        F = 0.0
        if (fed_batch and t >= feed_start_h and V < V_max
                and (feed_stop_h is None or t < feed_stop_h)):
            demand = max(feed_S_setpoint - S, 0.0) * feed_kp      # g/L/h
            demand = min(demand, feed_max_g_L_h)
            F = demand * V / feed_S_conc                          # L/h
        D = F / V                                                  # 희석률 1/h

        # --- 반응속도 ---
        logistic = max(1.0 - X / strain.X_max, 0.0)
        mu = strain.mu_max * monod * logistic * f_O2
        lp_beta = (strain.beta
                   * (0.30 + 0.70 * (1.0 - logistic))   # 정체기에 강해짐
                   * (0.25 + 0.75 * monod)
                   * (0.25 + 0.75 * f_O2))
        qp = strain.alpha * mu + lp_beta                 # g P / g X / h
        qs = (mu / strain.Yxs) + (qp / strain.Yps) + strain.ms

        dX = (mu - D) * X
        dP = qp * X - D * P
        dS = F * (feed_S_conc - S) / V - qs * X

        X = max(X + dX * dt, 1e-9)
        P = max(P + dP * dt, 0.0)
        S = max(S + dS * dt, 0.0)
        V += F * dt
        fed_kg += F * dt * feed_S_conc / 1000.0
        R.S_total_g_L += F * dt * feed_S_conc / V
        t += dt

        R.t.append(round(t, 3)); R.X.append(X); R.S.append(S); R.P.append(P)
        R.V.append(V); R.DO.append(do_pct * 100.0); R.mu.append(mu)
        R.visc.append(mu_app); R.kla.append(kla_now); R.OTR.append(OTR)
        R.OUR.append(OUR); R.rpm.append(rpm); R.vvm.append(vvm)
        R.fO2.append(f_O2); R.fed_kg.append(fed_kg)

        if stop_on_substrate and not fed_batch and S <= 0.05 and t > 24:
            stop = f"탄소원 고갈로 조기 종료 (t={t:.1f} h)"
            break
        if t > max(feed_start_h, 48.0) and qp * X < dP_stop:
            stop = f"생산속도 소진 -> 수확 판단 (t={t:.1f} h)"
            break

    R.t_end = t
    R.stop_reason = stop
    return R


def productivity(R: FermentationResult) -> Dict[str, float]:
    P, V = R.P[-1], R.V[-1]
    dt = R.t[1] - R.t[0] if len(R.t) > 1 else 0.25
    return {
        "titer_g_L": P,
        "broth_volume_L": V,
        "product_in_broth_kg": P * V / 1000.0,
        "biomass_g_L": R.X[-1],
        "residual_S_g_L": R.S[-1],
        "total_S_fed_g_L": R.S_total_g_L,
        "glucose_fed_kg": R.fed_kg[-1],
        "vol_productivity_g_L_d": P / (R.t_end / 24.0),
        "Yps_overall_g_g": P / max(R.S_total_g_L, 1e-9),
        "specific_yield_g_g": P / max(R.X[-1], 1e-9),
        "peak_viscosity_cP": max(R.visc),
        "peak_kLa_1_h": max(R.kla),
        "final_kLa_1_h": R.kla[-1],
        "min_DO_pct": min(R.DO),
        "hours_O2_limited": sum(1 for f in R.fO2 if f < 0.98) * dt,
        "mean_fO2_last_half": sum(R.fO2[len(R.fO2) // 2:]) / max(len(R.fO2) // 2, 1),
    }


# ============================================================================
# 다운스트림 수율/불순물 전파
# ============================================================================
def apply_unitop(mass_kg: float, impurities: Dict[str, float], op,
                 rng: random.Random) -> Dict[str, object]:
    y = max(min(rng.gauss(op.yield_mean, op.yield_sd), 1.0), 0.0)
    out_mass = mass_kg * y
    # 불순물은 (a) 목적물과 함께 물리적으로 손실되고 (b) 그 위에 선택적 제거가 걸린다
    imp = {}
    for k, v in impurities.items():
        carried = v * y if k not in ("endotoxin", "bioburden", "viable_cells", "odor") else v
        rem = op.impurity_removal.get(k, 0.0)
        if rem > 0.0:
            # 제거 효율도 로트마다 흔들린다 (효소 역가, 실제 DV 수, 활성탄 로트, 침전 온도)
            rem = min(max(rng.gauss(rem, rem * 0.09), 0.0), 0.999)
        imp[k] = carried * (1.0 - rem)
    return {"yield": y, "mass_kg": out_mass, "impurities": imp,
            "loss_kg": mass_kg - out_mass}
