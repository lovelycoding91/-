"""배치 간 변동 시뮬레이션 (Monte Carlo).

같은 SOP로 돌려도 배치는 똑같이 나오지 않는다. 어디서 흔들리는지,
규격 탈락률이 얼마인지를 본다.
"""
import random
import statistics as stats
from typing import Dict, List

from . import config as C
from . import kinetics as K
from . import report as R
from . import stages as S


def run_batches(n: int = 40, strain_key: str = "schizophyllum",
                seed: int = 7) -> Dict[str, object]:
    rng = random.Random(seed)
    recs: List[Dict[str, float]] = []
    sp = C.STRAINS[strain_key]
    f = C.EQUIPMENT["main"]

    for i in range(n):
        # 배치마다 흔들리는 실제 변수들
        mu = sp.mu_max * rng.gauss(1.0, 0.055)          # 접종원 활성
        X0 = max(rng.gauss(0.50, 0.09), 0.15)           # 종균 농도 편차
        S0 = rng.gauss(45.0, 1.4)                       # 칭량/멸균 손실
        T = sp.T_opt + rng.gauss(0.0, 0.35)             # 온도제어 편차
        strain = C.Strain(**{**sp.__dict__, "mu_max": mu})

        contaminated = rng.random() < 0.035             # 오염 발생률 3.5%
        if contaminated:
            recs.append({"batch": i + 1, "titer": 0.0, "powder": 0.0,
                         "purity": 0.0, "released": False, "fail": "오염"})
            continue

        Rr = K.simulate_fermentation(
            strain, f, S0=S0, t_max_h=144.0, X0=X0, dt=0.5, T_C=T,
            DO_setpoint=0.30, head_pressure_barg=0.3, fed_batch=True,
            feed_start_h=48.0, rng=rng)
        pr = K.productivity(Rr)

        # 다운스트림
        mass = pr["product_in_broth_kg"]
        jj = lambda v, cv: max(rng.gauss(v, v * cv), v * 0.4)
        imp = {"protein": jj(0.16, 0.14) * mass, "hmw_nonglucan": jj(0.17, 0.11) * mass,
               "pigment": jj(0.055, 0.22) * mass, "salt": jj(0.21, 0.08) * mass,
               "residual_sugar": pr["residual_S_g_L"] * pr["broth_volume_L"] / 1000.0}
        for op in C.DOWNSTREAM:
            if op.key == "dry":      # 건조 직전 pH 조정·대이온 염 투입
                imp["salt"] = imp.get("salt", 0.0) + rng.uniform(0.020, 0.032) * mass
            res = K.apply_unitop(mass, imp, op, rng)
            mass, imp = res["mass_kg"], res["impurities"]
        dry = mass + sum(imp.get(k, 0.0) for k in
                         ("hmw_nonglucan", "protein", "pigment", "salt", "residual_sugar"))
        purity = mass / dry * 100.0
        moisture = rng.gauss(5.2, 0.6)
        powder = dry / (1.0 - moisture / 100.0)

        tip = K.tip_speed_m_s(f, max(Rr.rpm))
        shear_h = sum(1 for r in Rr.rpm if r > f.max_rpm * 0.9) * 0.5
        import math
        Mw = 2600.0 * math.exp(-0.00042 * tip ** 2 * shear_h) * 0.95 * 0.92 * rng.uniform(0.95, 1.03)

        fails = []
        if purity < 80.0:
            fails.append("순도")
        if not (sp.mw_target_kda[0] <= Mw <= sp.mw_target_kda[1]):
            fails.append("분자량")
        if moisture > 7.0:
            fails.append("수분")
        if rng.random() < 0.02:
            fails.append("미생물")

        recs.append({"batch": i + 1, "titer": pr["titer_g_L"], "powder": powder,
                     "purity": purity, "Mw": Mw, "moisture": moisture,
                     "released": not fails, "fail": ",".join(fails) or "-"})
    return {"records": recs}


def render(res: Dict[str, object], n_show: int = 12) -> str:
    recs = res["records"]
    ok = [r for r in recs if r["released"]]
    out = [R.banner("MC", "배치 간 변동 시뮬레이션 (Monte Carlo)",
                    f"동일 SOP로 {len(recs)} 배치 반복 - 재현성과 규격 탈락률 확인")]

    out.append(R.section(f"1) 배치별 결과 (앞 {n_show} 배치)"))
    out.append(R.table(["배치", "역가 g/L", "분말 kg", "순도 %", "Mw kDa", "수분 %", "판정"],
                       [[r["batch"], f"{r['titer']:.2f}", f"{r['powder']:.1f}",
                         f"{r['purity']:.1f}", f"{r.get('Mw',0):,.0f}",
                         f"{r.get('moisture',0):.1f}",
                         "출하" if r["released"] else f"탈락({r['fail']})"]
                        for r in recs[:n_show]], "<>>>>><", indent=4))

    def summ(key):
        v = [r[key] for r in ok]
        if not v:
            return ("-",) * 5
        m, sd = stats.mean(v), (stats.stdev(v) if len(v) > 1 else 0.0)
        return (f"{m:.2f}", f"{sd:.2f}", f"{sd/m*100:.1f}%", f"{min(v):.2f}", f"{max(v):.2f}")

    out.append(R.section("2) 통계 (출하 배치 기준)"))
    out.append(R.table(["지표", "평균", "표준편차", "CV", "최소", "최대"],
                       [["발효 역가 (g/L)", *summ("titer")],
                        ["분말 수량 (kg)", *summ("powder")],
                        ["순도 (%)", *summ("purity")],
                        ["분자량 (kDa)", *summ("Mw")]], "<>>>>>", indent=4))

    fails: Dict[str, int] = {}
    for r in recs:
        if not r["released"]:
            for k in r["fail"].split(","):
                fails[k] = fails.get(k, 0) + 1
    out.append(R.section("3) 탈락 원인"))
    out.append(R.table(["원인", "건수", "비율"],
                       [[k, v, f"{v/len(recs)*100:.1f}%"] for k, v in
                        sorted(fails.items(), key=lambda x: -x[1])]
                       or [["없음", 0, "0.0%"]], "<>>", indent=4))
    rate = len(ok) / len(recs) * 100
    out.append(R.kv([("출하 성공률", f"{rate:.1f} %  ({len(ok)}/{len(recs)})"),
                     ("연간 유효 생산량",
                      f"{sum(r['powder'] for r in ok)/len(recs)*C.DEFAULT_COST.batches_per_year:,.0f} kg/yr"),
                     ("해석", "오염 1건의 손실 = 배치 원가 전액 + 발효조 1주 점유. 탈락률보다 '오염률'이 돈을 좌우한다"),
                     ("개선", "분자량 탈락은 전단 이력 관리(교반 상한·펌프 선정)로, 순도 탈락은 UF DV 증량으로 대응")]))
    out.append("")
    out.append(R.sparkline([r["titer"] for r in recs], width=58, height=7,
                           label="배치별 발효 역가 (g/L) - 0 은 오염 배치"))
    return "\n".join(out)
