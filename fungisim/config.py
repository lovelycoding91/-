"""공정 파라미터 정의.

수치는 담자균 액침배양 EPS(세포외다당) 생산의 공개 문헌 범위를 바탕으로 한
'대표값'이다. 실제 공정 설계 전에는 반드시 자사 균주로 실측 보정해야 한다.
"""
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Tuple


# ----------------------------------------------------------------------------
# 균주 (Stage 0~1)
# ----------------------------------------------------------------------------
@dataclass
class Strain:
    name: str                    # 학명
    product: str                 # 목표 대사산물
    inci: str                    # 목표 INCI 명
    mu_max: float                # 최대비성장속도 [1/h]
    Ks: float                    # 탄소원 반포화상수 [g/L]
    Yxs: float                   # 균체수율 [g DCW / g 당]
    Yps: float                   # 산물수율 [g EPS / g 당]
    alpha: float                 # Luedeking-Piret 성장연계 계수 [g P/g X]
    beta: float                  # Luedeking-Piret 비성장연계 계수 [g P/g X/h]
    ms: float                    # 유지대사 계수 [g 당/g DCW/h]
    qO2_max: float               # 최대비산소소비속도 [mmol O2/g DCW/h]
    X_max: float                 # 균체 수용능력 [g DCW/L]
    T_opt: float                 # 최적온도 [degC]
    pH_opt: float                # 최적 pH
    mw_target_kda: Tuple[float, float]   # 목표 분자량 범위 [kDa]
    morphology: str              # 액침배양 형태
    visc_scale: float = 1.0      # 산물의 점도 기여 (고분자 다당 1.0, 저분자 ~0)
    notes: str = ""


STRAINS: Dict[str, Strain] = {
    # 기본 트랙: 치마버섯 -> 스키조필란(β-1,3/1,6-글루칸). 보습·진정·장벽.
    "schizophyllum": Strain(
        name="Schizophyllum commune KCTC-xxxx",
        product="Schizophyllan (β-1,3/1,6-glucan, EPS)",
        inci="Schizophyllan",
        mu_max=0.055, Ks=1.2, Yxs=0.42, Yps=0.26,
        alpha=0.35, beta=0.011, ms=0.018,
        qO2_max=3.2, X_max=16.0,
        T_opt=28.0, pH_opt=5.3,
        mw_target_kda=(450.0, 2500.0),
        morphology="pellet + filamentous (혼합형, 점도 상승 큼)",
        visc_scale=1.00,
        notes="삼중나선 구조. 고점도 비뉴턴 배양액 -> 산소전달이 율속.",
    ),
    # 대안 트랙 A: 흰목이 -> 히알루론산 대체 보습 다당
    "tremella": Strain(
        name="Tremella fuciformis TF-02",
        product="Tremella polysaccharide (glucuronoxylomannan)",
        inci="Tremella Fuciformis Polysaccharide",
        mu_max=0.070, Ks=0.9, Yxs=0.38, Yps=0.22,
        alpha=0.28, beta=0.009, ms=0.020,
        qO2_max=3.8, X_max=13.0,
        T_opt=25.0, pH_opt=5.8,
        mw_target_kda=(800.0, 1800.0),
        morphology="yeast-like 분아 + 균사 (이형성)",
        visc_scale=0.85,
        notes="효모형 단계에서 EPS 생산성 높음. HA 대체 보습소재.",
    ),
    # 대안 트랙 B: 누룩곰팡이 -> 코직산(저분자, 다운스트림 전혀 다름)
    "aspergillus": Strain(
        name="Aspergillus oryzae AO-11",
        product="Kojic acid (저분자 미백소재)",
        inci="Kojic Acid",
        mu_max=0.110, Ks=0.5, Yxs=0.45, Yps=0.30,
        alpha=0.05, beta=0.026, ms=0.025,
        qO2_max=5.5, X_max=22.0,
        T_opt=30.0, pH_opt=3.0,
        mw_target_kda=(0.142, 0.142),
        morphology="pellet (0.5-2 mm)",
        visc_scale=0.02,   # 저분자 결정성 물질 - 배양액 점도는 균체에만 좌우
        notes="저분자라 UF 대신 결정화/활성탄 정제 경로를 탄다.",
    ),
}


# ----------------------------------------------------------------------------
# 배지 (Stage 2)
# ----------------------------------------------------------------------------
@dataclass
class MediumComponent:
    name: str
    g_per_L: float
    unit_cost_usd_per_kg: float
    role: str
    carbon_frac: float = 0.0     # 당으로 환산되는 비율 (탄소원만 >0)
    nitrogen_frac: float = 0.0   # 질소 함량 비율


# 화장품 원료는 '동물유래 무첨가 / non-GMO / allergen-free' 요구가 강하다.
# -> 효모추출물 대신 식물성 펩톤을 쓰고, 글루코스는 non-GMO 옥수수당 사용.
MEDIUM_MAIN: List[MediumComponent] = [
    MediumComponent("Glucose (non-GMO corn)", 45.0, 0.85, "탄소원", carbon_frac=1.00),
    MediumComponent("Soy peptone (식물성)",     4.0, 9.50, "유기질소원", nitrogen_frac=0.105),
    MediumComponent("Yeast extract (비동물성)", 2.0, 14.00, "성장인자/비타민", nitrogen_frac=0.100),
    MediumComponent("KH2PO4",                   1.0, 2.20, "인산/완충"),
    MediumComponent("MgSO4.7H2O",               0.5, 0.95, "보조인자 Mg"),
    MediumComponent("CaCl2.2H2O",               0.1, 0.80, "세포벽 안정"),
    MediumComponent("Thiamine-HCl",            0.001, 420.0, "필수 비타민"),
    MediumComponent("Antifoam (silicone-free)", 0.25, 11.0, "소포제"),
]

MEDIUM_SEED: List[MediumComponent] = [
    MediumComponent("Glucose (non-GMO corn)", 28.0, 0.85, "탄소원", carbon_frac=1.00),
    MediumComponent("Soy peptone (식물성)",     5.0, 9.50, "유기질소원", nitrogen_frac=0.105),
    MediumComponent("Yeast extract (비동물성)", 3.0, 14.00, "성장인자"),
    MediumComponent("KH2PO4",                   1.0, 2.20, "완충"),
    MediumComponent("MgSO4.7H2O",               0.5, 0.95, "보조인자"),
]


# ----------------------------------------------------------------------------
# 설비 (Stage 3~5)
# ----------------------------------------------------------------------------
@dataclass
class Fermenter:
    tag: str
    total_volume_L: float
    working_volume_L: float
    D_T_ratio: float = 0.38      # 임펠러경/조경
    n_impeller: int = 3
    impeller: str = "Rushton x1 (하단) + 2 x pitched-blade (상단)"
    max_rpm: float = 350.0
    max_vvm: float = 1.0
    Po: float = 4.5              # 비통기 동력수 (Rushton 기준 합성값)
    kLa_coeff: float = 0.026     # kLa = C*(P/V)^a*(vs)^b
    kLa_a: float = 0.45
    kLa_b: float = 0.52
    max_P_V_W_m3: float = 2500.0


EQUIPMENT: Dict[str, Fermenter] = {
    "flask":   Fermenter("250 mL 삼각플라스크", 0.25, 0.10, max_rpm=150, max_vvm=0.0,
                         impeller="orbital shaker 150 rpm", kLa_coeff=0.015),
    "seed1":   Fermenter("30 L 시드 발효조", 30.0, 20.0, max_rpm=450,
                         max_P_V_W_m3=4000.0),
    "seed2":   Fermenter("500 L 시드 발효조", 500.0, 350.0, max_rpm=280,
                         max_P_V_W_m3=2800.0),
    "main":    Fermenter("10,000 L 생산 발효조", 10000.0, 7000.0, max_rpm=140,
                         max_P_V_W_m3=2200.0),
}

# 접종률 (v/v)
INOCULUM_RATIO = 0.08

# ----------------------------------------------------------------------------
# 멸균 (Stage 3)
# ----------------------------------------------------------------------------
@dataclass
class SterilizationSpec:
    T_hold_C: float = 121.0
    z_value_C: float = 10.0
    D121_min: float = 1.8        # G. stearothermophilus 포자 D값
    spore_load_per_L: float = 1.0e5
    target_survivors: float = 1.0e-3   # 비무균 확률 1/1000 batch
    heatup_min: float = 28.0     # 7 kL 급 증기 승온 시간
    cooldown_min: float = 42.0


# ----------------------------------------------------------------------------
# 다운스트림 (Stage 6~8)
# ----------------------------------------------------------------------------
@dataclass
class UnitOp:
    key: str
    name: str
    purpose: str
    yield_mean: float            # 목적물 회수율
    yield_sd: float
    impurity_removal: Dict[str, float] = field(default_factory=dict)
    # 불순물 제거율: {'protein':0.6, 'pigment':0.3, ...}
    duration_h: float = 2.0
    consumables_usd_per_batch: float = 0.0
    notes: str = ""


DOWNSTREAM: List[UnitOp] = [
    UnitOp("kill", "가열 실활 (80 degC / 30 min)",
           "미생물 사멸 + 자가분해효소 실활, 점도 일부 저하",
           0.995, 0.002, {"viable_cells": 0.99999}, 1.0, 40.0,
           "화장품 원료는 생균 불검출이 전제. 과열 시 다당 해중합 주의."),
    UnitOp("dilute", "희석 (1 : 1.5 정제수)",
           "비뉴턴 고점도 배양액의 펌핑/분리성 확보",
           0.999, 0.001, {}, 0.5, 18.0,
           "희석 없이 원심분리 시 균체-다당 상분리 불량."),
    UnitOp("separate", "균체 분리 (디스크스택 원심 + 규조토 여과)",
           "균사체/세포파편 제거, 목적 EPS는 상등액에 잔류",
           0.915, 0.022, {"biomass": 0.985, "pigment": 0.25}, 6.0, 310.0,
           "최대 손실 공정. 균체 케이크에 EPS가 물리적으로 포획됨."),
    UnitOp("enzyme", "효소 처리 (중성 protease, 50 degC / 4 h)",
           "잔류 단백/펩타이드 분해 -> 후단 제거 용이",
           0.985, 0.005, {"protein": 0.72, "hmw_nonglucan": 0.15}, 5.0, 260.0,
           "Sevag 법(클로로포름)은 잔류용매 이슈로 화장품에선 비선호."),
    UnitOp("uf", "한외여과/투석여과 (100 kDa MWCO, 8 DV)",
           "저분자 불순물(당/염/색소/펩타이드) 제거 + 농축",
           0.885, 0.030, {"protein": 0.88, "pigment": 0.80, "salt": 0.985,
                          "residual_sugar": 0.99, "endotoxin": 0.70}, 10.0, 520.0,
           "정제의 핵심. 분획으로 분자량 분포도 함께 제어."),
    UnitOp("carbon", "활성탄 흡착 (0.5% w/v, 45 degC / 1 h)",
           "색소/이취/내독소 추가 제거",
           0.955, 0.012, {"pigment": 0.85, "odor": 0.90, "endotoxin": 0.80}, 2.0, 150.0,
           "과량 사용 시 다당 흡착 손실 급증."),
    UnitOp("polish", "정밀여과 (0.2 um) + 에탄올 침전 (3 vol, 4 degC)",
           "무균화 + 최종 분획 회수",
           0.900, 0.028, {"protein": 0.45, "residual_sugar": 0.95, "bioburden": 0.999,
                          "hmw_nonglucan": 0.25},
           8.0, 680.0,
           "에탄올 분별침전으로 공생성 고분자 불순물을 일부 분리. 용매는 회수증류 재사용."),
    UnitOp("dry", "분무건조 (inlet 170 / outlet 85 degC)",
           "분말화, 수분 <= 7%",
           0.930, 0.018, {"water": 0.93}, 6.0, 290.0,
           "동결건조는 회수율 97%이나 kg당 비용 4~6배. 용도별 선택."),
]


# ----------------------------------------------------------------------------
# 최종 규격 (Stage 10)
# ----------------------------------------------------------------------------
@dataclass
class Spec:
    item: str
    method: str
    lo: float
    hi: float
    unit: str
    critical: bool = True


SPECS: List[Spec] = [
    Spec("베타글루칸 순도 (건조물 기준)", "Phenol-H2SO4 / 효소법(Megazyme)", 80.0, 100.0, "%"),
    Spec("평균 분자량 Mw",               "GPC-MALS (pullulan 표준)",        450.0, 2500.0, "kDa"),
    Spec("잔류 단백질",                   "Bradford / Kjeldahl",              0.0, 2.0, "%"),
    Spec("건조감량 (수분)",               "105 degC 3 h",                     0.0, 7.0, "%"),
    Spec("강열잔분 (회분)",               "600 degC",                         0.0, 6.0, "%"),
    Spec("pH (1% 수용액)",                "pH meter 25 degC",                 5.0, 7.5, "-"),
    Spec("잔류 에탄올",                   "HS-GC/FID",                        0.0, 1000.0, "ppm"),
    Spec("총호기성생균수",                "대한민국화장품안전기준 통칙",       0.0, 100.0, "CFU/g"),
    Spec("대장균/녹농균/황색포도상구균",  "특정세균 시험",                     0.0, 0.0, "검출/g"),
    Spec("내독소 (주사제급 요구 시)",      "LAL kinetic chromogenic",          0.0, 50.0, "EU/g", False),
    Spec("납 Pb",                         "ICP-MS",                           0.0, 20.0, "ppm"),
    Spec("비소 As",                       "ICP-MS",                           0.0, 10.0, "ppm"),
    Spec("수은 Hg",                       "수은분석기",                        0.0, 1.0, "ppm"),
    Spec("카드뮴 Cd",                     "ICP-MS",                           0.0, 5.0, "ppm"),
]


# ----------------------------------------------------------------------------
# 원가 (Stage 13)
# ----------------------------------------------------------------------------
@dataclass
class CostModel:
    electricity_usd_per_kWh: float = 0.11
    steam_usd_per_ton: float = 32.0
    cw_usd_per_m3: float = 0.35          # 냉각수
    wfi_usd_per_m3: float = 2.8          # 정제수
    labor_usd_per_h: float = 38.0        # 교대 2인 기준 합산
    operators: int = 2
    qc_usd_per_batch: float = 2400.0
    waste_usd_per_m3: float = 9.5
    capex_usd: float = 4_200_000.0
    depreciation_years: int = 10
    batches_per_year: int = 38
    overhead_rate: float = 0.22          # 제조간접 (직접비 대비)


DEFAULT_COST = CostModel()


# 랜덤 시드 (재현성)
DEFAULT_SEED = 20261027
