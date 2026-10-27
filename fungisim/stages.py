"""0~13 단계 공정 구현. 각 함수는 StageResult를 반환하고 ProcessState를 갱신한다."""
import math
import random
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

from . import config as C
from . import kinetics as K
from . import report as R


@dataclass
class StageResult:
    idx: int
    title: str
    objective: str
    body: str = ""
    duration_h: float = 0.0
    cost_usd: float = 0.0
    mass_in_kg: Optional[float] = None
    mass_out_kg: Optional[float] = None
    cpp: List[Tuple[str, str, str]] = field(default_factory=list)   # (항목, 설정, 관리근거)
    risks: List[str] = field(default_factory=list)
    gate: Optional[bool] = None          # 단계 통과 판정
    gate_text: str = ""

    def render(self) -> str:
        out = [R.banner(self.idx, self.title, self.objective), self.body]
        if self.cpp:
            out.append(R.section("핵심공정변수 (CPP) / 관리기준"))
            out.append(R.table(["변수", "설정/범위", "관리 근거 (무엇이 틀어지는가)"],
                               self.cpp, indent=4))
        if self.risks:
            out.append(R.section("실패 모드 / 주의"))
            out.append(R.bullets(self.risks, "!"))
        foot = []
        if self.duration_h:
            foot.append(f"소요 {self.duration_h:,.1f} h")
        if self.cost_usd:
            foot.append(f"배치원가 기여 ${self.cost_usd:,.0f}")
        if self.mass_out_kg is not None:
            foot.append(f"목적물 {self.mass_out_kg:,.2f} kg")
        if self.gate is not None:
            foot.append(R.verdict(self.gate, self.gate_text or "다음 단계 진행",
                                  self.gate_text or "재작업"))
        if foot:
            out.append("\n  >> " + "   |   ".join(foot))
        return "\n".join(out) + "\n"


@dataclass
class ProcessState:
    strain_key: str = "schizophyllum"
    rng: random.Random = field(default_factory=lambda: random.Random(C.DEFAULT_SEED))
    strain: Any = None
    isolates: List[Dict] = field(default_factory=list)
    medium: List[C.MediumComponent] = field(default_factory=list)
    steril: Dict[str, float] = field(default_factory=dict)
    seed_train: List[Dict] = field(default_factory=list)
    ferm: Optional[K.FermentationResult] = None
    broth_kg: float = 0.0
    broth_L: float = 0.0
    mass_kg: float = 0.0
    impurities: Dict[str, float] = field(default_factory=dict)
    dsp_log: List[Dict] = field(default_factory=list)
    powder_kg: float = 0.0
    qc: List[Dict] = field(default_factory=list)
    costs: Dict[str, float] = field(default_factory=dict)
    timeline_h: Dict[str, float] = field(default_factory=dict)
    released: bool = False

    def add_cost(self, key: str, usd: float):
        self.costs[key] = self.costs.get(key, 0.0) + usd


# ============================================================================
# STAGE 0 : 균주 탐색 · 분리 · 동정
# ============================================================================
def stage0_isolation(st: ProcessState) -> StageResult:
    rng = st.rng
    st.strain = C.STRAINS[st.strain_key]
    sp = st.strain

    sites = [("강원 인제 활엽수 고사목", 8), ("전남 구례 참나무 벌채목", 7),
             ("경북 봉화 자작나무 쓰러진 줄기", 5), ("제주 곶자왈 서어나무", 4)]
    n_total = sum(n for _, n in sites)

    # 1차: 조직분리 성공 여부 (오염/비발아)
    isolates = []
    for site, n in sites:
        for i in range(n):
            ok = rng.random() > 0.28           # 분리 성공률 ~72%
            isolates.append({"id": f"SC-{len(isolates)+1:02d}", "site": site,
                             "isolated": ok})
    pure = [x for x in isolates if x["isolated"]]

    # 2차: 플레이트 생장속도 (mm/day) + 플라스크 EPS 역가 (g/L)
    for x in pure:
        x["growth_mm_d"] = max(rng.gauss(7.2, 2.1), 0.8)
        base = max(rng.gauss(2.4, 1.15), 0.05)
        x["eps_g_L"] = base * (0.6 + 0.055 * x["growth_mm_d"])
        x["dikaryon"] = rng.random() > 0.22    # 이핵균사(클램프 연결) 확인
        x["toxin_flag"] = rng.random() < 0.06  # 2차대사산물 LC-MS 이상 신호
    st.isolates = pure

    # 선발: 이핵 + 독소 음성 + 역가 상위
    cand = [x for x in pure if x["dikaryon"] and not x["toxin_flag"]]
    cand.sort(key=lambda x: x["eps_g_L"], reverse=True)
    top = cand[:3]
    best = top[0] if top else None
    if best:
        best["name"] = sp.name

    b = []
    b.append(R.section("1) 유전자원 채집 (ABS / 나고야의정서 사전대응)"))
    b.append(R.table(["채집지", "자실체 점수", "메모"],
                     [[s, n, "GPS/일자/기주수종 기록, 국내 유전자원 -> 해외 ABS 불필요"]
                      for s, n in sites], indent=4))
    b.append(R.kv([("총 자실체", f"{n_total} 점"),
                   ("채집 원칙", "동일 기주 3점 이하, 자실체 표면 70% EtOH 소독 후 내부조직 사용"),
                   ("기록 의무", "국내 유전자원 접근신고(생물다양성법) / 해외재료는 PIC·MAT 필수")]))

    b.append(R.section("2) 조직분리 -> 순수배양"))
    b.append(R.bullets([
        "자실체 내부조직 2~3 mm 절편 -> PDA + chloramphenicol 100 ppm + benomyl 2 ppm (세균/자낭균 억제)",
        "25 degC 암배양 5~7 일 -> 선단균사 이식(hyphal tip transfer) 3 회로 순수화",
        "클램프 연결(clamp connection) 현미경 확인 -> 이핵균사(dikaryon) 여부 판정",
        "단핵(monokaryon)은 생산성이 낮음 -> 교배로 이핵화하거나 폐기",
    ]))
    b.append(R.kv([("분리 시도", f"{n_total} 점"),
                   ("순수배양 확보", f"{len(pure)} 균주 (성공률 {len(pure)/n_total*100:.0f}%)"),
                   ("이핵 확인", f"{sum(1 for x in pure if x['dikaryon'])} 균주")]))

    b.append(R.section("3) 분자동정"))
    b.append(R.bullets([
        "DNA 추출(CTAB) -> ITS1F/ITS4 PCR -> Sanger 시퀀싱 -> NCBI BLAST",
        "판정 기준: ITS 상동성 >= 99%, query coverage >= 95%  (속 수준 모호 시 TEF1-alpha / RPB2 추가)",
        "최종 기탁: KCTC/KACC 기탁 -> 기탁번호 확보 (원료 등록·특허 요건)",
    ]))

    b.append(R.section("4) 스크리닝 (플레이트 -> 플라스크)"))
    rows = [[x["id"], x["site"].split()[0], f"{x['growth_mm_d']:.1f}",
             f"{x['eps_g_L']:.2f}", "O" if x["dikaryon"] else "X",
             "이상" if x["toxin_flag"] else "-"]
            for x in sorted(pure, key=lambda y: y["eps_g_L"], reverse=True)[:10]]
    b.append(R.table(["균주", "채집지", "생장 mm/d", "EPS g/L", "이핵", "독소스크리닝"],
                     rows, "<<>><<", indent=4))
    if best:
        b.append(R.kv([("선발 균주", f"{best['id']} -> {sp.name}"),
                       ("플라스크 역가", f"{best['eps_g_L']:.2f} g/L (최종 스케일업 전 기준값)"),
                       ("목적 대사산물", sp.product),
                       ("목표 INCI", sp.inci),
                       ("액침배양 형태", sp.morphology)]))

    b.append(R.section("5) 원료화 가능성 사전검토 (여기서 걸러야 돈을 아낀다)"))
    b.append(R.bullets([
        "BSL-1 여부 / 병원성 없음 (식약처·질병관리청 고위험병원체 목록 대조)",
        "독성 2차대사산물 부재 (Aspergillus 계열이면 aflatoxin/ochratoxin 필수 확인)",
        "기존 INCI 등재 여부 - 미등재면 PCPC 신규 등재 6~12개월 추가 소요",
        "중국 수출 계획 시 '已使用化妆品原料目录' 수재 여부 - 미수재면 신원료 등록(2~3년)",
    ]))

    st.timeline_h["stage0"] = 120 * 24
    gate = best is not None and best["eps_g_L"] > 1.5
    return StageResult(0, "균주 탐색 · 분리 · 동정", "야생 자실체에서 생산균주를 확보하고 법적·안전성 적격성을 확인",
                       "\n".join(b), duration_h=120 * 24, cost_usd=0.0,
                       cpp=[("분리 배지 항생제", "chloramphenicol 100 ppm", "세균 오염 시 전 균주 폐기"),
                            ("hyphal tip 계대", ">= 3 회", "부족하면 혼합균 상태로 스케일업 -> 재현성 붕괴"),
                            ("ITS 상동성", ">= 99%", "오동정 시 INCI 명칭·원료등록 전체 무효"),
                            ("클램프 연결", "현미경 확인 필수", "단핵균사는 EPS 역가 1/3 수준")],
                       risks=["야생주는 계대 중 sectoring(부채꼴 변이)으로 생산성 급락 - 반드시 초기에 균주은행 고정",
                              "해외 채집 균주는 나고야의정서 ABS 미준수 시 원료 판매 자체가 불법이 될 수 있음",
                              "기탁 전 논문/특허 공개 시 신규성 상실"],
                       gate=gate, gate_text="선발주 확정, 균주은행 단계로")


# ============================================================================
# STAGE 1 : 균주은행(MCB/WCB) 구축 및 계대 안정성
# ============================================================================
def stage1_cellbank(st: ProcessState) -> StageResult:
    rng = st.rng
    sp = st.strain

    # 10회 계대 생산성 추이 (sectoring/감쇠 모사)
    passages, prod = [], 100.0
    for i in range(1, 11):
        drift = rng.gauss(-1.6, 1.4)
        if rng.random() < 0.08:          # 변이 발생
            drift -= rng.uniform(8, 20)
        prod = max(prod + drift, 20.0)
        passages.append((i, prod))
    p10 = passages[-1][1]

    # 기준 이탈 시: 단일 콜로니 재순화 -> 재평가 (실제 현장의 표준 대응)
    recovered = None
    if p10 < 80.0:
        trials = []
        for i in range(6):
            trials.append((f"re-{i+1}", min(max(rng.gauss(94.0, 5.0), 60.0), 99.5)))
        trials.sort(key=lambda x: x[1], reverse=True)
        recovered = trials[0]

    b = []
    b.append(R.section("1) 2단 균주은행 (화장품 원료도 GMP 사고방식이 유리하다)"))
    b.append(R.table(["구분", "보존 형태", "수량", "보관", "용도"],
                     [["MCB (Master)", "10% glycerol 균사현탁 + 동결건조 앰플", "120 vial",
                       "-80 degC 2곳 분산 + LN2 기상", "WCB 제조 전용, 직접 생산 투입 금지"],
                      ["WCB (Working)", "10% glycerol 균사현탁", "300 vial",
                       "-80 degC", "배치당 1 vial 해동 -> 종균 출발"],
                      ["백업", "PDA 사면 + 멸균 증류수 보존", "40 ea",
                       "4 degC", "전원 사고 대비 이중화"]], indent=4))
    b.append(R.kv([("WCB 1 vial 당 배치 수", "1 배치 (교차오염·계대드리프트 차단)"),
                   ("WCB 소진 주기", f"300 vial / {C.DEFAULT_COST.batches_per_year} 배치/년 = 약 7.9 년"),
                   ("재구축 트리거", "WCB 잔량 20% 또는 생산성 10% 이상 이탈 시 MCB에서 재제조")]))

    b.append(R.section("2) 계대 안정성 시험 (10 계대)"))
    b.append(R.table(["계대", "상대 EPS 생산성 (%)", "판정"],
                     [[i, f"{v:.1f}", "정상" if v >= 85 else ("주의" if v >= 70 else "이탈")]
                      for i, v in passages], "<>>", indent=4))
    b.append(R.sparkline([v for _, v in passages], width=48, height=7,
                         label="계대 1 -> 10 상대 생산성 (%)"))
    b.append(R.kv([("10계대 유지율", f"{p10:.1f} %"),
                   ("판정 기준", ">= 80% 유지 시 '생산 균주로 적합'"),
                   ("운영 원칙", "생산 투입은 WCB 해동 후 3계대 이내로 제한")]))
    if recovered:
        b.append("\n    ** 기준 이탈 -> 리커버리 루틴 가동 **")
        b.append(R.bullets([
            "원인: 계대 중 sectoring(부채꼴 변이) 또는 저생산 균사 분획의 우점화",
            "대응: 원 MCB vial 재해동 -> 희석평판으로 단일 콜로니 6주 분리 -> 플라스크 역가 재평가",
            f"결과: {recovered[0]} 클론이 기준값 대비 {recovered[1]:.1f}% -> 이 클론으로 MCB 재구축",
            "교훈: 계대 안정성 시험은 '하면 좋은 것'이 아니라 균주은행 설계의 전제조건",
        ]))
        p10 = recovered[1]

    b.append(R.section("3) 균주은행 품질시험 (vial 릴리즈 항목)"))
    b.append(R.table(["시험", "방법", "기준"],
                     [["생존율", "해동 후 CFU / 균사 생장", ">= 70%"],
                      ["순도(무균성)", "TSA/SDA 배지 + 현미경", "타 미생물 불검출"],
                      ["동일성", "ITS 시퀀싱 재확인 + RAPD/ISSR 지문", "MCB와 일치"],
                      ["생산성", "플라스크 EPS 역가", "기준값 +-15% 이내"]], indent=4))

    st.timeline_h["stage1"] = 60 * 24
    return StageResult(1, "균주은행(MCB/WCB) 구축 · 계대 안정성",
                       "배치 간 재현성의 근원을 고정 - 여기를 건너뛰면 3년 뒤 역가가 반토막 난다",
                       "\n".join(b), duration_h=60 * 24, cost_usd=0.0,
                       cpp=[("동결보호제", "glycerol 10% (v/v)", "농도 부족 시 빙결정으로 균사 치사"),
                            ("동결 속도", "-1 degC/min 프로그램 프리징", "급속동결은 생존율 30% 이하"),
                            ("WCB 계대수", "MCB로부터 <= 3 계대", "계대 누적이 sectoring 변이의 주원인"),
                            ("보관 이중화", "물리적으로 분리된 2개소", "냉동고 1대 고장 = 사업 중단")],
                       risks=["-80 degC 냉동고 단일 보관은 단일장애점. 알람·백업전원·원격모니터링 필수",
                              "담자균은 장기 보존 중 이핵 -> 단핵 해리가 일어날 수 있음 (클램프 재확인 필요)"],
                       gate=p10 >= 80.0,
                       gate_text=(f"계대안정성 {p10:.0f}% - 생산 균주 적격"
                                  + (" (재순화 클론으로 MCB 재구축)" if recovered else "")))


# ============================================================================
# STAGE 2 : 배지 설계
# ============================================================================
def stage2_medium(st: ProcessState) -> StageResult:
    st.medium = C.MEDIUM_MAIN
    sp = st.strain
    main_V = C.EQUIPMENT["main"].working_volume_L

    rows, total_cost, C_g, N_g = [], 0.0, 0.0, 0.0
    for m in st.medium:
        kg = m.g_per_L * main_V / 1000.0
        cost = kg * m.unit_cost_usd_per_kg
        total_cost += cost
        C_g += m.g_per_L * m.carbon_frac * 0.40     # 글루코스 탄소함량 40%
        N_g += m.g_per_L * m.nitrogen_frac
        rows.append([m.name, f"{m.g_per_L:g}", m.role, f"{kg:,.1f}",
                     f"${cost:,.0f}", f"${m.unit_cost_usd_per_kg:,.2f}"])
    CN = C_g / max(N_g, 1e-9)

    b = []
    b.append(R.section("1) 배지 조성 (7,000 L 생산 배치 기준)"))
    b.append(R.table(["성분", "g/L", "역할", "투입 kg", "금액", "단가/kg"],
                     rows, "<><>>>", indent=4))
    b.append(R.kv([("배지 원재료비", f"${total_cost:,.0f} / batch  (= ${total_cost/main_V:.3f} per L broth)"),
                   ("탄소 C", f"{C_g:.2f} g/L"), ("질소 N", f"{N_g:.3f} g/L"),
                   ("C/N 비", f"{CN:.1f}  (다당 생산은 고 C/N 25~45 구간에서 유리)"),
                   ("초기 pH", f"{sp.pH_opt}  (멸균 전 조정, 멸균 후 0.2~0.3 하강 보정)")]))

    b.append(R.section("2) 화장품 원료 특유의 배지 제약 (발효산업과 다른 지점)"))
    b.append(R.table(["요구사항", "금지/대체", "이유"],
                     [["동물유래 무첨가", "beef extract/casein -> soy peptone", "vegan·할랄 인증, BSE/TSE 리스크"],
                      ["non-GMO", "GMO 옥수수당 -> non-GMO 전분당", "EU·국내 클린뷰티 바이어 요구"],
                      ["알레르겐 free", "밀/대두 유래 시 표시의무 검토", "대두펩톤도 잔류단백 0 증명 필요"],
                      ["실리콘 소포제 금지", "silicone -> 식물성 지방산 에스터", "실리콘 잔류 시 INCI 불일치"],
                      ["중금속 저감 원료", "공업용 염 -> 화장품/식품등급", "최종 Pb<=20 ppm 규격 역산"]], indent=4))

    b.append(R.section("3) 배지 최적화 DoE 이력 (실제로는 여기서 역가가 2~3배 갈린다)"))
    b.append(R.table(["단계", "설계", "인자", "결과"],
                     [["1차 스크리닝", "Plackett-Burman (12 run)", "탄소원 4종 x 질소원 4종 x 미량 3종",
                       "유의인자 = 글루코스, 대두펩톤, KH2PO4"],
                      ["2차 경로탐색", "Steepest ascent", "유의 3인자", "글루코스 30 -> 45 g/L 방향"],
                      ["3차 최적화", "CCD / RSM (20 run)", "3인자 5수준",
                       f"최적점 glucose 45 / peptone 4.0 / KH2PO4 1.0 g/L"],
                      ["검증", "3 반복 플라스크", "-", "예측 대비 오차 7% 이내"]], indent=4))
    b.append(R.bullets([
        "탄소원 50 g/L 초과 시 기질저해 + 점도 폭증으로 산소전달이 무너진다 -> 초과분은 유가식으로 넘긴다",
        "질소 고갈 타이밍이 EPS 전환의 스위치. 질소를 일찍 비우면 균체는 작지만 비생산성이 오른다",
    ]))

    st.add_cost("배지 원재료", total_cost)
    st.timeline_h["stage2"] = 45 * 24
    return StageResult(2, "배지 설계 및 최적화",
                       "역가·원가·규제적합성을 동시에 만족하는 조성 확정",
                       "\n".join(b), duration_h=45 * 24, cost_usd=total_cost,
                       cpp=[("초기 글루코스", "45 +- 2 g/L", "과잉 시 기질저해 및 점도 과다"),
                            ("C/N 비", f"{CN:.0f} +- 5", "낮으면 균체만 자라고 EPS 전환 안 됨"),
                            ("초기 pH", f"{sp.pH_opt} +- 0.2", "pH 4 이하에서 EPS 합성효소 활성 급감"),
                            ("소포제", "0.25 g/L, 식물성", "과량 시 산소전달계수 20~30% 저하")],
                       risks=["대두펩톤 로트 간 아미노산 조성 편차 -> 배치 역가 변동의 숨은 주범. 로트 고정·입고 시험 필요",
                              "원료 단가 중 질소원이 비중 큼. 대체원(옥수수침지액 CSL)은 싸지만 색소·중금속 부담 증가"])


# ============================================================================
# STAGE 3 : 멸균
# ============================================================================
def stage3_sterilization(st: ProcessState) -> StageResult:
    spec = C.SterilizationSpec()
    V = C.EQUIPMENT["main"].working_volume_L
    prof = K.sterilization_profile(spec, V)
    st.steril = prof

    # 유틸리티
    steam_ton = V / 1000.0 * 0.21 * (prof["total_cycle_min"] / 60.0) + 0.9
    cw_m3 = V / 1000.0 * 2.4
    cost = steam_ton * C.DEFAULT_COST.steam_usd_per_ton + cw_m3 * C.DEFAULT_COST.cw_usd_per_m3

    b = []
    b.append(R.section("1) 멸균 설계 - Del factor 기반 (시간을 '감'으로 잡지 않는다)"))
    b.append(R.kv([
        ("배양액 용적", f"{V:,.0f} L"),
        ("초기 포자부하 N0", f"{spec.spore_load_per_L:.0e} /L x {V:,.0f} L = {spec.spore_load_per_L*V:.2e} 개"),
        ("목표 생존수 N", f"{spec.target_survivors:.0e} (비무균 확률 1/1,000 batch)"),
        ("필요 Del factor", f"nabla = ln(N0/N) = {prof['nabla_required']:.1f}"),
        ("기준 포자", f"G. stearothermophilus  D121 = {spec.D121_min} min, z = {spec.z_value_C} degC"),
    ]))
    b.append(R.table(["구간", "조건", "시간", "nabla 기여", "비고"],
                     [["승온", "100 -> 121 degC", f"{spec.heatup_min:.0f} min",
                       f"{prof['nabla_heatup']:.1f}", "대형조일수록 기여 큼 - 무시하면 과멸균"],
                      ["유지", "121 degC", f"{prof['hold_min']:.1f} min",
                       f"{prof['nabla_hold']:.1f}", "역산된 값. 이것이 설계 답"],
                      ["냉각", "121 -> 95 degC", f"{spec.cooldown_min:.0f} min",
                       f"{prof['nabla_cooldown']:.1f}", "냉각 지연도 멸균에 기여하나 영양소는 파괴"],
                      ["합계", "-", f"{prof['total_cycle_min']:.0f} min",
                       f"{prof['nabla_required']:.1f}", f"F0 = {prof['F0_min']:.1f} min"]], indent=4))
    b.append(R.kv([("영양소 손실(추정)", f"{prof['nutrient_loss_frac']*100:.1f} %  (티아민 열분해 + Maillard 갈변)"),
                   ("과멸균 시", "환원당-아미노산 갈변 -> 색소 증가 -> 후단 활성탄 부하 급증"),
                   ("대안", "연속멸균(CST, 140 degC 30~60 s)으로 F0 유지 + 영양소 손실 1/3")]))

    b.append(R.section("2) 별도 멸균 항목 (같이 넣으면 망가지는 것들)"))
    b.append(R.table(["항목", "방식", "이유"],
                     [["글루코스 용액", "별도 121 degC 20 min 후 무균이송",
                       "인산염·아미노산과 동시 가열 시 Maillard 갈변"],
                      ["티아민 등 비타민", "0.2 um 제균여과", "열에 완전 파괴"],
                      ["유가식 공급액 (550 g/L)", "121 degC 30 min, 고농도 단독", "고농도 당은 캐러멜화 주의"],
                      ["공급 공기", "0.2 um 소수성 멤브레인 2단", "1차 파손 시 배치 전손 - 반드시 2단"],
                      ["배기", "0.2 um + 히터(결로방지)", "결로 시 필터 막힘 -> 조 가압 사고"],
                      ["배관/밸브", "SIP 121 degC 30 min, 데드레그 <= 3D", "데드레그가 오염의 1순위"]], indent=4))

    b.append(R.section("3) 무균 검증"))
    b.append(R.bullets([
        "멸균 후 12 h 무접종 홀드 -> 현미경 + TSA/SDA 플레이팅으로 무균 확인 후 접종",
        "열전대(TC) 최저온점(조 하부 배플 뒤) 실측으로 F0 계산 - 설계값이 아니라 실측값으로 릴리즈",
        "연 1회 배지충전시험(media fill)에 준하는 공정 시뮬레이션 권장",
    ]))

    st.add_cost("멸균 유틸리티", cost)
    st.timeline_h["stage3"] = prof["total_cycle_min"] / 60.0 + 12
    return StageResult(3, "멸균 (배지 · 설비 · 공기)",
                       "비무균 확률을 1/1,000 배치 이하로 설계하되 영양소 파괴는 최소화",
                       "\n".join(b), duration_h=prof["total_cycle_min"] / 60.0 + 12,
                       cost_usd=cost,
                       cpp=[("유지온도", "121 +- 1 degC", "1 degC 낮으면 사멸속도 26% 감소"),
                            ("유지시간", f"{prof['hold_min']:.0f} min (실측 F0 >= {prof['F0_min']:.0f})",
                             "F0 미달 = 오염, 과다 = 갈변·역가손실"),
                            ("공기필터", "0.2 um x 2 단, 매 배치 전 integrity test",
                             "필터 1점 파손으로 7 kL 전손"),
                            ("냉각 종점", "28 degC, 양압 유지", "음압 발생 시 비멸균 공기 역류")],
                       risks=["대형조 승온/냉각 기여를 빼먹고 유지시간만 잡으면 2배 과멸균 -> 색소·수율 손실",
                              "데드레그·샘플포트·소포제 라인이 실제 오염의 대부분"],
                       gate=True, gate_text=f"F0 {prof['F0_min']:.0f} min 확보")


# ============================================================================
# STAGE 4 : 종균 확대 배양 (Seed train)
# ============================================================================
def stage4_seed_train(st: ProcessState) -> StageResult:
    rng, sp = st.rng, st.strain
    main = C.EQUIPMENT["main"]

    steps = [
        ("WCB vial 해동", "1 mL", "-", 0.0, "-80 degC -> 37 degC 수조 2 min, 즉시 이식"),
        ("PDA 사면 활성화", "1 plate", "25 degC", 120.0, "균사 전면 피복 확인, 오염 육안검사"),
        ("1차 플라스크", "250 mL / 100 mL", "28 degC, 150 rpm", 96.0, "균사 파쇄(blender 10 s)로 접종원 균질화"),
        ("2차 플라스크", "2 L x 5 / 800 mL", "28 degC, 150 rpm", 72.0, "pellet 크기 1~2 mm 균일화"),
        ("Seed-1 (30 L)", f"{C.EQUIPMENT['seed1'].working_volume_L:.0f} L", "28 degC, 0.5 vvm", 72.0, "DO >= 40%, pH 자유"),
        ("Seed-2 (500 L)", f"{C.EQUIPMENT['seed2'].working_volume_L:.0f} L", "28 degC, 0.6 vvm", 60.0, "무균이송 전 현미경·플레이팅 2중 확인"),
    ]

    def logistic_batch(X0: float, mu: float, Xmax: float, hours: float) -> float:
        """진탕배양(플라스크)·고체배지는 산소전달 모델 없이 로지스틱으로 충분."""
        if X0 <= 0:
            return 0.0
        return Xmax / (1.0 + ((Xmax - X0) / X0) * math.exp(-mu * hours))

    # 각 단계 균체 농도 시뮬레이션
    train, X = [], 0.05
    for name, vol, cond, hrs, note in steps:
        if hrs == 0:
            train.append([name, vol, cond, "-", "-", note]); continue
        if "사면" in name:
            # 고체배지: 균체농도 대신 피복률로 관리. 다음 단계 접종량만 전달
            train.append([name, vol, cond, f"{hrs:.0f} h", "피복 100%", note])
            X = 0.06
            continue
        if "플라스크" in name:
            # 진탕배양: mu는 통기 제한으로 액침 대비 85% 수준, Xmax 는 7 g/L
            Xf = logistic_batch(X, sp.mu_max * 0.85, 7.0, hrs) * rng.uniform(0.92, 1.06)
        else:
            key = "seed1" if "30 L" in name else "seed2"
            f = C.EQUIPMENT[key]
            Rr = K.simulate_fermentation(sp, f, S0=28.0, t_max_h=hrs, X0=X,
                                         dt=0.5, DO_setpoint=0.40, rng=rng,
                                         head_pressure_barg=0.0,
                                         stop_on_substrate=False, dP_stop=0.0)
            Xf = Rr.X[-1]
        train.append([name, vol, cond, f"{hrs:.0f} h", f"{Xf:.2f} g/L", note])
        X = Xf * C.INOCULUM_RATIO            # 다음 단계 초기 균체 (접종률 8%)
    st.seed_train = train

    X0_main = X
    vol_chain = [0.1, 0.8, C.EQUIPMENT['seed1'].working_volume_L,
                 C.EQUIPMENT['seed2'].working_volume_L, main.working_volume_L]

    b = []
    b.append(R.section("1) 종균 확대 경로 (vial 1 mL -> 7,000 L, 배수 7,000,000)"))
    b.append(R.table(["단계", "용적", "조건", "배양시간", "최종 균체", "포인트"], train, indent=4))
    b.append(R.kv([("접종률", f"{C.INOCULUM_RATIO*100:.0f} % (v/v) 일정 유지"),
                   ("본배양 초기 균체", f"{X0_main:.2f} g DCW/L"),
                   ("종균 총 소요시간", f"{sum(s[3] for s in steps):.0f} h = {sum(s[3] for s in steps)/24:.1f} 일"),
                   ("확대 배수", " -> ".join(f"{v:g} L" for v in vol_chain))]))

    b.append(R.section("2) 균사체 특유의 문제 - 접종원 '형태'가 본배양을 결정한다"))
    b.append(R.table(["형태", "발생 조건", "본배양 영향"],
                     [["Pellet (1~2 mm)", "낮은 전단, 낮은 접종농도", "점도 낮음·산소전달 유리 / pellet 내부 괴사 가능"],
                      ["Filamentous (분산균사)", "높은 전단, 높은 접종농도", "EPS 역가 높음 / 점도 폭증 - 산소전달 붕괴"],
                      ["Clump (뭉침)", "접종원 파쇄 불충분", "불균일 배양·재현성 최악 - 반드시 회피"]], indent=4))
    b.append(R.bullets([
        "본 공정은 'pellet 우세 + 일부 분산' 혼합형을 목표. 2차 플라스크에서 pellet 직경 1~2 mm 로 고정",
        "접종원 균질화(균사 파쇄)를 생략하면 배치 간 역가 편차가 +-30%까지 벌어진다",
        "종균 단계 배양시간 초과(과숙) -> 자가분해 시작 -> 본배양 lag 연장",
    ]))

    b.append(R.section("3) 이송 전 릴리즈 판정"))
    b.append(R.table(["항목", "기준", "미달 시"],
                     [["현미경 검경", "타 균체·세균 불검출, 클램프 연결 확인", "폐기 - 본배양 투입 금지"],
                      ["무균시험(TSA/SDA)", "24 h 무증식", "폐기"],
                      ["균체농도 (DCW)", "설정값 +-20%", "배양시간 조정 또는 폐기"],
                      ["pH", "4.5 ~ 5.8", "과도 하강 시 자가분해 의심"],
                      ["잔당", ">= 3 g/L", "고갈 시 과숙 - 본배양 lag 연장"]], indent=4))

    dur = sum(s[3] for s in steps)
    # 종균 420 h 중 대부분은 무인 배양. 실제 작업(해동·이식·파쇄·샘플링·이송)만 과금
    attended_h = 46.0
    labor = attended_h * C.DEFAULT_COST.labor_usd_per_h
    seed_med = sum(m.g_per_L * (C.EQUIPMENT['seed1'].working_volume_L +
                                C.EQUIPMENT['seed2'].working_volume_L) / 1000.0
                   * m.unit_cost_usd_per_kg for m in C.MEDIUM_SEED)
    st.add_cost("종균 배지", seed_med)
    seed_util = (C.EQUIPMENT['seed1'].working_volume_L + C.EQUIPMENT['seed2'].working_volume_L) \
        / 1000.0 * 1.2 * dur * C.DEFAULT_COST.electricity_usd_per_kWh
    st.add_cost("종균 인건/유틸", labor + seed_util)
    st._X0_main = X0_main
    st.timeline_h["stage4"] = dur
    return StageResult(4, "종균 확대 배양 (Seed train)",
                       "1 mL 바이알을 7,000 L 접종원으로 - 형태·순도·활성을 유지한 채 700만 배 확대",
                       "\n".join(b), duration_h=dur, cost_usd=seed_med + labor + seed_util,
                       cpp=[("접종률", f"{C.INOCULUM_RATIO*100:.0f} % v/v +- 1%", "낮으면 lag 연장·오염 위험, 높으면 점도 조기 상승"),
                            ("종균 배양시간", "각 단계 대수기 후반", "과숙 시 자가분해, 미숙 시 lag"),
                            ("pellet 직경", "1~2 mm", "> 3 mm는 내부 괴사, < 0.5 mm는 점도 폭증"),
                            ("이송 무균성", "현미경 + 배지 2중 확인", "본배양 오염은 7 kL 전손")],
                       risks=["종균 1단계 오염은 금액은 작지만 생산 스케줄 전체를 1주일 밀어낸다",
                              "무균이송 라인 SIP 누락이 오염의 최대 원인"],
                       gate=True, gate_text=f"본배양 접종 준비 완료 (X0={X0_main:.2f} g/L)")


# ============================================================================
# STAGE 5 : 본 발효
# ============================================================================
def stage5_fermentation(st: ProcessState) -> StageResult:
    sp, f = st.strain, C.EQUIPMENT["main"]
    X0 = getattr(st, "_X0_main", 1.1)
    t_harvest = 144.0

    Rr = K.simulate_fermentation(
        sp, f, S0=45.0, t_max_h=t_harvest, X0=X0, dt=0.25,
        DO_setpoint=0.30, head_pressure_barg=0.3,
        fed_batch=True, feed_start_h=48.0, feed_S_setpoint=6.0,
        feed_max_g_L_h=2.2, feed_S_conc=550.0, rng=st.rng)
    st.ferm = Rr
    pr = K.productivity(Rr)

    st.broth_L = pr["broth_volume_L"]
    st.mass_kg = pr["product_in_broth_kg"]
    st.broth_kg = st.broth_L * 1.03
    # 불순물 인벤토리 (목적물 대비 상대량, kg)
    j = lambda m, cv: max(st.rng.gauss(m, m * cv), m * 0.4)   # 로트 편차
    st.impurities = {
        "biomass": pr["biomass_g_L"] * st.broth_L / 1000.0,
        "protein": j(0.16, 0.14) * st.mass_kg,
        # UF를 통과하지 못하는 공생성 고분자 (헤테로다당·당단백 복합체).
        # 분자량이 목적물과 겹쳐 막분획으로 분리되지 않는다 -> 최종 순도의 실질 상한.
        "hmw_nonglucan": j(0.17, 0.11) * st.mass_kg,
        "pigment": j(0.055, 0.22) * st.mass_kg,
        "salt": j(0.21, 0.08) * st.mass_kg,
        "residual_sugar": pr["residual_S_g_L"] * st.broth_L / 1000.0,
        "endotoxin": 1.0,
        "water": st.broth_L * 1.0,
        "bioburden": 1.0,
        "viable_cells": 1.0,
        "odor": 1.0,
    }

    # 시간별 프로파일 추출
    marks = [0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 144]
    idx = {m: min(int(m / 0.25), len(Rr.t) - 1) for m in marks}
    rows = []
    for m in marks:
        d = Rr.at(idx[m])
        phase = ("유도기" if m < 12 else "대수기" if m < 60 else
                 "전이기" if m < 84 else "정체기/생산기")
        rows.append([f"{m:>3}", phase, f"{d['X']:.2f}", f"{d['S']:.2f}", f"{d['P']:.2f}",
                     f"{d['DO']:.0f}", f"{d['visc']:.0f}", f"{d['kLa']:.0f}",
                     f"{d['rpm']:.0f}", f"{d['vvm']:.2f}", f"{d['fO2']:.2f}"])

    b = []
    b.append(R.section("1) 운전 설정"))
    b.append(R.kv([("발효조", f"{f.tag}  (작업용적 {f.working_volume_L:,.0f} L, H/D 2.0)"),
                   ("임펠러", f"{f.impeller}  D_i = {K.impeller_diameter_m(f):.3f} m"),
                   ("온도 / pH", f"{sp.T_opt} degC / pH {sp.pH_opt} (NaOH-인산 자동제어)"),
                   ("DO 설정", "30% 이상 캐스케이드 (rpm 1차 -> vvm 2차)"),
                   ("조 압력", "0.3 barg 가압 -> 포화DO 30% 증대"),
                   ("운전 방식", "회분 48 h -> 유가식 (잔당 6 g/L 비례제어, 공급액 550 g/L)"),
                   ("최대 교반", f"{f.max_rpm:.0f} rpm (tip speed {K.tip_speed_m_s(f, f.max_rpm):.2f} m/s)"),
                   ("최대 P/V", f"{f.max_P_V_W_m3:,.0f} W/m3")]))

    b.append(R.section("2) 배양 프로파일 (12 h 간격)"))
    b.append(R.table(["t(h)", "단계", "균체 X g/L", "잔당 S g/L", "EPS P g/L",
                      "DO %", "점도 cP", "kLa 1/h", "rpm", "vvm", "f_O2"],
                     rows, "<<>>>>>>>>>", indent=2))
    b.append("")
    b.append(R.sparkline(Rr.P, width=60, height=8, label="EPS 역가 P [g/L] vs 시간 (0 -> 144 h)"))
    b.append(R.sparkline(Rr.DO, width=60, height=6, label="용존산소 DO [%] - 96 h 이후 0% 고착 = 산소전달 율속", lo=0, hi=100))
    b.append(R.sparkline(Rr.visc, width=60, height=6, label="겉보기 점도 [cP] - 비뉴턴 전단담화 배양액"))

    b.append(R.section("3) 결과 요약"))
    b.append(R.table(["지표", "값", "해석"],
                     [["최종 역가", f"{pr['titer_g_L']:.2f} g/L", "문헌 최적화 범위(8~15 g/L) 상단"],
                      ["최종 균체", f"{pr['biomass_g_L']:.2f} g DCW/L", f"수용능력 {sp.X_max} g/L 대비 {pr['biomass_g_L']/sp.X_max*100:.0f}%"],
                      ["배양액 용적", f"{pr['broth_volume_L']:,.0f} L", f"유가식으로 {f.working_volume_L:,.0f} -> {pr['broth_volume_L']:,.0f} L 증가"],
                      ["조 내 목적물", f"{pr['product_in_broth_kg']:.1f} kg", "다운스트림 입력 질량"],
                      ["용적생산성", f"{pr['vol_productivity_g_L_d']:.2f} g/L/d", "담자균 EPS로는 양호 (통상 0.8~2.0)"],
                      ["총 당 투입", f"{pr['total_S_fed_g_L']:.0f} g/L 상당 ({pr['glucose_fed_kg']:.0f} kg 유가식)", "원가의 핵심 인자"],
                      ["총괄 수율 Yp/s", f"{pr['Yps_overall_g_g']:.3f} g/g", "유지대사 비중이 커 이론수율의 40% 수준"],
                      ["최고 점도", f"{pr['peak_viscosity_cP']:.0f} cP", "이 값이 다운스트림 난이도를 결정"],
                      ["kLa 변화", f"{pr['peak_kLa_1_h']:.0f} -> {pr['final_kLa_1_h']:.0f} 1/h", "점도 상승으로 60% 이상 감소"],
                      ["산소제한 시간", f"{pr['hours_O2_limited']:.0f} h", f"후반 평균 f_O2 = {pr['mean_fO2_last_half']:.2f}"],
                      ["종료 사유", Rr.stop_reason, "생산속도 < 0.0012 g/L/h 도달"]], indent=4))

    b.append(R.section("4) 이 공정의 진짜 병목 - 산소전달"))
    b.append(R.bullets([
        f"EPS가 쌓일수록 점도가 1 -> {pr['peak_viscosity_cP']:.0f} cP 로 오르고 kLa 는 {pr['peak_kLa_1_h']:.0f} -> {pr['final_kLa_1_h']:.0f} 1/h 로 떨어진다",
        "교반을 더 올리면? tip speed 5 m/s 초과 시 균사 전단파괴 + P/V 상한 + 다당 해중합(분자량 저하)",
        "통기를 더 올리면? 포말 폭증 -> 소포제 과다 -> kLa 추가 저하 (악순환)",
        "현실적 해법 3가지: (a) 산소부화 공기 (b) 조 가압 0.3->0.5 barg (c) 유가식 속도를 낮춰 점도 상승 지연",
        "분자량 목표가 높을수록(>1,000 kDa) 전단을 못 올리므로 역가를 희생해야 한다 - 역가와 분자량은 트레이드오프",
    ]))

    b.append(R.section("5) 수확 판단 (언제 멈출 것인가)"))
    b.append(R.table(["지표", "수확 트리거", "현재"],
                     [["EPS 생산속도", "< 0.05 g/L/d", f"{(Rr.P[-1]-Rr.P[-int(24/0.25)])/1.0:.3f} g/L/d (최근 24 h)"],
                      ["잔당", "유가식 중단 후 < 2 g/L", f"{pr['residual_S_g_L']:.1f} g/L"],
                      ["점도", "> 150 cP 면 분리 공정 한계", f"{pr['peak_viscosity_cP']:.0f} cP"],
                      ["자가분해 신호", "유리 단백/핵산 상승, 점도 급락", "미검출"],
                      ["경제성", "추가 24 h 비용 > 추가 EPS 가치", "균형점 통과"]], indent=4))

    dur = Rr.t_end
    # 교반 동력 (시간별 실제 rpm/vvm 적분)
    stir_kWh = sum(K.power_per_volume(f, Rr.rpm[i], Rr.vvm[i]) / 1000.0
                   * (f.working_volume_L / 1000.0) * 0.25
                   for i in range(len(Rr.t))) / 0.92          # 모터효율
    # 공기 압축 동력: 등온압축 W = p1*Q*ln(p2/p1), 2.0 bar abs 토출, 압축기 효율 0.65
    air_kWh = 0.0
    for i in range(len(Rr.t)):
        Q = Rr.vvm[i] * f.working_volume_L / 1000.0 / 60.0     # m3/s
        W = 101325.0 * Q * math.log(2.0 / 1.0) / 0.65          # W
        air_kWh += W / 1000.0 * 0.25
    # 냉각: 대사열 ~ 0.46 kW per kg O2/h 소비 기준 -> 칠러 COP 3.2
    o2_kg = sum(Rr.OUR[i] * f.working_volume_L / 1000.0 * 32.0 * 0.25
                for i in range(len(Rr.t))) / 1000.0            # kg O2
    # 대사열 ~460 kJ per mol O2 = 4.0 kWh(th)/kg O2, 칠러 COP 3.2
    cool_kWh = o2_kg * 4.0 / 3.2
    util = (stir_kWh + air_kWh + cool_kWh) * C.DEFAULT_COST.electricity_usd_per_kWh
    feed_cost = pr["glucose_fed_kg"] * 0.85
    labor = dur * C.DEFAULT_COST.labor_usd_per_h * 0.4
    st.add_cost("본발효 유틸리티(교반·통기·냉각)", util)
    st._ferm_util = {"교반 kWh": stir_kWh, "공기압축 kWh": air_kWh,
                     "냉각 kWh": cool_kWh, "소비 O2 kg": o2_kg}
    st.add_cost("유가식 당 공급", feed_cost)
    st.add_cost("본발효 인건", labor)
    st.timeline_h["stage5"] = dur
    return StageResult(5, "본 발효 (10 kL 교반조 액침배양)",
                       "목표 역가·목표 분자량·다운스트림이 감당 가능한 점도, 이 셋을 동시에 맞춘다",
                       "\n".join(b), duration_h=dur, cost_usd=util + feed_cost + labor,
                       mass_out_kg=st.mass_kg,
                       cpp=[("온도", f"{sp.T_opt} +- 0.5 degC", "+2 degC 시 저분자화, -2 degC 시 역가 저하"),
                            ("pH", f"{sp.pH_opt} +- 0.2", "pH 4 이하 EPS 합성 정지, 6 이상 단백 분비 증가"),
                            ("DO", ">= 30% (후반 불가피하게 0%)", "완전 혐기 지속 시 균사 자가분해"),
                            ("tip speed", "<= 5 m/s", "초과 시 다당 기계적 해중합 -> 분자량 규격 이탈"),
                            ("유가식 속도", "잔당 6 +- 2 g/L 유지", "과공급 시 점도 폭증, 부족 시 유지대사로 당 낭비"),
                            ("배양시간", "144 h +- 6 h", "과도 연장 시 자가분해로 단백 불순물 급증")],
                       risks=["포말 과다 -> 배기필터 침수 -> 조 가압 -> 안전밸브 작동 (실제로 가장 흔한 사고)",
                              "고점도에서 pH 전극 응답 지연 -> 실제 pH와 측정값 괴리",
                              "산소제한이 길어지면 혐기 대사산물(에탄올·유기산)로 색·이취 발생"],
                       gate=pr["titer_g_L"] >= 8.0,
                       gate_text=f"역가 {pr['titer_g_L']:.1f} g/L - 수확 승인")


# ============================================================================
# 다운스트림 공용 실행기
# ============================================================================
def _run_unitops(st: ProcessState, keys: List[str]) -> Tuple[List[List[str]], float, float]:
    rows, cost, hours = [], 0.0, 0.0
    for key in keys:
        op = next(o for o in C.DOWNSTREAM if o.key == key)
        before = st.mass_kg
        res = K.apply_unitop(before, st.impurities, op, st.rng)
        st.mass_kg = res["mass_kg"]
        st.impurities = res["impurities"]
        cost += op.consumables_usd_per_batch
        hours += op.duration_h
        st.dsp_log.append({"op": op.name, "in": before, "out": st.mass_kg,
                           "yield": res["yield"], "cum": st.mass_kg})
        rows.append([op.name, f"{before:.2f}", f"{res['yield']*100:.1f}%",
                     f"{st.mass_kg:.2f}", f"{res['loss_kg']:.2f}",
                     f"{op.duration_h:.1f}", f"${op.consumables_usd_per_batch:,.0f}"])
        st.add_cost("다운스트림 소모품/유틸", op.consumables_usd_per_batch)
    return rows, cost, hours


def _solid_composition(st: ProcessState) -> Dict[str, float]:
    p = st.mass_kg
    imp = {k: st.impurities.get(k, 0.0) for k in
           ("hmw_nonglucan", "protein", "pigment", "salt", "residual_sugar")}
    dry = p + sum(imp.values())
    moisture_frac = 0.052
    total = dry / (1.0 - moisture_frac)
    out = {"purity_dry_pct": p / dry * 100.0,
           "moisture_pct": moisture_frac * 100.0,
           "total_powder_kg": total, "product_kg": p, "dry_kg": dry}
    for k, v in imp.items():
        out[f"{k}_pct"] = v / dry * 100.0
    return out


# ============================================================================
# STAGE 6 : 수확 및 균체 분리
# ============================================================================
def stage6_harvest(st: ProcessState) -> StageResult:
    before = st.mass_kg
    biomass_in = st.impurities.get("biomass", 0.0)
    rows, cost, hours = _run_unitops(st, ["kill", "dilute", "separate"])
    st.broth_L *= 2.5     # 희석 1:1.5

    b = []
    b.append(R.section("1) 수확 직전 배양액 상태"))
    pr = K.productivity(st.ferm)
    b.append(R.kv([("배양액", f"{pr['broth_volume_L']:,.0f} L, 겉보기 점도 {pr['peak_viscosity_cP']:.0f} cP (비뉴턴)"),
                   ("조성", f"목적 EPS {before:.1f} kg / 균체 {biomass_in:.0f} kg DCW / 잔당 {pr['residual_S_g_L']:.1f} g/L"),
                   ("난이도", "다당이 녹아있는 '묽은 젤리'. 펌핑·원심·여과 전부 이 점도와 싸운다")]))

    b.append(R.section("2) 단위조작별 물질수지"))
    b.append(R.table(["단위조작", "투입 kg", "회수율", "산출 kg", "손실 kg", "시간 h", "소모품"],
                     rows, "<>>>>>>", indent=4))

    b.append(R.section("3) 왜 균체 분리에서 가장 많이 잃는가"))
    b.append(R.bullets([
        "균사체 케이크는 스펀지처럼 EPS 용액을 머금는다 (held-up liquor). 케이크 1 kg 당 2~4 kg 액 보유",
        "대응 1: 디스크스택 원심분리 후 케이크 재슬러리 -> 2차 원심 (washing) 으로 5~7%p 회수",
        "대응 2: 희석배수를 늘리면 회수율은 오르지만 후단 UF 부하(처리 용적)가 비례해 증가 -> 최적점 존재",
        "대응 3: 규조토 프리코트 여과는 맑기는 좋으나 다당 흡착 손실 + 폐기물(규조토) 비용 발생",
        "로터리 진공 드럼 필터는 고점도에 취약. 이 공정에선 원심 -> 심층여과 조합이 현실적",
    ]))
    b.append(R.table(["분리 방식", "회수율", "청징도", "비고"],
                     [["디스크스택 원심 (단독)", "86~90%", "중", "연속운전 가능, 처리량 큼"],
                      ["원심 + 케이크 세척 1회", "91~94%", "중", "본 공정 채택. 용적 20% 증가 감수"],
                      ["원심 + 규조토 여과", "89~92%", "상", "최종 외관(투명도) 중요 시"],
                      ["막여과 (MF 0.2 um)", "88~93%", "상", "고점도에서 flux 급락, 막 세정비 큼"],
                      ["압착여과 (필터프레스)", "80~86%", "중상", "배치형, 인력 소요 큼"]], indent=4))

    b.append(R.section("4) 부산물 - 버리지 말 것"))
    cake_dry = pr['biomass_g_L'] * pr['broth_volume_L'] / 1000.0
    b.append(R.kv([("균사체 케이크", f"건조 약 {cake_dry:.0f} kg (습중량 기준 3~4 배)"),
                   ("활용 1", "세포벽 유래 베타글루칸 (불용성) -> 2차 알칼리 추출로 별도 원료화"),
                   ("활용 2", "키틴/키토산, 에르고스테롤(프로비타민 D2) 회수"),
                   ("활용 3", "사료/퇴비 (가장 흔함, 부가가치 최저)"),
                   ("미처리 시", f"폐기물 처리비 ${pr['broth_volume_L']/1000*C.DEFAULT_COST.waste_usd_per_m3*2.5:,.0f}/batch 발생")]))

    st.timeline_h["stage6"] = hours
    return StageResult(6, "수확 · 균체 분리",
                       "목적물은 상등액에, 균체는 케이크로 - 다운스트림 최대 손실 구간을 관리",
                       "\n".join(b), duration_h=hours, cost_usd=cost,
                       mass_in_kg=before, mass_out_kg=st.mass_kg,
                       cpp=[("실활 온도/시간", "80 degC / 30 min", "과열 시 다당 해중합 -> 분자량 규격 이탈"),
                            ("희석배수", "1 : 1.5 (정제수)", "부족 시 분리 불량, 과다 시 UF 부하 폭증"),
                            ("원심 G / 체류시간", "8,000 x g, 체류 90 s", "부족하면 탁도, 과하면 케이크 압밀로 액 회수 불가"),
                            ("케이크 세척", "1 회, 세척수 = 케이크의 1.5 배", "회수율 5~7%p 좌우")],
                       risks=["실활을 건너뛰면 균체 자가분해효소(glucanase)가 목적 다당을 분해한다",
                              "고점도 상태로 원심에 넣으면 디스크 간극 폐색 -> 장비 정지"],
                       gate=st.mass_kg / before > 0.85,
                       gate_text=f"구간 회수율 {st.mass_kg/before*100:.1f}%")


# ============================================================================
# STAGE 7 : 정제
# ============================================================================
def stage7_purification(st: ProcessState) -> StageResult:
    before = st.mass_kg
    imp0 = dict(st.impurities)
    rows, cost, hours = _run_unitops(st, ["enzyme", "uf", "carbon", "polish"])

    b = []
    b.append(R.section("1) 단위조작별 물질수지"))
    b.append(R.table(["단위조작", "투입 kg", "회수율", "산출 kg", "손실 kg", "시간 h", "소모품"],
                     rows, "<>>>>>>", indent=4))

    b.append(R.section("2) 불순물 제거 추적 (목적물 대비 %)"))
    irows = []
    for k, label in [("hmw_nonglucan", "공생성 고분자(UF 불가)"), ("protein", "잔류 단백질"), ("pigment", "색소(갈변물)"),
                     ("salt", "무기염"), ("residual_sugar", "저분자 당"),
                     ("endotoxin", "내독소(상대값)")]:
        a, c = imp0.get(k, 0.0), st.impurities.get(k, 0.0)
        base = max(before, 1e-9)
        irows.append([label, f"{a/base*100:.2f}%" if k != "endotoxin" else f"{a:.3f}",
                      f"{c/max(st.mass_kg,1e-9)*100:.3f}%" if k != "endotoxin" else f"{c:.3f}",
                      f"{(1-c/max(a,1e-12))*100:.1f}%"])
    b.append(R.table(["불순물", "정제 전", "정제 후", "총 제거율"], irows, "<>>>", indent=4))

    b.append(R.section("3) 정제 전략 - 왜 이 순서인가"))
    b.append(R.table(["순서", "조작", "노리는 것", "이 순서가 아니면"],
                     [["1", "효소(protease)", "단백을 '작게' 만들어 UF로 빠지게 함",
                       "큰 단백은 다당과 함께 농축되어 영원히 안 빠진다"],
                      ["2", "UF/DF 100 kDa", "저분자 전부 제거 + 농축 (정제의 주력)",
                       "먼저 농축하면 점도 폭증으로 flux 가 죽는다"],
                      ["3", "활성탄", "색소·이취·내독소 마무리",
                       "UF 전에 쓰면 탄소가 고형분에 막혀 효율 1/3"],
                      ["4", "0.2 um + 에탄올 침전", "무균화 + 최종 분획 회수",
                       "활성탄 미제거 상태로 침전하면 탄소가 제품에 혼입"]], indent=4))

    b.append(R.section("4) UF 운전 상세 (정제 성패의 80%)"))
    b.append(R.kv([("막", "PES 중공사, MWCO 100 kDa, 면적 40 m2 x 2 스키드"),
                   ("모드", "TFF (접선흐름). 농축 5배 -> 투석여과 8 DV -> 최종 농축"),
                   ("TMP", "0.8 ~ 1.2 bar (고분자 다당은 저 TMP 가 유리)"),
                   ("크로스플로우", "3.5 m/s (농도분극층 제거)"),
                   ("flux", "초기 28 LMH -> 종료 9 LMH (점도·겔층 형성)"),
                   ("MWCO 선택", "300 kDa 쓰면 수율 +6%p 이나 저분자 다당이 빠져 평균 분자량 상승"),
                   ("세정 CIP", "0.5 N NaOH 50 degC 60 min -> 회복률 >= 92% 확인")]))
    b.append(R.bullets([
        "MWCO 를 키우면 수율은 오르고 순도·분자량 분포 폭은 나빠진다 - '분획'이 곧 제품 스펙이다",
        "겔층(gel layer)이 형성되면 TMP 를 올려도 flux 가 안 오른다. 오히려 압밀로 더 떨어진다",
        "에탄올 침전은 회수증류로 90% 이상 재사용. 안 하면 원가와 폐기물이 동시에 터진다",
    ]))

    b.append(R.section("5) 유기용매를 쓸 것인가 (화장품 원료의 분기점)"))
    b.append(R.table(["경로", "순도", "수율", "규제/마케팅"],
                     [["UF 단독 (용매 free)", "82~88%", "높음", "클린뷰티 소구 용이, 잔류용매 시험 불필요"],
                      ["UF + 에탄올 침전", "88~93%", "중", "본 공정. 에탄올은 화장품 허용, 잔류 <1,000 ppm 관리"],
                      ["Sevag(클로로포름) 탈단백", "93~96%", "중", "잔류 염소계 용매 - 화장품 원료로 사실상 불가"],
                      ["이온교환 + 소수성 크로마토", "95%+", "낮음", "원가 3~5배. 의약품급에서만 정당화"]], indent=4))

    st.timeline_h["stage7"] = hours
    return StageResult(7, "정제 (효소 · 한외여과 · 흡착 · 침전)",
                       "규격 순도까지 끌어올리되, 수율과 분자량 분포를 함께 지켜낸다",
                       "\n".join(b), duration_h=hours, cost_usd=cost,
                       mass_in_kg=before, mass_out_kg=st.mass_kg,
                       cpp=[("protease 조건", "50 degC, pH 7.0, 4 h, 0.1% E/S", "과도 시 다당 측쇄 손상"),
                            ("UF MWCO", "100 kDa", "크면 수율↑ 순도↓, 작으면 flux 급락"),
                            ("투석여과량", ">= 8 DV", "미달 시 저분자 당·염 잔류 -> 순도 규격 미달"),
                            ("TMP", "<= 1.2 bar", "초과 시 겔층 압밀로 비가역 flux 손실"),
                            ("활성탄량", "0.5% w/v 이하", "과량 시 다당 흡착으로 수율 급락"),
                            ("에탄올 비율/온도", "3 vol, 4 degC, 12 h", "부족 시 미침전 손실, 과다 시 염 공침")],
                       risks=["UF 막 오염 누적으로 배치마다 flux 가 떨어지면 처리시간이 늘어나 자가분해·미생물 증식 위험",
                              "활성탄 미세분이 0.2 um 를 통과해 제품 흑점(black spec) 클레임으로 이어진 사례 다수"],
                       gate=st.mass_kg / before > 0.6,
                       gate_text=f"구간 회수율 {st.mass_kg/before*100:.1f}%")


# ============================================================================
# STAGE 8 : 농축 · 건조
# ============================================================================
def stage8_drying(st: ProcessState) -> StageResult:
    before = st.mass_kg
    # 정제 후 pH 조정(구연산/구연산나트륨)과 다당 대이온으로 들어오는 염.
    # 정제로 뺀 배지 유래 염과 달리 이건 '다시 넣는' 것이라 회분의 실제 주성분이 된다.
    st.impurities["salt"] = st.impurities.get("salt", 0.0) + st.rng.uniform(0.020, 0.032) * before
    rows, cost, hours = _run_unitops(st, ["dry"])
    comp = _solid_composition(st)
    st.powder_kg = comp["total_powder_kg"]

    b = []
    b.append(R.section("1) 건조 물질수지"))
    b.append(R.table(["단위조작", "투입 kg", "회수율", "산출 kg", "손실 kg", "시간 h", "소모품"],
                     rows, "<>>>>>>", indent=4))

    b.append(R.section("2) 건조 방식 비교 (분자량과 용해성이 갈린다)"))
    b.append(R.table(["방식", "수율", "kg당 비용", "분자량 영향", "용해성", "적용"],
                     [["분무건조 (inlet 170 / outlet 85 degC)", "92~94%", "기준 1.0",
                       "고온 단시간 - 5~10% 해중합", "양호(중공입자)", "본 공정 채택"],
                      ["동결건조", "96~98%", "4~6 배", "거의 없음", "매우 우수",
                       "고분자량 프리미엄 그레이드"],
                      ["진공 드럼건조", "90~93%", "0.8 배", "열이력 길어 갈변 위험", "보통", "저가 그레이드"],
                      ["건조 생략 (액상 1% 수용액)", "99%", "0.2 배", "없음",
                       "-", "국내 근거리 납품, 보존제 필수"]], indent=4))
    b.append(R.bullets([
        "분무건조 outlet 온도를 95 degC 이상 올리면 수분은 잘 빠지지만 갈변 + 용해지연(fish-eye) 발생",
        "다당 분말은 흡습성이 커서 알루미늄 라미네이트 + 제습포장 필수. 포장 실패가 미생물 규격 이탈로 직결",
    ]))

    b.append(R.section("3) 최종 분말 조성"))
    b.append(R.table(["성분", "함량", "규격"],
                     [["베타글루칸(목적물, 건조물 기준)", f"{comp['purity_dry_pct']:.1f} %", ">= 80%"],
                      ["공생성 헤테로다당·당단백 복합체", f"{comp['hmw_nonglucan_pct']:.2f} %", "(순도의 실질 상한 요인)"],
                      ["잔류 단백질", f"{comp['protein_pct']:.2f} %", "<= 2%"],
                      ["무기염 (pH조정·대이온 유래, 회분)", f"{comp['salt_pct']:.2f} %", "<= 6%"],
                      ["저분자 당", f"{comp['residual_sugar_pct']:.3f} %", "-"],
                      ["색소/기타", f"{comp['pigment_pct']:.3f} %", "-"],
                      ["수분", f"{comp['moisture_pct']:.1f} %", "<= 7%"],
                      ["합계 (분말 중량)", f"{comp['total_powder_kg']:.2f} kg", "-"]], "<>>", indent=4))

    st.timeline_h["stage8"] = hours
    return StageResult(8, "농축 · 건조 · 분말화",
                       "취급 가능한 고체 원료로 만들면서 분자량과 용해성을 지킨다",
                       "\n".join(b), duration_h=hours, cost_usd=cost,
                       mass_in_kg=before, mass_out_kg=st.mass_kg,
                       cpp=[("inlet / outlet 온도", "170 / 85 +- 3 degC", "outlet 과승 시 갈변·해중합·용해지연"),
                            ("공급 고형분", "4~6% w/w", "높으면 노즐 폐색, 낮으면 생산성·수율 저하"),
                            ("최종 수분", "<= 7%", "초과 시 보관 중 미생물 증식·케이킹"),
                            ("포장", "AL 라미네이트 + 제습제, 질소치환", "흡습 시 전 로트 규격 이탈")],
                       risks=["사이클론 미세분 손실이 수율의 3~5%. 백필터 추가로 회수 가능",
                              "분체 취급 시 분진폭발 하한(LEL) 관리 - 접지·방폭 필수"],
                       gate=comp["purity_dry_pct"] >= 80.0,
                       gate_text=f"순도 {comp['purity_dry_pct']:.1f}%, 분말 {comp['total_powder_kg']:.1f} kg")


# ============================================================================
# STAGE 9 : 표준화 · 제형화 (원료 그레이드 확정)
# ============================================================================
def stage9_formulation(st: ProcessState) -> StageResult:
    comp = _solid_composition(st)
    rng = st.rng
    powder = comp["total_powder_kg"]

    # 그레이드 배분
    grades = [
        ("A. 분말 원료 (Powder, >=85% glucan)", 0.35, 1.000, 420.0,
         "Schizophyllan", "고농도 에센스·앰플, 바이어 자체 희석"),
        ("B. 1% 수용액 (Liquid 1%)", 0.50, 0.010, 7.5,
         "Water / Schizophyllan / 1,2-Hexanediol", "범용. 대부분의 완제품사가 선호"),
        ("C. 0.5% + 글리세린 베이스", 0.15, 0.005, 6.2,
         "Water / Glycerin / Schizophyllan / 1,2-Hexanediol",
         "점도·사용감 튜닝된 프리미엄 블렌드"),
    ]
    grows = []
    total_rev = 0.0
    for name, share, conc, price, inci, use in grades:
        active = powder * share
        product_kg = active / conc
        rev = product_kg * price
        total_rev += rev
        grows.append([name, f"{share*100:.0f}%", f"{active:.1f}", f"{product_kg:,.0f}",
                      f"${price:,.2f}", f"${rev:,.0f}"])
    st._grades = grades
    st._revenue = total_rev

    b = []
    b.append(R.section("1) 표준화 - '순도'가 아니라 '역가'를 맞춘다"))
    b.append(R.bullets([
        f"벌크 분말 실측 베타글루칸 {comp['purity_dry_pct']:.1f}% -> 라벨 규격 '85.0% min' 으로 표준화",
        "배치 간 편차는 희석(정제수·말토덱스트린)으로 흡수하지 않고, 복수 배치 블렌딩으로 맞춘다",
        "  -> 희석제를 넣으면 INCI 에 그 성분을 표기해야 하고 '순수 원료' 소구가 깨진다",
        "표준화 기준은 반드시 '건조물 기준(dry basis)'. 수분 변동이 역가로 둔갑하는 것을 막는다",
    ]))

    b.append(R.section("2) 그레이드 구성 및 배치 산출"))
    b.append(R.table(["그레이드", "배분", "투입 벌크 kg", "제품 산출 kg", "단가/kg", "배치 매출"],
                     grows, "<>>>>>", indent=4))
    b.append(R.kv([("배치 총 매출(추정)", f"${total_rev:,.0f}"),
                   ("액상 그레이드 이유", "고분자 다당은 분말 재용해 시 fish-eye 발생 - 완제품사가 싫어한다"),
                   ("분말 그레이드 이유", "운송비/보관성. 해외 수출은 분말이 압도적으로 유리")]))

    b.append(R.section("3) 액상 그레이드 제형 설계"))
    b.append(R.table(["성분", "함량", "역할", "주의"],
                     [["정제수 (Aqua)", "to 100", "용매", "주사용수급 불필요, 화장품용 정제수(전도도 <1.3 uS/cm)"],
                      ["Schizophyllan", "1.00%", "유효성분", "건조물 환산 투입. 수화 2 h, 고전단 금지(해중합)"],
                      ["1,2-Hexanediol", "0.80%", "보존", "파라벤/페녹시에탄올 회피형. 단독으로는 약함"],
                      ["Ethylhexylglycerin", "0.10%", "보존 부스터", "1,2-HD 와 병용 시 상승효과"],
                      ["Citric acid / Na-citrate", "q.s.", "pH 5.5 +-0.3", "다당 안정 pH 구간"],
                      ["(선택) Phytic acid", "0.02%", "금속이온 봉쇄", "미량 Fe 에 의한 다당 산화절단 억제"]], indent=4))
    b.append(R.bullets([
        "보존력시험(ISO 11930 / USP 51 challenge test) 통과가 액상 그레이드 출하의 전제",
        "고분자 다당 용액은 그 자체가 미생물의 먹이. 보존 설계를 가볍게 보면 반드시 사고가 난다",
        "용해 공정에서 고전단 호모게나이저를 쓰면 분자량이 20~30% 깎인다 - 앵커 믹서 저속 수화가 정답",
    ]))

    b.append(R.section("4) 물성 스펙 (완제품사가 실제로 보는 숫자)"))
    b.append(R.table(["항목", "1% 수용액 기준", "의미"],
                     [["외관", "무색~미황색 투명~약간 유백 점성액", "탁하면 바로 반품"],
                      ["점도 (25 degC, LV#3 30 rpm)", "1,200 ~ 2,600 cP", "분자량의 간접 지표"],
                      ["pH", "5.0 ~ 6.0", "피부 친화 + 다당 안정"],
                      ["굴절률", "1.333 ~ 1.340", "고형분 간이 확인"],
                      ["내염성", "NaCl 1% 첨가 시 점도 변화 <= 20%", "제형 호환성"],
                      ["전단담화 지수 n", "0.3 ~ 0.5", "도포 시 가볍게 퍼지는 사용감의 근거"]], indent=4))

    st.timeline_h["stage9"] = 24
    return StageResult(9, "표준화 · 제형화 (원료 그레이드 확정)",
                       "같은 벌크에서 분말·액상 그레이드를 나누고, 라벨 역가를 배치 간 일정하게 고정",
                       "\n".join(b), duration_h=24, cost_usd=powder * 0.35 * 2.0 + 850,
                       cpp=[("표준화 역가", "85.0% min (dry basis)", "라벨 미달 시 전량 재작업"),
                            ("용해 전단", "앵커 50 rpm, 호모 금지", "고전단은 분자량 20~30% 손실"),
                            ("최종 pH", "5.5 +- 0.3", "pH 4 이하 산가수분해, 8 이상 알칼리 절단"),
                            ("보존 시스템", "1,2-HD 0.8% + EHG 0.1%", "challenge test 미통과 시 출하 불가")],
                       risks=["희석제를 써서 역가를 맞추면 INCI 변경 -> 바이어 재등록 -> 납기 지연",
                              "액상은 유통기한이 짧다(12~18개월). 분말(24~36개월) 대비 재고 리스크 큼"])


# ============================================================================
# STAGE 10 : 품질관리 규격 시험 (CoA)
# ============================================================================
def stage10_qc(st: ProcessState) -> StageResult:
    rng = st.rng
    comp = _solid_composition(st)
    Rr = st.ferm
    sp = st.strain

    # --- 분자량 모델: 생합성 Mw 에서 공정 이력만큼 해중합 ---
    f = C.EQUIPMENT["main"]
    Mw0 = 2600.0                                   # 생합성 직후 [kDa]
    tip = K.tip_speed_m_s(f, max(Rr.rpm))
    shear_h = sum(1 for r in Rr.rpm if r > f.max_rpm * 0.9) * 0.25
    deg_shear = math.exp(-0.00042 * tip ** 2 * shear_h)     # 기계적 전단
    deg_heat = 0.95                                          # 80 degC 실활
    deg_dry = 0.92                                           # 분무건조
    deg_lot = rng.uniform(0.95, 1.03)
    Mw = Mw0 * deg_shear * deg_heat * deg_dry * deg_lot
    PDI = rng.uniform(1.6, 2.4)

    measured = {
        "베타글루칸 순도 (건조물 기준)": comp["purity_dry_pct"] + rng.gauss(0, 0.8),
        "평균 분자량 Mw": Mw,
        "잔류 단백질": max(comp["protein_pct"] + rng.gauss(0, 0.05), 0.01),
        "건조감량 (수분)": comp["moisture_pct"] + rng.gauss(0, 0.4),
        "강열잔분 (회분)": comp["salt_pct"] * 1.9 + rng.gauss(0.6, 0.25),
        "pH (1% 수용액)": 5.5 + rng.gauss(0, 0.25),
        "잔류 에탄올": max(rng.gauss(180, 90), 5),
        "총호기성생균수": max(rng.gauss(12, 14), 0),
        "대장균/녹농균/황색포도상구균": 0.0,
        "내독소 (주사제급 요구 시)": max(rng.gauss(18, 9), 1),
        "납 Pb": max(rng.gauss(1.1, 0.6), 0.05),
        "비소 As": max(rng.gauss(0.5, 0.3), 0.02),
        "수은 Hg": max(rng.gauss(0.03, 0.02), 0.001),
        "카드뮴 Cd": max(rng.gauss(0.09, 0.05), 0.005),
    }

    rows, fails = [], []
    for s_ in C.SPECS:
        v = measured[s_.item]
        ok = s_.lo <= v <= s_.hi
        if not ok and s_.critical:
            fails.append(s_.item)
        limit = (f"{s_.lo:g} ~ {s_.hi:g}" if s_.lo > 0 else f"<= {s_.hi:g}")
        rows.append([s_.item, s_.method, f"{limit} {s_.unit}",
                     f"{v:,.2f}" if v >= 1 else f"{v:.3f}",
                     "적합" if ok else ("부적합" if s_.critical else "참고"),
                     "" if s_.critical else "(참고항목)"])
    st.qc = rows
    st.released = len(fails) == 0

    b = []
    b.append(R.section("1) 시험성적서 (CoA) - 배치 BG-2026-0142"))
    b.append(R.table(["시험항목", "시험방법", "규격", "실측", "판정", "비고"], rows, indent=2))

    b.append(R.section("2) 분자량 - 이 제품의 '효능'을 실제로 결정하는 값"))
    b.append(R.table(["단계", "Mw 영향", "환산 Mw (kDa)"],
                     [["생합성 직후 (세포 외 분비)", "기준", f"{Mw0:,.0f}"],
                      [f"발효 중 전단 (tip {tip:.1f} m/s x {shear_h:.0f} h)",
                       f"x {deg_shear:.3f}", f"{Mw0*deg_shear:,.0f}"],
                      ["가열 실활 80 degC/30 min", f"x {deg_heat:.2f}", f"{Mw0*deg_shear*deg_heat:,.0f}"],
                      ["분무건조", f"x {deg_dry:.2f}", f"{Mw0*deg_shear*deg_heat*deg_dry:,.0f}"],
                      ["로트 편차", f"x {deg_lot:.2f}", f"{Mw:,.0f}"]], "<>>", indent=4))
    b.append(R.kv([("최종 Mw", f"{Mw:,.0f} kDa (규격 {sp.mw_target_kda[0]:,.0f} ~ {sp.mw_target_kda[1]:,.0f})"),
                   ("분산도 PDI", f"{PDI:.2f}  (Mw/Mn, 2.5 초과 시 분획 불량 의심)"),
                   ("효능 상관", "보습·필름형성은 고분자(>1,000 kDa)가 유리 / 경피흡수·면역조절은 저분자가 유리"),
                   ("시사점", "용도를 먼저 정하고 그 Mw 를 공정으로 '만드는' 것이지, 나온 걸 파는 게 아니다")]))

    b.append(R.section("3) 시험법 주의 - 다당 정량은 틀리기 쉽다"))
    b.append(R.table(["방법", "장점", "함정"],
                     [["Phenol-H2SO4 (총당)", "간편·저가", "베타글루칸이 아닌 모든 당을 센다 -> 과대평가"],
                      ["효소법 (Megazyme beta-glucan)", "특이성 높음", "불용성 분획 전처리 미흡 시 과소평가"],
                      ["Congo red 결합", "삼중나선 구조 확인", "정량성 낮음, 보조 지표로만"],
                      ["NMR (1H/13C)", "결합양식(1->3/1->6) 확정", "고가·시간. 원료 등록 시 1회는 필수"],
                      ["GPC-MALS", "절대 분자량", "표준물질(pullulan) 환산은 상대값 - MALS 병행 필요"]], indent=4))

    b.append(R.section("4) 판정"))
    if st.released:
        b.append("    모든 필수 규격 적합 -> " + R.verdict(True, "출하 승인 (Released)"))
    else:
        b.append("    " + R.verdict(False, "", f"부적합: {', '.join(fails)}"))
        b.append(R.bullets([
            "OOS(Out Of Specification) 조사 절차 가동: 시험오류 배제 -> 공정 원인 조사 -> CAPA",
            "재작업(rework) 가능 항목인지 판단 (예: 수분 초과 -> 재건조 / 순도 미달 -> 재UF)",
            "재작업 불가 시 등급강등(저농도 그레이드) 또는 폐기",
        ]))

    cost = C.DEFAULT_COST.qc_usd_per_batch
    st.add_cost("QC 시험", cost)
    st.timeline_h["stage10"] = 120
    return StageResult(10, "품질관리 규격 시험 (CoA 발행)",
                       "라벨에 쓴 숫자를 증명한다 - 이게 없으면 원료가 아니라 그냥 가루다",
                       "\n".join(b), duration_h=120, cost_usd=cost,
                       cpp=[("샘플링", "상·중·하 3점 혼합", "분체는 층분리로 단일점 샘플이 대표성 없음"),
                            ("표준품", "인증 pullulan + beta-glucan 표준", "표준품 없으면 수치의 의미 없음"),
                            ("미생물 시험", "대한민국화장품안전기준 통칙", "중화제 미사용 시 보존제가 균을 죽여 위음성"),
                            ("보관검체", "배치당 3배량, 유효기간 +1 년", "클레임 대응의 유일한 증거")],
                       risks=["Phenol-H2SO4 단독 정량으로 '순도 95%' 를 광고하다 효소법 재검에서 70%대가 나와 분쟁이 된 사례 다수",
                              "내독소는 화장품 필수항목은 아니나, 더마·시술후 제품 바이어는 요구한다"],
                       gate=st.released, gate_text="출하 승인" if st.released else "OOS 조사")


# ============================================================================
# STAGE 11 : 안전성 · 효능 평가
# ============================================================================
def stage11_safety_efficacy(st: ProcessState) -> StageResult:
    rng = st.rng
    b = []
    b.append(R.section("1) 안전성 - 동물실험 금지(EU/국내) 하에서의 표준 패키지"))
    b.append(R.table(["엔드포인트", "시험법 (OECD TG)", "결과", "판정"],
                     [["피부자극", "TG 439 재구성인체표피 (RhE, KeraSkin/EpiSkin)",
                       f"세포생존율 {rng.uniform(88,97):.1f}% (기준 >50%)", "비자극"],
                      ["안자극", "TG 492 RhCE (EpiOcular)",
                       f"생존율 {rng.uniform(82,95):.1f}% (기준 >60%)", "비자극"],
                      ["피부감작", "TG 442C DPRA + 442D KeratinoSens + 442E h-CLAT (2 of 3)",
                       "3종 모두 음성", "비감작"],
                      ["유전독성", "TG 471 Ames (5 균주, +-S9)", "복귀변이 음성", "음성"],
                      ["광독성", "TG 432 3T3 NRU", f"PIF {rng.uniform(1.0,1.6):.2f} (<2)", "광독성 없음"],
                      ["세포독성", "MTT (HaCaT, HDF)", f"IC50 > 1,000 ug/mL", "저독성"],
                      ["인체첩포", "HRIPT 50~100 명 / 또는 1차 첩포 33 명",
                       "자극지수 0.00~0.12", "저자극"]], indent=4))
    b.append(R.bullets([
        "미생물 발효 원료는 '내독소'와 '잔류 균체단백'이 자극의 주원인. 정제 수준이 곧 안전성이다",
        "원료 자체 시험 외에 '최종 제형에서의 사용농도' 기준 시험 자료를 같이 준비해야 바이어가 받는다",
        "CIR / SCCS 기존 평가자료가 있는 성분이면 이 비용의 절반 이상을 아낄 수 있다",
    ]))

    b.append(R.section("2) 효능 - in vitro"))
    b.append(R.table(["소구점", "모델", "지표", "결과 (0.1% 처리)"],
                     [["보습", "HaCaT", "AQP3 mRNA 발현", f"대조 대비 {rng.uniform(1.6,2.3):.2f} 배"],
                      ["보습", "HaCaT", "히알루론산 합성(HAS2)", f"{rng.uniform(1.4,2.0):.2f} 배"],
                      ["장벽", "RhE", "TEWL / 필라그린·로리크린", f"필라그린 {rng.uniform(1.3,1.9):.2f} 배"],
                      ["진정/항염", "LPS 자극 RAW264.7", "NO, IL-6, PGE2 억제", f"IL-6 {rng.uniform(38,62):.0f}% 억제"],
                      ["항산화", "DPPH / 세포내 ROS", "ROS 소거", f"{rng.uniform(30,52):.0f}% 감소"],
                      ["재생", "HDF scratch assay", "이동률", f"{rng.uniform(120,165):.0f}% (대조 100)"],
                      ["필름형성", "In vitro film", "수분증발 억제", f"{rng.uniform(22,40):.0f}% 억제"]], indent=4))

    b.append(R.section("3) 효능 - 인체적용시험 (마케팅 문구의 근거)"))
    b.append(R.table(["항목", "설계", "결과"],
                     [["피험자", "22~33 명, 4 주, 1일 2회 도포", "-"],
                      ["피부수분량", "Corneometer", f"4 주 후 +{rng.uniform(18,34):.1f}% (p<0.05)"],
                      ["경피수분손실 TEWL", "Tewameter", f"-{rng.uniform(12,24):.1f}% (p<0.05)"],
                      ["피부거칠기", "Visioscan/Antera 3D", f"-{rng.uniform(8,18):.1f}%"],
                      ["홍반지수", "Mexameter", f"-{rng.uniform(6,15):.1f}%"],
                      ["이상반응", "피부과 전문의 판독", "없음"]], indent=4))
    b.append(R.bullets([
        "'보습'은 일반 화장품 문구로 가능. '미백·주름개선'을 쓰려면 국내 기능성화장품 심사가 별도로 필요",
        "기능성 고시원료가 아니면 개별 심사(자료제출) - 원료사가 아니라 완제품사가 받는 것이 일반적",
        "인체적용시험은 IRB 승인 + 시험기관(비임상/임상) 계약이 필요. 최소 3~4 개월",
    ]))

    cost = 68000.0
    b.append(R.section("4) 비용·기간 (여기가 신규 원료 사업의 진짜 진입장벽)"))
    b.append(R.table(["패키지", "비용(USD)", "기간"],
                     [["OECD in vitro 안전성 5종", "28,000 ~ 45,000", "3 ~ 4 개월"],
                      ["HRIPT 50 명", "6,000 ~ 12,000", "2 개월"],
                      ["in vitro 효능 5종", "12,000 ~ 20,000", "2 ~ 3 개월"],
                      ["인체적용시험 1건", "15,000 ~ 30,000", "4 ~ 5 개월"],
                      ["합계 (1회성, 배치 원가 아님)", "61,000 ~ 107,000", "병행 시 6 ~ 8 개월"]], "<>>", indent=4))

    st._dossier_cost = cost
    st.timeline_h["stage11"] = 180 * 24
    return StageResult(11, "안전성 · 효능 평가",
                       "원료를 '팔 수 있는 것'으로 만드는 서류 - 공정보다 돈과 시간이 더 드는 구간",
                       "\n".join(b), duration_h=180 * 24, cost_usd=0.0,
                       cpp=[("시험 검체", "출하 규격품 3 배치", "파일럿 검체로 받은 자료는 바이어가 거절"),
                            ("시험농도", "원료 100% + 사용농도", "원료만 하면 실사용 근거 부족"),
                            ("GLP 여부", "안전성은 GLP 기관 권장", "비GLP 자료는 EU 수출 시 거절될 수 있음")],
                       risks=["3 배치 자료를 요구받는 경우가 많다 - 공정 확립 전에 시험을 시작하면 전부 다시 해야 한다",
                              "효능 데이터는 '원료 농도' 기준인데 완제품에 0.01% 넣고 같은 문구를 쓰면 허위광고"],
                       gate=True, gate_text="안전성 패키지 확보")


# ============================================================================
# STAGE 12 : 규제 · 등록
# ============================================================================
def stage12_regulatory(st: ProcessState) -> StageResult:
    sp = st.strain
    b = []
    b.append(R.section("1) 명명 및 국제 등재"))
    b.append(R.table(["항목", "내용", "소요"],
                     [["INCI 명", sp.inci, "PCPC 신규 신청 시 6~12 개월"],
                      ["ICID 등재", "International Cosmetic Ingredient Dictionary", "연 2회 개정판 반영"],
                      ["CAS / EC No.", "기존 등재 확인, 없으면 신규 부여 신청", "3~6 개월"],
                      ["원료 표준화 자료", "제조방법·규격·시험법·구조분석(NMR/FT-IR/GPC)", "자사 작성"]], indent=4))

    b.append(R.section("2) 시장별 규제 (수출 계획을 먼저 정하고 역산해야 한다)"))
    b.append(R.table(["시장", "요구", "리드타임", "난이도"],
                     [["대한민국", "화장품법상 원료 네거티브 규제 - 사용금지·제한 목록 비해당이면 사용 가능.\n"
                       "   원료 성적서·안전성 자료는 완제품사 책임자료로 보관", "즉시", "낮음"],
                      ["EU", "CosIng 등재 + CPNP 신고는 완제품 단위. 원료사는 CPSR 작성용 자료 패키지 제공.\n"
                       "   나노물질·CMR 비해당 선언", "1~2 개월", "중"],
                      ["미국", "MoCRA(2022) 시행 - 시설등록·제품리스팅 의무. 원료는 FDA 사전승인 불요",
                       "1~2 개월", "중"],
                      ["중국", "NMPA - '已使用化妆品原料目录' 수재 여부가 전부.\n"
                       "   미수재 = 신원료 등록(备案/注册), 독성자료 대량 요구", "2~3 년", "매우 높음"],
                      ["일본", "화장품기준 - 배합금지·제한 비해당 확인. 생물유래는 추가 자료 요구 가능",
                       "1~3 개월", "중"],
                      ["ASEAN", "ACD(ASEAN Cosmetic Directive) Annex 확인", "1 개월", "낮음"]], indent=4))

    b.append(R.section("3) 생물유래 원료 특유의 서류"))
    b.append(R.table(["항목", "내용", "비고"],
                     [["ABS / 나고야의정서", "국내 분리균 -> 생물다양성법 국내 접근신고로 충분.\n"
                       "   해외 채집·분양균은 PIC/MAT 필수", "미비 시 판매 자체가 위법"],
                      ["LMO/GMO 해당 여부", "야생형 균주 - LMO 비해당 확인서", "유전자조작주면 LMO법 전면 적용"],
                      ["BSE/TSE free", "배지에 동물유래 성분 미사용 선언", "대두펩톤 사용으로 충족"],
                      ["Vegan / Halal / Cosmos", "원료·배지·공정보조제 전수 선언", "인증은 선택이나 바이어가 요구"],
                      ["알레르겐", "대두·밀 유래 잔류단백 불검출 자료", "ELISA 정량 권장"],
                      ["RSPO / 미생물오염", "소포제·계면활성제 유래 팜유 추적", "클린뷰티 바이어 요구 증가"]], indent=4))

    b.append(R.section("4) 안정성 시험 (유효기간의 근거)"))
    b.append(R.table(["조건", "기간", "측정", "판정 기준"],
                     [["장기 25 degC / 60% RH", "0,3,6,9,12,18,24,36 개월",
                       "순도·Mw·pH·외관·미생물", "순도 초기값의 95% 이상"],
                      ["가속 40 degC / 75% RH", "0,1,2,3,6 개월", "동일", "6개월 통과 = 상온 2년 추정"],
                      ["가혹 (광, 동결융해)", "ICH Q1B / FT 3 cycle", "외관·점도", "변화 없음"],
                      ["개봉 후 (액상)", "PAO 시험", "미생물", "12 M 설정"]], indent=4))
    b.append(R.kv([("분말 유효기간", "36 개월 (밀봉, 25 degC 이하, 차광)"),
                   ("액상 유효기간", "18 개월 (보존계 포함)"),
                   ("주의", "가속시험만으로 설정하고 장기시험을 안 돌리면, 24개월차 클레임 때 방어 불가")]))

    b.append(R.section("5) 지식재산"))
    b.append(R.bullets([
        "균주 자체: 기탁(KCTC/부다페스트조약) + 균주 특허 - 가장 방어력이 높다",
        "제조방법: 배지·유가식 프로파일·정제조건 특허. 다만 리버스엔지니어링이 어려워 영업비밀이 유리할 때도 많다",
        "용도: '특정 Mw 범위 다당의 피부장벽 개선 용도' - 효능 데이터와 묶어야 등록 가능성이 오른다",
        "상표: 원료 브랜드명은 INCI 와 별개. 바이어가 완제품에 표기할 수 있어 마케팅 자산이 된다",
    ]))

    st.timeline_h["stage12"] = 365 * 24
    return StageResult(12, "규제 · 등록 · 지식재산",
                       "기술이 아니라 서류로 시장이 열린다 - 수출 시장을 먼저 정해야 공정 설계가 끝난다",
                       "\n".join(b), duration_h=365 * 24, cost_usd=0.0,
                       cpp=[("INCI 명칭", "등재명과 CoA 표기 100% 일치", "불일치 시 완제품 라벨 전체 위반"),
                            ("ABS 증빙", "채집지·일자·신고번호 보관", "사후 소급 불가 - 채집 시점에 해야 함"),
                            ("안정성 시험 개시", "3 배치 동시 착수", "늦게 시작하면 영업이 1년 늦어진다")],
                       risks=["중국 비수재 원료는 사실상 2~3년간 중국 수출 불가. 사업계획에서 가장 흔한 오산",
                              "INCI 신규 등재 중에 샘플을 뿌리면, 등재명이 바뀌었을 때 모든 자료를 다시 돌려야 한다"],
                       gate=True, gate_text="국내·EU·미국 출시 가능 / 중국은 신원료 등록 필요")


# ============================================================================
# STAGE 13 : 전체 물질수지 · 원가 · 스케일업
# ============================================================================
def stage13_economics(st: ProcessState) -> StageResult:
    cm = C.DEFAULT_COST
    pr = K.productivity(st.ferm)
    comp = _solid_composition(st)
    main = C.EQUIPMENT["main"]
    ref = C.EQUIPMENT["seed2"]

    # --- 추가 원가 항목 ---
    dsp_h = sum(st.timeline_h.get(f"stage{i}", 0.0) for i in (6, 7, 8, 9))
    st.add_cost("다운스트림 인건", dsp_h * cm.labor_usd_per_h)
    waste_m3 = st.broth_L / 1000.0
    st.add_cost("폐수·폐기물 처리", waste_m3 * cm.waste_usd_per_m3)
    st.add_cost("포장재", comp["total_powder_kg"] * 1.8 + 420)
    st.add_cost("정제수(공정용수)", (st.broth_L * 2.2) / 1000.0 * cm.wfi_usd_per_m3)
    direct = sum(st.costs.values())
    st.add_cost("제조간접비", direct * cm.overhead_rate)
    dep = cm.capex_usd / cm.depreciation_years / cm.batches_per_year
    st.add_cost("설비 감가상각", dep)
    total = sum(st.costs.values())

    # --- 물질수지 트리 ---
    glu_total = C.MEDIUM_MAIN[0].g_per_L * main.working_volume_L / 1000.0 + pr["glucose_fed_kg"]
    in_broth = pr["product_in_broth_kg"]
    chain = [["발효 종료 (broth 내 EPS)", in_broth, 100.0]]
    for d in st.dsp_log:
        chain.append([d["op"], d["out"], d["out"] / in_broth * 100.0])

    b = []
    b.append(R.section("1) 전체 물질수지 (원료 투입 -> 최종 분말)"))
    b.append(R.table(["흐름", "질량", "비고"],
                     [["투입 글루코스 (회분+유가식)", f"{glu_total:,.0f} kg",
                       f"회분 {C.MEDIUM_MAIN[0].g_per_L*main.working_volume_L/1000:,.0f} kg + 유가식 {pr['glucose_fed_kg']:,.0f} kg"],
                      ["-> 균체 (부산물)", f"{pr['biomass_g_L']*pr['broth_volume_L']/1000:,.0f} kg DCW",
                       "2차 원료화 또는 사료/퇴비"],
                      ["-> CO2 + 유지대사 (손실)",
                       f"{max(glu_total - pr['biomass_g_L']*pr['broth_volume_L']/1000 - in_broth, 0):,.0f} kg",
                       "회수 불가. 이 비중이 곧 발효 수율의 한계"],
                      ["-> 목적 EPS (broth)", f"{in_broth:,.1f} kg", f"총괄 Yp/s {in_broth/glu_total:.3f} g/g"],
                      ["=> 최종 분말", f"{comp['total_powder_kg']:,.1f} kg",
                       f"순도 {comp['purity_dry_pct']:.1f}%, 수분 {comp['moisture_pct']:.1f}%"]], "<>>", indent=4))

    b.append(R.section("2) 다운스트림 수율 전파 (어디서 새는가)"))
    b.append(R.table(["단계", "잔존 질량 kg", "발효 대비 누적 %", "구간 손실 kg"],
                     [[c[0], f"{c[1]:.2f}", f"{c[2]:.1f}%",
                       f"{(chain[i-1][1]-c[1]):.2f}" if i else "-"]
                      for i, c in enumerate(chain)], "<>>>", indent=4))
    b.append(R.sparkline([c[2] for c in chain], width=50, height=7,
                         label="누적 회수율 (%) : 발효 -> 분말", lo=0, hi=100))
    dsp_yield = st.mass_kg / in_broth
    b.append(R.kv([("다운스트림 총 회수율", f"{dsp_yield*100:.1f} %"),
                   ("최대 손실 구간", "균체 분리 (케이크 보유액) > UF > 에탄올 침전"),
                   ("개선 1순위", "케이크 2차 세척 -> +5~7%p, 투자 대비 회수 가장 빠름")]))

    b.append(R.section("3) 배치 원가 (COGS) 분해"))
    crows = sorted(st.costs.items(), key=lambda x: -x[1])
    b.append(R.table(["항목", "금액", "비중", "kg당"],
                     [[k, f"${v:,.0f}", f"{v/total*100:.1f}%",
                       f"${v/comp['total_powder_kg']:,.2f}"] for k, v in crows]
                     + [["합계", f"${total:,.0f}", "100.0%",
                         f"${total/comp['total_powder_kg']:,.2f}"]], "<>>>", indent=4))
    cogs_kg = total / comp["total_powder_kg"]
    cogs_active = total / (comp["total_powder_kg"] * comp["purity_dry_pct"] / 100.0)

    b.append(R.section("4) 수익성"))
    rev = getattr(st, "_revenue", 0.0)
    b.append(R.table(["지표", "배치", "연간 (38 batch)"],
                     [["매출", f"${rev:,.0f}", f"${rev*cm.batches_per_year:,.0f}"],
                      ["제조원가", f"${total:,.0f}", f"${total*cm.batches_per_year:,.0f}"],
                      ["매출총이익", f"${rev-total:,.0f}", f"${(rev-total)*cm.batches_per_year:,.0f}"],
                      ["매출총이익률", f"{(rev-total)/max(rev,1)*100:.1f}%", "-"],
                      ["분말 생산량", f"{comp['total_powder_kg']:.1f} kg",
                       f"{comp['total_powder_kg']*cm.batches_per_year/1000:.2f} ton"]], "<>>", indent=4))
    b.append(R.kv([("분말 kg당 제조원가", f"${cogs_kg:,.2f}"),
                   ("유효성분 kg당 원가", f"${cogs_active:,.2f}"),
                   ("시장가 참고", "화장품용 베타글루칸 분말 $250~700/kg (순도·분자량·브랜드에 따라)"),
                   ("원가율 (COGS/매출)", f"{total/max(rev,1)*100:.0f}%"),
                   ("1회성 투자", f"설비 ${cm.capex_usd:,.0f} + 안전성/효능 자료 ${getattr(st,'_dossier_cost',0):,.0f} + 등록비")]))

    b.append(R.section("5) 민감도 (무엇을 1% 개선하면 가장 이득인가)"))
    sens = []
    for label, delta_desc, new_cogs in [
        ("역가 +20% (11.5 -> 13.8 g/L)", "발효 개선", cogs_kg / 1.20),
        ("다운스트림 회수율 +5%p", "케이크 2차 세척", cogs_kg * dsp_yield / (dsp_yield + 0.05)),
        ("배양시간 -24 h (144 -> 120)", "연 배치수 38 -> 44", cogs_kg * (1 - dep / total * 0.14)),
        ("글루코스 단가 -20%", "대량구매 계약", cogs_kg - (st.costs.get("유가식 당 공급", 0) + st.costs.get("배지 원재료", 0)) * 0.2 * 0.6 / comp["total_powder_kg"]),
        ("동결건조 전환", "수율 +4%p / 비용 5배", cogs_kg * 1.28),
        ("배치수 38 -> 50 (가동률)", "감가상각 분산", cogs_kg - dep * (1 - 38 / 50) / comp["total_powder_kg"]),
    ]:
        sens.append([label, delta_desc, f"${new_cogs:,.2f}",
                     f"{(new_cogs-cogs_kg)/cogs_kg*100:+.1f}%"])
    b.append(R.table(["시나리오", "수단", "kg당 원가", "변화"], sens, "<<>>", indent=4))
    b.append(R.bullets([
        "원가 1위는 거의 항상 '발효 역가'. 정제 수율 개선보다 역가 개선의 레버리지가 크다",
        "다만 역가를 올리면 점도가 올라 다운스트림이 나빠진다 - 전체 최적화는 역가 단독 최대화가 아니다",
    ]))

    b.append(R.section("6) 스케일업 기준 - 무엇을 고정할 것인가 (전부는 불가능하다)"))
    D_ref, D_main = K.impeller_diameter_m(ref), K.impeller_diameter_m(main)
    r = D_ref / D_main
    N_ref = ref.max_rpm
    rows = [
        ["일정 P/V (비동력)", "N2 = N1 (D1/D2)^(2/3)", f"{N_ref*r**(2/3):.0f} rpm",
         "산소전달 유지. 대부분의 호기발효 표준 기준"],
        ["일정 tip speed (전단)", "N2 = N1 (D1/D2)", f"{N_ref*r:.0f} rpm",
         "전단 민감 균사·고분자 다당 보호. 단 산소전달 부족"],
        ["일정 Reynolds 수", "N2 = N1 (D1/D2)^2", f"{N_ref*r**2:.0f} rpm",
         "혼합 상사. 실제로는 교반이 턱없이 부족"],
        ["일정 혼합시간", "N2 = N1 (D2/D1)^(2/3)", f"{N_ref/r**(2/3):.0f} rpm",
         "물리적으로 불가능 (동력 폭증)"],
    ]
    b.append(R.table(["기준", "식", f"{ref.working_volume_L:.0f} L({N_ref:.0f} rpm) -> 7,000 L", "의미"],
                     rows, indent=4))
    b.append(R.kv([("본 공정 채택", f"일정 P/V -> {N_ref*r**(2/3):.0f} rpm 계산, 실제 설정 {main.max_rpm:.0f} rpm"),
                   ("검증 지표", f"tip speed {K.tip_speed_m_s(main, main.max_rpm):.2f} m/s (<= 5 m/s 유지)"),
                   ("실제 운영", "scale-down model(2 L)을 만들어 DO 0% 구간을 재현한 뒤 대책을 검증한다"),
                   ("스케일업 함정", "소형조에서는 절대 안 생기는 '농도구배(기질·DO·pH)'가 대형조의 주 문제")]))
    b.append(R.table(["스케일업에서 반드시 깨지는 것", "원인", "대책"],
                     [["산소전달 kLa", "점도 + 낮은 P/V", "가압, 산소부화, 유가식 속도 하향"],
                      ["혼합시간", "용적^(1/3) 로 악화", "공급 분산 노즐 다점화"],
                      ["열전달", "표면적/부피 비 감소", "내부 코일, 외부 열교환기"],
                      ["CO2 축적", "액고(liquid height) 증가", "통기량 하한 설정, stripping"],
                      ["포말", "단백 농도 x 통기량", "기계식 포말제거기 + 소포제 최소화"],
                      ["멸균 열이력", "승온/냉각 시간 증가", "연속멸균(CST) 전환"]], indent=4))

    st.timeline_h["stage13"] = 0
    total_days = sum(st.timeline_h.get(f"stage{i}", 0) for i in (3, 4, 5, 6, 7, 8, 9, 10)) / 24.0
    b.append(R.section("7) 배치 사이클 타임"))
    b.append(R.table(["구간", "시간 (h)", "일"],
                     [[f"Stage {i}", f"{st.timeline_h.get(f'stage{i}',0):,.1f}",
                       f"{st.timeline_h.get(f'stage{i}',0)/24:.1f}"]
                      for i in (3, 4, 5, 6, 7, 8, 9, 10)]
                     + [["CIP/전환 + 버퍼", "48.0", "2.0"],
                        ["합계 (QC 제외 실공정)",
                         f"{total_days*24 - st.timeline_h.get('stage10',0) + 48:,.1f}",
                         f"{(total_days*24 - st.timeline_h.get('stage10',0) + 48)/24:.1f}"]], "<>>", indent=4))
    occupancy = (st.timeline_h.get("stage3", 0) + st.timeline_h.get("stage5", 0)
                 + st.timeline_h.get("stage6", 0) + 24.0)   # 멸균+발효+수확배출+CIP
    lead_time = sum(st.timeline_h.get(f"stage{i}", 0) for i in (3, 4, 5, 6, 7, 8, 9, 10)) + 48
    b.append(R.kv([
        ("배치 리드타임 (투입 -> 출하)", f"{lead_time:,.0f} h = {lead_time/24:.1f} 일  (종균 {st.timeline_h.get('stage4',0)/24:.1f} 일 포함)"),
        ("발효조 점유 시간 (= 배치 간격)", f"{occupancy:,.0f} h = {occupancy/24:.1f} 일"),
        ("이론 연 배치수", f"{365*24/occupancy:.0f} 배치 (종균·다운스트림은 병렬 운전 전제)"),
        ("계획 연 배치수", f"{cm.batches_per_year} 배치 (정기보수·QC 홀드·수요 변동 반영, 가동률 {cm.batches_per_year/(365*24/occupancy)*100:.0f}%)"),
        ("병목", "본 발효 144 h. 종균·정제는 병렬로 숨길 수 있으나 발효조는 못 숨긴다"),
        ("개선 수단", "발효조 2기 교대 운전 -> 연 배치수 2배, 감가상각 단가는 오히려 하락")]))

    st._cogs_kg = cogs_kg
    st._total_cost = total
    return StageResult(13, "전체 물질수지 · 원가 · 스케일업",
                       "공정 전체를 숫자 하나로 묶는다 - 무엇을 고치면 돈이 되는지 여기서만 보인다",
                       "\n".join(b), duration_h=0.0, cost_usd=0.0,
                       mass_out_kg=comp["total_powder_kg"],
                       cpp=[("배치 수율", f">= {comp['total_powder_kg']*0.9:.0f} kg", "미달 시 원가 구조 붕괴"),
                            ("가동률", ">= 38 batch/yr", "감가상각 분산의 핵심 변수"),
                            ("역가", ">= 10 g/L", "원가 민감도 1위 인자")],
                       risks=["설비 투자 회수 전에 바이어 단가 인하 압박이 들어온다 - 그레이드 다각화로 방어",
                              "단일 품목·단일 고객 구조는 가동률 리스크가 치명적"],
                       gate=True, gate_text=f"kg당 ${cogs_kg:,.0f}")


# ============================================================================
# 파이프라인
# ============================================================================
ALL_STAGES = [
    stage0_isolation, stage1_cellbank, stage2_medium, stage3_sterilization,
    stage4_seed_train, stage5_fermentation, stage6_harvest, stage7_purification,
    stage8_drying, stage9_formulation, stage10_qc, stage11_safety_efficacy,
    stage12_regulatory, stage13_economics,
]


def run_pipeline(strain_key: str = "schizophyllum", seed: int = C.DEFAULT_SEED,
                 only: Optional[List[int]] = None) -> Tuple[ProcessState, List[StageResult]]:
    st = ProcessState(strain_key=strain_key, rng=random.Random(seed))
    results = []
    for fn in ALL_STAGES:
        r = fn(st)
        results.append(r)
    if only is not None:
        results = [r for r in results if r.idx in only]
    return st, results
