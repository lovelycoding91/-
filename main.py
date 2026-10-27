#!/usr/bin/env python3
"""균류 기반 화장품 원료 제조 전공정 시뮬레이터.

사용 예:
    python3 main.py                       # 전체 0~13 단계
    python3 main.py --stages 5,6,7        # 특정 단계만
    python3 main.py --strain tremella     # 다른 균주 트랙
    python3 main.py --monte-carlo 40      # 배치 간 변동 시뮬레이션
    python3 main.py --compare             # 균주 3종 비교
    python3 main.py --out out/report.txt  # 파일로 저장
"""
import argparse
import sys
from io import StringIO

from fungisim import config as C
from fungisim import kinetics as K
from fungisim import montecarlo as MC
from fungisim import report as R
from fungisim import stages as S


def header(strain_key: str) -> str:
    sp = C.STRAINS[strain_key]
    out = ["", "#" * R.WIDTH,
           R.pad(" 균류 활용 화장품 원료 제조 - 전공정 시뮬레이션 ", R.WIDTH, "^"),
           "#" * R.WIDTH, ""]
    out.append(R.kv([
        ("생산 균주", sp.name),
        ("목적 산물", sp.product),
        ("목표 INCI", sp.inci),
        ("생산 방식", "액침 교반배양 (submerged fed-batch) -> 세포외 분비물 회수"),
        ("생산 규모", f"{C.EQUIPMENT['main'].tag} (작업용적 {C.EQUIPMENT['main'].working_volume_L:,.0f} L)"),
        ("균주 특성", sp.notes),
    ], indent=2, width=12))
    out.append("")
    out.append(R.table(
        ["#", "단계", "성격", "핵심 질문"],
        [["0", "균주 탐색·분리·동정", "R&D", "팔 수 있는 균을 찾았는가"],
         ["1", "균주은행 구축", "R&D", "3년 뒤에도 같은 균인가"],
         ["2", "배지 설계", "R&D", "싸고 규제에 맞는 밥인가"],
         ["3", "멸균", "생산", "F0 는 충분하고 영양소는 살았는가"],
         ["4", "종균 확대", "생산", "700만 배 키우는 동안 순수한가"],
         ["5", "본 발효", "생산", "역가·분자량·점도를 동시에 잡았는가"],
         ["6", "수확·균체분리", "생산", "케이크에 얼마나 빼앗겼는가"],
         ["7", "정제", "생산", "규격 순도까지 몇 % 남았는가"],
         ["8", "건조·분말화", "생산", "분자량과 용해성을 지켰는가"],
         ["9", "표준화·제형화", "생산", "라벨 역가를 배치마다 맞출 수 있는가"],
         ["10", "QC 규격시험", "품질", "숫자를 증명할 수 있는가"],
         ["11", "안전성·효능", "사업", "바이어가 받을 서류가 있는가"],
         ["12", "규제·등록", "사업", "어느 나라에 팔 수 있는가"],
         ["13", "물질수지·원가", "사업", "돈이 남는가"]], indent=2))
    out.append("")
    return "\n".join(out)


def final_summary(st: S.ProcessState) -> str:
    comp = S._solid_composition(st)
    pr = K.productivity(st.ferm)
    out = ["", "#" * R.WIDTH, R.pad(" 전공정 요약 ", R.WIDTH, "^"), "#" * R.WIDTH, ""]
    out.append(R.table(
        ["단계", "입력", "출력", "핵심 수치"],
        [["0-1 균주", "야생 자실체 24 점", f"생산균주 1 주 + MCB/WCB",
          f"순수배양 {len(st.isolates)} 주 중 선발"],
         ["2-3 배지·멸균", "글루코스 45 g/L 외", "무균 배지 7,000 L",
          f"F0 {st.steril['F0_min']:.0f} min, 영양소 손실 {st.steril['nutrient_loss_frac']*100:.1f}%"],
         ["4 종균", "WCB 1 vial (1 mL)", "접종원 560 L",
          f"X0 = {getattr(st,'_X0_main',0):.2f} g/L, 확대배수 7.0e6"],
         ["5 발효", "당 총 {:.0f} kg".format(
             C.MEDIUM_MAIN[0].g_per_L * 7.0 + pr["glucose_fed_kg"]),
          f"배양액 {pr['broth_volume_L']:,.0f} L", f"역가 {pr['titer_g_L']:.1f} g/L = {pr['product_in_broth_kg']:.0f} kg"],
         ["6 분리", f"{pr['product_in_broth_kg']:.1f} kg", f"{st.dsp_log[2]['out']:.1f} kg",
          f"회수율 {st.dsp_log[2]['out']/pr['product_in_broth_kg']*100:.1f}%"],
         ["7 정제", f"{st.dsp_log[2]['out']:.1f} kg", f"{st.dsp_log[6]['out']:.1f} kg",
          f"회수율 {st.dsp_log[6]['out']/st.dsp_log[2]['out']*100:.1f}%"],
         ["8 건조", f"{st.dsp_log[6]['out']:.1f} kg", f"분말 {comp['total_powder_kg']:.1f} kg",
          f"순도 {comp['purity_dry_pct']:.1f}%, 수분 {comp['moisture_pct']:.1f}%"],
         ["9-10 제품·QC", "벌크 분말", "3 그레이드 + CoA",
          "출하 승인" if st.released else "OOS 조사"],
         ["11-13 사업화", "안전성·효능·등록", "판매 가능 원료",
          f"COGS ${getattr(st,'_cogs_kg',0):,.0f}/kg"]], indent=2))
    out.append("")
    out.append(R.kv([
        ("발효 -> 분말 총 회수율", f"{st.mass_kg/pr['product_in_broth_kg']*100:.1f} %"),
        ("당 -> 제품 총괄 전환율",
         f"{pr['product_in_broth_kg']/(C.MEDIUM_MAIN[0].g_per_L*7.0+pr['glucose_fed_kg'])*100:.1f} %"),
        ("배치 생산량", f"{comp['total_powder_kg']:.1f} kg 분말"),
        ("배치 제조원가", f"${getattr(st,'_total_cost',0):,.0f}  (kg당 ${getattr(st,'_cogs_kg',0):,.0f})"),
        ("배치 매출(추정)", f"${getattr(st,'_revenue',0):,.0f}"),
        ("기술 개발 ~ 첫 출하", "약 30 ~ 42 개월 (균주 탐색 4 + 공정개발 10 + 안전성/효능 8 + 등록 12, 일부 병행)"),
    ], indent=2, width=24))
    out.append("")
    out.append(R.section("이 공정을 실제로 망하게 하는 것 TOP 5"))
    out.append(R.bullets([
        "1. 균주은행을 늦게 만든다 -> 2년 뒤 역가가 반토막, 원인 추적 불가",
        "2. 발효 역가만 쫓는다 -> 점도 폭증으로 다운스트림 회수율이 더 크게 떨어져 총수율 악화",
        "3. 수출 시장을 안 정하고 공정을 확정한다 -> 중국 신원료 등록 2~3년, 사업계획 전면 수정",
        "4. 안전성·효능 자료를 파일럿 검체로 받는다 -> 생산 배치로 전부 재시험",
        "5. 가동률을 과대 추정한다 -> 감가상각이 원가의 40%, 배치 수가 줄면 그대로 적자",
    ]))
    out.append("")
    return "\n".join(out)


def compare_strains() -> str:
    out = ["", "#" * R.WIDTH, R.pad(" 균주 트랙 비교 ", R.WIDTH, "^"), "#" * R.WIDTH, ""]
    rows = []
    for key in C.STRAINS:
        st, _ = S.run_pipeline(key)
        comp = S._solid_composition(st)
        pr = K.productivity(st.ferm)
        rows.append([C.STRAINS[key].name.split()[0] + " " + C.STRAINS[key].name.split()[1],
                     C.STRAINS[key].inci, f"{pr['titer_g_L']:.1f}",
                     f"{pr['peak_viscosity_cP']:.0f}",
                     f"{st.mass_kg/pr['product_in_broth_kg']*100:.0f}%",
                     f"{comp['total_powder_kg']:.0f}",
                     f"${getattr(st,'_cogs_kg',0):,.0f}"])
    out.append(R.table(["균주", "목표 INCI", "역가 g/L", "최고점도 cP",
                        "DSP 회수율", "분말 kg/batch", "COGS/kg"], rows, "<<>>>>>", indent=2))
    out.append("")
    out.append(R.bullets([
"Aspergillus(코직산)는 저분자라 배양액 점도가 안 오르고(산소전달 여유) 역가가 가장 높다",
        "  단 정제 경로는 UF 가 아니라 활성탄 탈색 -> 냉각결정화 -> 재결정으로 완전히 달라지고,",
        "  국내 배합한도(2%)·유럽 SCCS 안전성 논란 등 규제 리스크가 세 트랙 중 가장 크다",
        "Tremella 는 보습 소구가 강하고 HA 대체재 포지션이 선명하나, 이형성(효모형<->균사형) 제어가 어렵다",
        "Schizophyllum 은 역가·안전성·규제 균형이 가장 좋아 본 시뮬레이션의 기본 트랙으로 채택",
        "주의: 세 트랙 모두 동일한 다당 다운스트림 모델로 돌린 비교다. 코직산의 결정화 경로는",
        "  실제로는 회수율·원가가 모두 더 유리하므로 위 표의 COGS 는 보수적으로 읽어야 한다",
    ]))
    out.append("")
    return "\n".join(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description="균류 화장품 원료 전공정 시뮬레이터")
    ap.add_argument("--strain", default="schizophyllum", choices=list(C.STRAINS))
    ap.add_argument("--stages", default=None,
                    help="표시할 단계 (예: 5,6,7). 생략 시 전체")
    ap.add_argument("--seed", type=int, default=C.DEFAULT_SEED)
    ap.add_argument("--monte-carlo", type=int, default=0, metavar="N",
                    help="N 배치 반복 시뮬레이션")
    ap.add_argument("--compare", action="store_true", help="균주 3종 비교")
    ap.add_argument("--no-header", action="store_true")
    ap.add_argument("--out", default=None, help="결과를 파일로도 저장")
    a = ap.parse_args(argv)

    buf = StringIO()

    def w(x=""):
        buf.write(str(x) + "\n")

    if a.compare:
        w(compare_strains())
    elif a.monte_carlo:
        w(MC.render(MC.run_batches(a.monte_carlo, a.strain, seed=a.seed)))
    else:
        only = [int(x) for x in a.stages.split(",")] if a.stages else None
        st, results = S.run_pipeline(a.strain, a.seed, only)
        if not a.no_header and only is None:
            w(header(a.strain))
        for r in results:
            w(r.render())
        if only is None:
            w(final_summary(st))

    text = buf.getvalue()
    sys.stdout.write(text)
    if a.out:
        with open(a.out, "w", encoding="utf-8") as fh:
            fh.write(text)
        sys.stderr.write(f"\n[saved] {a.out}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
