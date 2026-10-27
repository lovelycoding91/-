"""터미널 리포트 유틸. 한글(동아시아 전각) 폭을 계산해 표를 정렬한다."""
import unicodedata
from typing import Iterable, List, Optional, Sequence

WIDTH = 96


def dw(s: str) -> int:
    """문자열의 표시 폭."""
    return sum(2 if unicodedata.east_asian_width(c) in "WF" else 1 for c in str(s))


def pad(s: str, n: int, align: str = "<") -> str:
    s = str(s)
    gap = max(n - dw(s), 0)
    if align == "<":
        return s + " " * gap
    if align == ">":
        return " " * gap + s
    left = gap // 2
    return " " * left + s + " " * (gap - left)


def rule(ch: str = "-", n: int = WIDTH) -> str:
    return ch * n


def banner(idx, title: str, objective: str = "") -> str:
    head = f" STAGE {idx} | {title} "
    line = "=" * WIDTH
    out = [line, pad(head, WIDTH, "^"), line]
    if objective:
        out.append(f"  목적 : {objective}")
        out.append(rule())
    return "\n".join(out)


def section(title: str) -> str:
    return f"\n  [{title}]\n  " + rule("-", WIDTH - 2)


def table(headers: Sequence[str], rows: Iterable[Sequence], aligns: Optional[str] = None,
          indent: int = 2) -> str:
    rows = [[("" if c is None else str(c)) for c in r] for r in rows]
    ncol = len(headers)
    aligns = aligns or "<" * ncol
    w = [dw(h) for h in headers]
    for r in rows:
        for i in range(ncol):
            w[i] = max(w[i], dw(r[i]))
    sp = " " * indent
    sep = sp + "  ".join("-" * x for x in w)
    out = [sp + "  ".join(pad(headers[i], w[i], "^") for i in range(ncol)), sep]
    for r in rows:
        out.append(sp + "  ".join(pad(r[i], w[i], aligns[i]) for i in range(ncol)))
    return "\n".join(out)


def kv(pairs, indent: int = 4, width: int = 34) -> str:
    sp = " " * indent
    return "\n".join(f"{sp}{pad(k, width)}: {v}" for k, v in pairs)


def bullets(items: Iterable[str], mark: str = "-", indent: int = 4) -> str:
    sp = " " * indent
    return "\n".join(f"{sp}{mark} {t}" for t in items)


def sparkline(values: Sequence[float], width: int = 56, height: int = 9,
              label: str = "", lo: Optional[float] = None,
              hi: Optional[float] = None) -> str:
    """ASCII 라인 플롯."""
    if not values:
        return ""
    lo = min(values) if lo is None else lo
    hi = max(values) if hi is None else hi
    if hi - lo < 1e-12:
        hi = lo + 1.0
    n = len(values)
    cols = [values[int(i * (n - 1) / (width - 1))] for i in range(width)]
    grid = [[" "] * width for _ in range(height)]
    for x, v in enumerate(cols):
        y = int(round((v - lo) / (hi - lo) * (height - 1)))
        grid[height - 1 - min(max(y, 0), height - 1)][x] = "*"
    out = []
    for i, row in enumerate(grid):
        tick = hi if i == 0 else (lo if i == height - 1 else None)
        axis = f"{tick:8.2f} |" if tick is not None else " " * 8 + " |"
        out.append("    " + axis + "".join(row))
    out.append("    " + " " * 8 + " +" + "-" * width)
    if label:
        out.append("    " + " " * 9 + " " + label)
    return "\n".join(out)


def verdict(ok: bool, ok_text: str = "적합", ng_text: str = "부적합") -> str:
    return f"[PASS] {ok_text}" if ok else f"[FAIL] {ng_text}"
