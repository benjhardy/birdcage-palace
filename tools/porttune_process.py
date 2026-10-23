"""Process diagram "Tuning with circuitRF": the port-tuning flow beside the eigenmode loop, with the measured costs.

    python tools/porttune_process.py                      # docs/img/circuitrf-tuning-process.svg + .png (1600 px wide)
    python tools/porttune_process.py --portrait out.png   # 1080 x 1350 copy for a post

Every number is from the README (sections "Port tuning with circuitRF" and "Results"), results/lp16_3T.json
(eigenmode tune history) and results/lp16_porttune.json; the 24 s circuitRF time is a re-run of
`birdcage porttune tune` on results/birdcage16_unloaded.s16p (2026-10-23, same machine).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
INK, INK2, MUTED, LINE = "#0b0b0b", "#52514e", "#8a8984", "#d6d4ce"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
FILL = {"palace": "#eef4fc", "crf": "#fdf0ea", "result": "#f4f3f0"}
EDGE = {"palace": BLUE, "crf": ORANGE, "result": INK2}

PORT = [
    ("palace", "1  Palace builds the bare coil",
     "Eighth model. Each rung gap is a 500 Ω lumped port;\nno capacitors in the 3D model.", "built once"),
    ("palace", "2  Four symmetry runs (driven sweep)",
     "x = 0 and y = 0 planes: PEC or PMC (em, me, mm, ee);\nz = 0 PEC. 4 gap excitations per run.",
     "0.16 M unknowns each; 33 min in all"),
    ("palace", "3  Assemble the 16-port S-matrix",
     "The four runs give every entry (character orthogonality).\nLossless to |S Sᴴ − I| < 5.1e-5.", "seconds"),
    ("crf", "4  circuitRF closes the gaps",
     "C + 0.68 Ω ESR on every gap. Capacitive-divider feed\n(Cd, Cp, Cm) on rungs 0 and 4 to two 50 Ω ports.", ""),
    ("crf", "5  Minimax optimiser",
     "Goals: |S11|, |S22|, |S21| ≤ −20 dB\nat 127.73 MHz ± 2 kHz.", "102 + 28 evaluations; 24 s in all"),
    ("result", "6  Tuned C values",
     "Ct 1.844 pF (12 rungs), Ci 1.841 pF (rungs 2, 10),\nCd 1.768 pF + Cp 24.70 pF, Cm 20.6 pF.", ""),
    ("palace", "7  Verify in Palace: z-half eigen solve",
     "All 16 rungs, feeds disconnected. Shift from uniform C:\n+0.169 MHz in Palace, +0.169 MHz from the 16-port.",
     "0.58 M unknowns, 2.3 min"),
]
EIGEN = [
    ("palace", "A  Palace builds the full coil",
     "Ladder first guess in every gap:\nC = 2.807 pF.", ""),
    ("palace", "B  Eigen solve, k = 1 pair",
     "Full model, 1.12 M unknowns;\n13.6–20.4 min per solve.", ""),
    ("result", "C  Secant step on 1/f²",
     "k = 1: 111.08 → 128.04\n→ 127.73 MHz.", ""),
    ("result", "D  Tuned: 1.828 pF in every gap",
     "3 solves, 48 min of Palace time.\nNo feed or match network.", ""),
]
TITLE = "Tuning with circuitRF"
SUB = ("16-rung low-pass birdcage, D 28 cm, 3 T, unloaded; production mesh density (5 / 2.5 / 35 mm, order 2), "
       "4 MPI ranks. Simulation only.")
FOOT = ("Palace wall times; the port runs and the loop solves ran on a CPU shared at times with another job. "
        "Sources: README, results/lp16_3T.json, results/lp16_porttune.json.")


def box(ax, x, y, w, h, kind, head, body, cost, fs):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.012", fc=FILL[kind], ec=EDGE[kind],
                                lw=1.4, transform=ax.transAxes, mutation_aspect=ax.figure.get_figwidth() / ax.figure.get_figheight()))
    px = 0.012
    H = ax.figure.get_figheight()
    ax.text(x + px, y + h - 0.12 / H, head, fontsize=fs * 1.08, fontweight="bold", color=INK, va="top", transform=ax.transAxes)
    ax.text(x + px, y + h - (0.12 + fs * 0.024) / H, body, fontsize=fs * 0.93, color=INK2, va="top", linespacing=1.35,
            transform=ax.transAxes)
    if cost:
        ax.text(x + w - px, y + h - 0.12 / H, cost, fontsize=fs * 0.9, color=EDGE[kind], fontweight="bold", ha="right", va="top",
                transform=ax.transAxes)


def arrow(ax, p, q, color=MUTED, rad=0.0, lw=1.6, ls="-"):
    ax.add_patch(FancyArrowPatch(p, q, transform=ax.transAxes, arrowstyle="-|>", mutation_scale=14, color=color, lw=lw,
                                 connectionstyle=f"arc3,rad={rad}", linestyle=ls, shrinkA=0, shrinkB=0))


def column(ax, items, x, w, y_top, h, gap, fs):
    ys = []
    for i, (kind, head, body, cost) in enumerate(items):
        y = y_top - h - i * (h + gap)
        box(ax, x, y, w, h, kind, head, body, cost, fs)
        ys.append(y)
    for i in range(len(items) - 1):
        arrow(ax, (x + w / 2, ys[i]), (x + w / 2, ys[i + 1] + h))
    return ys


def legend(ax, x, y, fs, dx=0.15):
    for i, (k, lab) in enumerate((("palace", "Palace (3D FEM)"), ("crf", "circuitRF (circuit)"), ("result", "values"))):
        xx = x + i * dx
        ax.add_patch(FancyBboxPatch((xx, y), 0.018, 0.016, boxstyle="round,pad=0,rounding_size=0.004", fc=FILL[k], ec=EDGE[k],
                                    lw=1.2, transform=ax.transAxes))
        ax.text(xx + 0.024, y + 0.008, lab, fontsize=fs * 0.9, color=INK2, va="center", transform=ax.transAxes)


def landscape(out_svg: Path, out_png: Path):
    fig = plt.figure(figsize=(16, 10), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    fs = 12
    ax.text(0.03, 0.965, TITLE, fontsize=24, fontweight="bold", color=INK, va="top")
    ax.text(0.03, 0.915, SUB, fontsize=11.5, color=INK2, va="top")
    legend(ax, 0.03, 0.865, fs)
    # port tuning, left
    xl, wl, h, gap, top = 0.03, 0.50, 0.09, 0.016, 0.83
    ax.text(xl, top + 0.004, "Port tuning (circuitRF)", fontsize=14, fontweight="bold", color=ORANGE, va="bottom")
    ys = column(ax, PORT, xl, wl, top - 0.008, h, gap, fs)
    # optimiser loop: 5 -> 4
    xr = xl + wl
    arrow(ax, (xr, ys[4] + h * 0.5), (xr, ys[3] + h * 0.5), color=ORANGE, rad=0.55, lw=2.0)
    ax.text(xr + 0.035, (ys[3] + ys[4] + h) / 2, "new C values:\none circuit\nevaluation", fontsize=fs * 0.95,
            color=ORANGE, fontweight="bold", va="center", ha="left")
    # eigenmode loop, right
    xe, we, he = 0.66, 0.31, 0.098
    ax.text(xe, top + 0.004, "Eigenmode loop (default)", fontsize=14, fontweight="bold", color=BLUE, va="bottom")
    ye = column(ax, EIGEN, xe, we, top - 0.008, he, 0.045, fs)
    arrow(ax, (xe, ye[2] + he * 0.5), (xe, ye[1] + he * 0.5), color=BLUE, rad=-0.45, lw=2.0)
    ax.text(xe - 0.05, (ye[1] + ye[2] + he) / 2, "new C:\none 3D\nsolve", fontsize=fs * 0.95, color=BLUE,
            fontweight="bold", va="center", ha="right")
    # cost comparison card
    yc = 0.075
    ax.add_patch(FancyBboxPatch((xe, yc), we, 0.175, boxstyle="round,pad=0,rounding_size=0.012", fc="white", ec=LINE, lw=1.2,
                                transform=ax.transAxes, mutation_aspect=1.6))
    ax.text(xe + 0.012, yc + 0.162, "Cost of a capacitor change", fontsize=fs * 1.05, fontweight="bold", color=INK, va="top")
    rows = [("Eigenmode loop", "one 3D solve, 13.6–20.4 min", BLUE),
            ("Port tuning", "one circuit evaluation;", ORANGE),
            ("", "130 evaluations in 24 s", ORANGE),
            ("Port tuning, once", "4 eighth runs, 33 min", INK2)]
    for i, (a, b, c) in enumerate(rows):
        yy = yc + 0.122 - i * 0.03
        ax.text(xe + 0.012, yy, a, fontsize=fs * 0.95, color=c, fontweight="bold", va="top")
        ax.text(xe + 0.135, yy, b, fontsize=fs * 0.95, color=INK2, va="top")
    ax.text(0.03, 0.022, FOOT, fontsize=10, color=MUTED, va="bottom")
    plt.rcParams["svg.fonttype"] = "path"
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    fig.savefig(out_svg, facecolor="white")
    fig.savefig(out_png, dpi=100, facecolor="white")
    plt.close(fig)


def portrait(out_png: Path):
    fig = plt.figure(figsize=(10.8, 13.5), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    fs = 12.5
    ax.text(0.05, 0.972, TITLE, fontsize=27, fontweight="bold", color=INK, va="top")
    ax.text(0.05, 0.932, "16-rung low-pass birdcage, D 28 cm, 3 T, unloaded.\nProduction mesh density, 4 MPI ranks. Simulation only.",
            fontsize=12, color=INK2, va="top", linespacing=1.4)
    legend(ax, 0.05, 0.868, fs, dx=0.25)
    xl, wl, h, gap, top = 0.05, 0.76, 0.078, 0.014, 0.85
    ys = column(ax, PORT, xl, wl, top, h, gap, fs)
    xr = xl + wl
    arrow(ax, (xr, ys[4] + h * 0.5), (xr, ys[3] + h * 0.5), color=ORANGE, rad=0.6, lw=2.0)
    ax.text(xr + 0.05, (ys[3] + ys[4] + h) / 2, "new C:\none circuit\nevaluation", fontsize=fs * 0.95, color=ORANGE,
            fontweight="bold", va="center", ha="left")
    yb = ys[-1] - 0.03
    ax.plot([0.05, 0.95], [yb, yb], color=LINE, lw=1.0, transform=ax.transAxes)
    ax.text(0.05, yb - 0.012, "The eigenmode loop it replaces", fontsize=14, fontweight="bold", color=BLUE, va="top")
    lines = [("Each C change:", "one full 3D eigen solve, 1.12 M unknowns, 13.6–20.4 min"),
             ("16-rung tune:", "2.807 → 1.813 → 1.828 pF; 111.08 → 128.04 → 127.73 MHz"),
             ("", "3 solves, 48 min of Palace time, no feed network"),
             ("Port tuning:", "4 eighth runs once (33 min), then 130 circuit evaluations in 24 s")]
    for i, (a, b) in enumerate(lines):
        yy = yb - 0.048 - i * 0.026
        ax.text(0.05, yy, a, fontsize=fs, color=INK, fontweight="bold", va="top")
        ax.text(0.22, yy, b, fontsize=fs, color=INK2, va="top")
    ax.text(0.05, 0.012, "Palace wall times, CPU shared at times. Sources: birdcage-palace README and results/.", fontsize=10,
            color=MUTED, va="bottom")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    fig.savefig(out_png, dpi=100, facecolor="white")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "docs" / "img" / "circuitrf-tuning-process"))
    ap.add_argument("--portrait", default=None, help="also write a 1080 x 1350 PNG here")
    a = ap.parse_args()
    plt.rcParams.update({"font.family": "DejaVu Sans"})
    o = Path(a.out)
    o.parent.mkdir(parents=True, exist_ok=True)
    landscape(o.with_suffix(".svg"), o.with_suffix(".png"))
    print("wrote", o.with_suffix(".svg"), o.with_suffix(".png"))
    if a.portrait:
        portrait(Path(a.portrait))
        print("wrote", a.portrait)


if __name__ == "__main__":
    main()
