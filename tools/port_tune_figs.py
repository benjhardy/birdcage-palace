"""Figures for the port-tuning step: drive-port S-parameters before/after, the circuitRF trajectory, the C values.

    python tools/port_tune_figs.py --unloaded runs/pt_e/circuitrf/tune_unloaded.json \
        --loaded runs/ptl_e/circuitrf/tune_loaded.json [--also <dir>]
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402,F401

from birdcage_palace import plots  # noqa: E402
import port_tune as pt  # noqa: E402

LAB = {"Ct": "C, 12 rungs", "Ci": "C, rungs 2 and 10\n(isolation)", "Cd1": "Cd, drive gaps (top)", "Cd2": "Cd, drive gap 2",
       "Cm1": "Cm, series to port", "Cm2": "Cm 2", "Cp1": "Cp, tap", "Cp2": "Cp 2"}


def fig_sparams(res: list[tuple[str, dict]], out: Path, figsize=(11.5, 4.2), fs=1.0) -> Path:
    fig, axs = plt.subplots(1, len(res), figsize=figsize)
    axs = np.atleast_1d(axs)
    for ax, (title, R) in zip(axs, res):
        f0 = R["f0_hz"] / 1e6
        for tag, ls, alpha in (("before", "--", 0.9), ("after", "-", 1.0)):
            sw = R[f"sweep_{tag}"]
            f = np.array(sw["f_hz"]) / 1e6
            ax.plot(f, sw["S11_db"], ls, color=plots.BLUE, lw=1.6, alpha=alpha, label=f"|S11| {tag}")
            ax.plot(f, sw["S21_db"], ls, color=plots.ORANGE, lw=1.3, alpha=alpha, label=f"|S21| {tag}")
        ax.axvline(f0, color=plots.MUTED, lw=0.8, ls=":")
        ax.axhline(-20, color=plots.MUTED, lw=0.8, ls=":")
        ax.set_xlabel("frequency (MHz)", fontsize=10 * fs)
        ax.set_ylabel("dB", fontsize=10 * fs)
        ax.set_ylim(-50, 2)
        ax.set_title(title, loc="left", fontsize=10.5 * fs)
        ax.grid(color=plots.GRID, lw=0.6)
        ax.tick_params(labelsize=8.5 * fs)
    axs[0].legend(frameon=False, fontsize=8 * fs, loc="lower left")
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight", dpi=170)
    plt.close(fig)
    return out


def fig_trajectory(R: dict, out: Path, title: str, figsize=(11.5, 4.0), fs=1.0) -> Path:
    tr = R["trajectory"]
    # drop the tail where the seeded run no longer changes (the optimizer has stopped moving)
    while len(tr) > 2 and abs(tr[-1]["cost"] - tr[-2]["cost"]) < 1e-9 and tr[-1]["vals"] == tr[-2]["vals"]:
        tr = tr[:-1]
    it = [t["iter"] for t in tr]
    fig, axs = plt.subplots(1, 2, figsize=figsize)
    ax = axs[0]
    for g, col in (("match1", plots.BLUE), ("match2", plots.AQUA), ("isol", plots.ORANGE)):
        ax.plot(it, [t["goals"][g] for t in tr], "-o", ms=4, color=col, label={"match1": "worst |S11| in band",
                                                                                "match2": "worst |S22| in band",
                                                                                "isol": "worst |S21| in band"}[g])
    ax.axhline(-20, color=plots.MUTED, ls=":", lw=0.8)
    ax.set_xlabel("circuitRF minimax iteration", fontsize=10 * fs)
    ax.set_ylabel("dB", fontsize=10 * fs)
    ax.legend(frameon=False, fontsize=8 * fs)
    ax.grid(color=plots.GRID, lw=0.6)
    ax.set_title(title, loc="left", fontsize=10.5 * fs)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    ax = axs[1]
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    def gap(v):  # capacitance across a drive gap: Cd in series with the tap Cp (divider feed)
        return 1 / (1 / v["Cd1"] + 1 / v["Cp1"]) if "Cp1" in v else v["Cd1"]
    series = (("C, 12 rungs", lambda v: v["Ct"], plots.INK2), ("C, rungs 2 and 10", lambda v: v["Ci"], plots.VIOLET),
              ("drive gaps, Cd in series with Cp", gap, plots.BLUE))
    for lab, fn, col in series:
        ax.plot(it, [fn(t["vals"]) * 1e12 for t in tr], "-o", ms=3, lw=1.5, color=col, label=lab)
    if "Cp1" in tr[0]["vals"]:
        ax3 = ax.twinx()
        ax3.plot(it, [t["vals"]["Cm1"] * 1e12 for t in tr], "--", lw=1.3, color=plots.ORANGE, label="Cm, series to port (right)")
        ax3.plot(it, [t["vals"]["Cp1"] * 1e12 for t in tr], ":", lw=1.6, color=plots.AQUA, label="Cp, tap (right)")
        ax3.set_ylabel("feed capacitors (pF)", fontsize=9 * fs)
        ax3.spines["right"].set_visible(True)
        ax3.legend(frameon=False, fontsize=7.5 * fs, loc="center right")
    ax.set_xlabel("circuitRF minimax iteration", fontsize=10 * fs)
    ax.set_ylabel("capacitance across the gap (pF)", fontsize=10 * fs)
    ax.legend(frameon=False, fontsize=7.5 * fs, loc="upper left")
    ax.grid(color=plots.GRID, lw=0.6)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight", dpi=170)
    plt.close(fig)
    return out


def fig_cvalues(res: list[tuple[str, dict]], out: Path, figsize=(12.5, 4.2)) -> Path:
    """Left: capacitance across each gap (drive gaps: Cd in series with the tap Cp). Right: the feed capacitors."""
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=figsize, gridspec_kw={"width_ratios": [2.2, 1]})
    w = 0.38
    labs = ["eigenmode loop\n(unloaded,\nall rungs)", "port tuning,\nall rungs\nequal", "12 rungs", "rungs 2, 10\n(isolation)",
            "drive gap 1\n(Cd in series\nwith Cp)", "drive gap 2"]
    flabs = ["Cd (top)", "Cp (tap)", "Cm (series\nto port)"]
    for j, (title, R) in enumerate(res):
        v = R["raw"]["vals"]
        ce = R["c_eigen"] * 1e12 if j == 0 else np.nan  # the eigenmode loop tuned the unloaded coil only
        gap = [1e12 / (1 / v[f"Cd{k}"] + 1 / v[f"Cp{k}"]) if "Cp1" in v else v[f"Cd{k}"] * 1e12 for k in (1, 2)]
        vals = [ce, R.get("c_resonance_only", np.nan) * 1e12, v["Ct"] * 1e12, v["Ci"] * 1e12] + gap
        x = np.arange(len(vals))
        col = [plots.BLUE, plots.ORANGE][j]
        ax.bar(x + (j - 0.5) * w, vals, width=w, color=col, label=title)
        for xi, vi in zip(x + (j - 0.5) * w, vals):
            if np.isfinite(vi):
                ax.text(xi, vi + 0.03, f"{vi:.3f}", ha="center", va="bottom", fontsize=7, rotation=90, color=plots.INK)
        if "Cp1" in v:
            fv = [v["Cd1"] * 1e12, v["Cp1"] * 1e12, v["Cm1"] * 1e12]
            x2 = np.arange(3)
            ax2.bar(x2 + (j - 0.5) * w, fv, width=w, color=col)
            for xi, vi in zip(x2 + (j - 0.5) * w, fv):
                ax2.text(xi, vi * 1.08, f"{vi:.2f}", ha="center", va="bottom", fontsize=7, rotation=90, color=plots.INK)
    ax.set_xticks(np.arange(len(labs)))
    ax.set_xticklabels(labs, fontsize=8)
    ax.set_ylabel("capacitance across the gap (pF)")
    ax.set_ylim(0, 2.6)
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")
    ax2.set_xticks(np.arange(3))
    ax2.set_xticklabels(flabs, fontsize=8)
    ax2.set_yscale("log")
    ax2.set_ylim(0.5, 200)
    ax2.set_ylabel("pF (drive port 1)")
    ax2.set_title("capacitive-divider feed", fontsize=9.5, loc="left")
    for a_ in (ax, ax2):
        a_.grid(axis="y", color=plots.GRID, lw=0.6)
        a_.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight", dpi=170)
    plt.close(fig)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--unloaded", required=True)
    ap.add_argument("--loaded", default=None)
    ap.add_argument("--also", action="append", default=[])
    a = ap.parse_args()
    img = ROOT / "docs" / "img"
    U = json.loads(Path(a.unloaded).read_text())
    res = [("unloaded (capacitor ESR 0.68 Ω)", U)]
    if a.loaded:
        res.append(("loaded (saline phantom)", json.loads(Path(a.loaded).read_text())))
    figs = [fig_sparams(res, img / "porttune_sparams.png"),
            fig_trajectory(res[-1][1], img / "porttune_iterations.png", f"circuitRF tuning, {res[-1][0]}"),
            fig_cvalues(res, img / "porttune_cvalues.png")]
    for d in a.also:
        Path(d).mkdir(parents=True, exist_ok=True)
        for f in figs:
            shutil.copy(f, Path(d) / f.name)
    print("wrote", ", ".join(f.name for f in figs))


if __name__ == "__main__":
    main()
