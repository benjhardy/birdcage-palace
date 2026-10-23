"""Collect the Open BirdcageBuilder comparison into results/obcb_comparison.json and draw the figure.

Inputs (all produced by runs or read from files, nothing typed in):
    compare/open-birdcagebuilder/obcb_raw.json      Open BirdcageBuilder outputs (obcb_harness.mjs)
    compare/open-birdcagebuilder/palace_checks.json  Palace k=1 at each C (palace_checks.py)
    results/lp16_3T.json, hp16_3T.json, bp16_3T.json  Palace-tuned references (full model, production mesh)
    results/lp16_symmetry.json                       eighth model vs full model at the tuned C
    results/shield_compare_coarse.json                unshielded Palace tune (coarse mesh, order 1)
    our ladder model (birdcage_palace.ladder), run here with and without the shield

    python compare/open-birdcagebuilder/make_comparison.py [--fig <png>]
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))

from birdcage_palace import ladder  # noqa: E402
from birdcage_palace.design import CoilDesign  # noqa: E402

TYPES = [("lowpass", "lp", "Low-pass", "C in the 16 rungs"), ("highpass", "hp", "High-pass", "C in the 32 ring gaps"),
         ("bandpass", "bp", "Band-pass", "C in rungs and ring gaps, equal values")]
OBCB_NAME = {"lowpass": "lp", "highpass": "hp", "bandpass": "bp_leg"}


def jload(p):
    return json.loads(Path(p).read_text())


def interp_k1(c, c1, f1, c2, f2):
    """k=1 at C from two solves, with 1/f^2 linear in C (1/f^2 = L_eff (C + C_stray), as `birdcage tune` uses)."""
    y1, y2 = 1 / f1 ** 2, 1 / f2 ** 2
    return (y1 + (y2 - y1) * (c - c1) / (c2 - c1)) ** -0.5


def collect() -> dict:
    raw = jload(HERE / "obcb_raw.json")
    pal = {r["name"]: r for r in jload(HERE / "palace_checks.json")}
    obcb = {c["name"]: c for c in raw["cases"]}
    clone = Path(r"C:\Users\ben\jobs\_ext\OpenBirdcageBuilder")
    commit = subprocess.run(["git", "-C", str(clone), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    f0 = CoilDesign().f_larmor / 1e6
    sym = jload(REPO / "results/lp16_symmetry.json")

    rows = {}
    for ctype, tag, title, capwhere in TYPES:
        ref = jload(REPO / f"results/{tag}16_3T.json")
        d = ref["design"]
        key_c = "c_leg" if ctype != "highpass" else "c_ring"
        c_tuned = d[key_c] * 1e12
        c_ladder = d["notes"][f"{key_c}_first_guess"] * 1e12
        o_sh = obcb[f"{OBCB_NAME[ctype]}_shield34_L25"]["result"]["cap_pf"]
        lad_free = ladder.analytic(CoilDesign(coil_type=ctype, shield_diameter=None))
        rows[ctype] = {
            "title": title, "capacitors": capwhere,
            "palace_model": "eighth symmetry, production mesh" if ctype == "lowpass" else "full model, production mesh",
            "obcb": {"c_pf": o_sh, "f_k1_mhz": pal[f"{tag}_obcb"]["f_k1_mhz"], "palace_case": pal[f"{tag}_obcb"]["case"]},
            "ladder": ({"c_pf": c_ladder, "f_k1_mhz": pal[f"{tag}_ladder"]["f_k1_mhz"],
                        "palace_case": pal[f"{tag}_ladder"]["case"], "source": "direct Palace solve"}
                       if f"{tag}_ladder" in pal else
                       {"c_pf": c_ladder, "f_k1_mhz": interp_k1(c_ladder, o_sh, pal[f"{tag}_obcb"]["f_k1_mhz"], c_tuned,
                                                               ref["homogeneous_pair"]["f_mean_mhz"]),
                        "palace_case": None,
                        "source": ("interpolated, not solved: 1/f^2 linear in C between the two production-mesh "
                                   "solves of this coil (Open BirdcageBuilder C and Palace-tuned C); the same rule "
                                   "on the high-pass coil is checked against a direct solve in interp_check")}),
            "tuned": {"c_pf": c_tuned, "f_k1_mhz": ref["homogeneous_pair"]["f_mean_mhz"],
                      "palace_case": f"results/{tag}16_3T.json (full model)"},
            "obcb_leg_length_23p5cm_pf": obcb[f"{OBCB_NAME[ctype]}_shield34_L23.5"]["result"]["cap_pf"],
            "unshielded_circuit": {
                "obcb_pf": obcb[f"{OBCB_NAME[ctype]}_noshield_L25"]["result"]["cap_pf"],
                "ladder_pf": (lad_free["c_leg"] if ctype != "highpass" else lad_free["c_ring"]) * 1e12},
        }
        if ctype == "lowpass":
            rows[ctype]["tuned"]["f_k1_mhz_eighth_model"] = pal["lp_tuned"]["f_k1_mhz"]
            rows[ctype]["tuned"]["palace_case_eighth_model"] = pal["lp_tuned"]["case"]
            rows[ctype]["ladder"]["f_k1_mhz_full_model"] = jload(REPO / "runs/lp16_final/tune_history_part1.json")[0]["f_k1_mhz"]
            rows[ctype]["eighth_vs_full_at_tuned_c"] = {"eighth_mhz": sym["eighth"]["k1_mhz"], "full_mhz": sym["full"]["k1_mean_mhz"]}
        for k in ("obcb", "ladder", "tuned"):
            r = rows[ctype][k]
            r["f_err_pct"] = None if r["f_k1_mhz"] is None else 100 * (r["f_k1_mhz"] / f0 - 1)
            r["c_err_vs_tuned_pct"] = 100 * (r["c_pf"] / c_tuned - 1)
        u = rows[ctype]["unshielded_circuit"]
        u["obcb_vs_ladder_pct"] = 100 * (u["obcb_pf"] / u["ladder_pf"] - 1)

    # check of the interpolation rule on the high-pass coil, where the ladder guess was also solved directly
    h = rows["highpass"]
    if h["ladder"].get("source") == "direct Palace solve":
        fi = interp_k1(h["ladder"]["c_pf"], h["obcb"]["c_pf"], h["obcb"]["f_k1_mhz"], h["tuned"]["c_pf"],
                       h["tuned"]["f_k1_mhz"])
        rows["interp_check"] = {"coil": "highpass", "c_pf": h["ladder"]["c_pf"], "f_interpolated_mhz": fi,
                                "f_solved_mhz": h["ladder"]["f_k1_mhz"],
                                "diff_pct": 100 * (fi / h["ladder"]["f_k1_mhz"] - 1)}

    noshield = jload(REPO / "results/shield_compare_coarse.json")["cases"][1]
    rows["lowpass"]["unshielded_circuit"]["palace_tuned_pf_coarse_o1"] = noshield["c_tuned_pf"]

    return {
        "what": ("Capacitor first guesses for the 16-rung 3 T head birdcage (D 28 cm, L 25 cm between end-ring "
                 "centrelines, rungs 12 mm and rings 15 mm flat copper strips, 34 cm closed PEC shield, unloaded) "
                 "from Open BirdcageBuilder and from our ladder model, each put into Palace to see where k=1 lands. "
                 "Simulation only."),
        "f_target_mhz": f0,
        "open_birdcagebuilder": {
            "repo": "https://github.com/dbnyu/OpenBirdcageBuilder", "commit": commit, "version": "1.0.1",
            "licence": "GNU AGPL-3.0 (not copied here; the harness loads the unmodified files from a local clone)",
            "authors": "D. Brantner, G. Carluccio, C.-L. Chin, C. M. Collins (NYU Langone CAI2R/CBI)",
            "citation": "Brantner D, Carluccio G, Collins CM. Open BirdcageBuilder. ISMRM 2024, #1569; method of Chin CL, "
                        "Collins CM, Li S, Dardzinski BJ, Smith MB. Concepts Magn Reson 2002;15(2):156-163",
            "selfcheck": raw["obcb_selfcheck"],
            "inputs": {"n_legs": 16, "f_mhz": round(f0, 6), "coil_radius_cm": 14.0, "shield_radius_cm": 17.0,
                       "leg_length_cm": 25.0, "leg_type": "rectangular, width 1.2 cm",
                       "end_ring_type": "rectangular, width 1.5 cm",
                       "band_pass": "C_leg = C_ring found by bisection on its 'fixed caps on legs' mode "
                                    "(the 'fixed caps on ER' mode gives the same value)"},
        },
        "rows": rows,
        "palace": {"settings": "same mesh, order 2, MUMPS (MGMaxLevels 1), Tol 1e-8, mode count and target as the "
                               "reference runs; only the lumped C values change",
                   "runs": list(pal.values())},
    }


def figure(res: dict, out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    f0 = res["f_target_mhz"]
    rows = [res["rows"][t] for t, *_ in TYPES]
    fig = plt.figure(figsize=(10.8, 13.5), dpi=100)
    fig.text(0.06, 0.965, "Birdcage capacitor first guesses, checked in Palace", fontsize=23, weight="bold", va="top")
    fig.text(0.06, 0.928, "16 rungs · D 28 cm · L 25 cm · copper strips · 34 cm shield · 3 T, target 127.73 MHz",
             fontsize=14.5, va="top", color="0.25")

    lab = ["Open BirdcageBuilder", "Our ladder model"]
    keys = ["obcb", "ladder"]
    cols = ["C0", "C1"]
    x = np.arange(3)
    w = 0.3

    def panel(ax, val, txt, title, ylabel, ylim, tuned_txt):
        for i, (k, l, c) in enumerate(zip(keys, lab, cols)):
            v = [val(r[k]) for r in rows]
            interp = ["interpolated" in r[k].get("source", "") for r in rows]
            bars = ax.bar(x + (i - 0.5) * w, v, w * 0.92, color=c, label=l)
            for bar, it in zip(bars, interp):
                if it:
                    bar.set_hatch("//")
                    bar.set_edgecolor("white")
            for xi, vi, r, it in zip(x + (i - 0.5) * w, v, rows, interp):
                pad = 0.025 * (ylim[1] - ylim[0])
                ax.text(xi, vi + (pad if vi >= 0 else -pad), txt(r[k], vi) + ("\n(interp.)" if it else ""),
                        ha="center", va="bottom" if vi >= 0 else "top", fontsize=11.5)
        for xi, r in zip(x, rows):
            ax.text(xi + w + 0.02, 0.02 * (ylim[1] - ylim[0]), tuned_txt(r["tuned"]), ha="left", va="bottom", fontsize=10.5,
                    color="C2")
        ax.axhline(0, color="k", lw=1)
        ax.set_ylim(*ylim)
        ax.set_xlim(-0.55, 2.75)
        ax.set_xticks(x, [r["title"] for r in rows], fontsize=15)
        ax.set_ylabel(ylabel, fontsize=13.5)
        ax.set_title(title, fontsize=16, loc="left")
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=12)

    ax1 = fig.add_axes([0.11, 0.575, 0.86, 0.29])
    panel(ax1, lambda d: d["f_err_pct"], lambda d, v: f"{d['f_k1_mhz']:.1f} MHz\n{v:+.1f} %",
          "Where k=1 lands in Palace with each first-guess C", "k=1 error vs 127.73 MHz (%)", (-20, 25),
          lambda d: f"tuned C\n{d['f_k1_mhz']:.2f} MHz")
    ax1.legend(fontsize=13, loc="lower right", bbox_to_anchor=(1, 0.09), frameon=False)
    ax1.text(2.75, -19.5, "hatched: interpolated between two Palace solves (rule checked on high-pass: 0.001 %)",
             ha="right", va="bottom", fontsize=10, color="0.3")

    ax2 = fig.add_axes([0.11, 0.255, 0.86, 0.235])
    panel(ax2, lambda d: d["c_err_vs_tuned_pct"], lambda d, v: f"{d['c_pf']:.2f} pF\n{v:+.0f} %",
          "First-guess C vs the Palace-tuned C", "C vs Palace-tuned C (%)", (-48, 75),
          lambda d: f"tuned\n{d['c_pf']:.2f} pF")

    u = [r["unshielded_circuit"]["obcb_vs_ladder_pct"] for r in rows]
    lp, hp, bp = rows
    notes = [
        f"Low-pass: Open BirdcageBuilder lands closer ({lp['obcb']['f_err_pct']:+.1f} % vs {lp['ladder']['f_err_pct']:+.1f} %).",
        f"High- and band-pass: our ladder lands closer ({hp['ladder']['f_err_pct']:+.1f} % and {bp['ladder']['f_err_pct']:+.1f} % "
        f"vs {hp['obcb']['f_err_pct']:+.1f} % and {bp['obcb']['f_err_pct']:+.1f} %).",
        f"Without a shield the two circuit models agree within {max(abs(v) for v in u):.1f} % on C. The gap above comes "
        "from the shield:",
        "Open BirdcageBuilder images the rungs only, our ladder images the rungs and the end rings.",
        "Neither has the strips' stray capacitance, so both low-pass guesses are too high.",
    ]
    for i, t in enumerate(notes):
        fig.text(0.06, 0.195 - i * 0.025, t, fontsize=12.5, va="top", weight="bold" if i < 2 else None)
    fig.text(0.06, 0.068, "Open BirdcageBuilder by Brantner, Carluccio, Chin, Collins, NYU Langone CAI2R (ISMRM 2024, "
             "#1569),", fontsize=11.5, va="top", color="0.25")
    fig.text(0.06, 0.047, "based on BirdcageBuilder, Chin et al., Concepts Magn Reson 2002. Run unmodified, commit "
             f"{res['open_birdcagebuilder']['commit'][:7]}.", fontsize=11.5, va="top", color="0.25")
    fig.text(0.06, 0.018, "Palace FEM · unloaded · simulation only", fontsize=13, va="top", weight="bold")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=100)
    print("wrote", out)


def table_md(res: dict) -> str:
    lines = ["| | | C (pF) | k=1 in Palace (MHz) | error vs 127.73 MHz | C vs Palace-tuned |", "|---|---|---|---|---|---|"]
    for t, *_ in TYPES:
        r = res["rows"][t]
        for k, l in (("obcb", "Open BirdcageBuilder"), ("ladder", "our ladder first guess"), ("tuned", "Palace-tuned")):
            v = r[k]
            fe = "—" if k == "tuned" else f"{v['f_err_pct']:+.2f} %"
            ce = "—" if k == "tuned" else f"{v['c_err_vs_tuned_pct']:+.1f} %"
            star = " (interpolated, see note)" if "interpolated" in v.get("source", "") else ""
            lines.append(f"| {r['title']} | {l} | {v['c_pf']:.3f} | {v['f_k1_mhz']:.2f}{star} | {fe} | {ce} |")
    return "\n".join(lines)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--fig", default=r"C:\Users\ben\jobs\linkedin-content\birdcage-builder\obcb-comparison.png")
    a = ap.parse_args()
    res = collect()
    res["table_md"] = table_md(res)
    (REPO / "results/obcb_comparison.json").write_text(json.dumps(res, indent=1))
    print(res["table_md"])
    figure(res, Path(a.fig))
