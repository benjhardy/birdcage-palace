"""Write results/<name>.json for a tuned high-pass / band-pass case (same fields as results/lp16_3T.json, plus the
coarse-mesh tuning that seeded it and the order of the modes).

    python tools/export_types.py --case runs/hp16_3T --coarse runs/hp16_c1 --name hp16_3T
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True)
    ap.add_argument("--coarse", required=True)
    ap.add_argument("--name", required=True)
    a = ap.parse_args()
    case, coarse = Path(a.case), Path(a.coarse)
    res = json.loads((case / "modes.json").read_text())
    rc = json.loads((coarse / "modes.json").read_text())
    an_c = json.loads((coarse / "design.json").read_text())["notes"]
    des = res["design"]
    # the production case was built with the coarse-tuned C, so its notes hold that C; restore the ladder first guess
    des["notes"] = {"c_leg_first_guess": an_c.get("c_leg_first_guess"), "c_ring_first_guess": an_c.get("c_ring_first_guess"),
                    "built_with": "C from the coarse-mesh tune (" + coarse.name + "), then retuned on this mesh"}
    keys = ("m", "f_mhz", "label", "k", "k_purity", "roi_bt_cv", "roi_transverse_frac", "cap_free", "homogeneous",
            "ring_current_frac", "ring_bot_over_top")
    summ = {"design": des, "homogeneous_pair": res["homogeneous_pair"], "comparison": res["comparison"],
            "solve_stats": res["solve_stats"],
            "modes": [{k: m.get(k) for k in keys} for m in res["modes"]]}
    if (case / "tune_history.json").exists():
        summ["tune_history"] = json.loads((case / "tune_history.json").read_text())
    ch = json.loads((coarse / "tune_history.json").read_text())
    summ["first_guess"] = {"c_leg": ch[0]["c_leg"], "c_ring": ch[0]["c_ring"], "f_k1_mhz": ch[0]["f_k1_mhz"],
                           "rel_err": ch[0]["rel_err"], "mesh": "coarse (8/4/40 mm, order 1)", "case": coarse.name,
                           "note": "Palace at the ladder first-guess C was solved on the coarse mesh only"}
    summ["coarse_tune"] = {"case": coarse.name, "history": ch, "homogeneous_pair": rc["homogeneous_pair"],
                           "comparison": rc["comparison"]}
    lad = [m for m in res["modes"] if m.get("k") is not None and not m.get("cap_free") and m.get("k_purity", 0) > 0.5]
    summ["mode_order_by_frequency"] = [f"{m['f_mhz']:.2f} MHz {m['label']}" + (" (homogeneous)" if m["homogeneous"] else "")
                                       for m in sorted(res["modes"], key=lambda m: m["f_mhz"])]
    k1 = res["homogeneous_pair"]["f_mean_mhz"]
    # ring_bot_over_top < 0 (port sense; the bottom ports point the other way) = rings co-rotating: not a ladder mode
    ladder_k = [m for m in lad if m.get("ring_bot_over_top", 1) >= 0 and m.get("k", 0) > 0]
    corot = [m for m in res["modes"] if m.get("ring_bot_over_top", 1) < 0 and not m.get("cap_free")
             and m.get("ring_current_frac", 1) > 0.5]
    summ["k1_rank"] = {"ladder_k_modes_below_k1": sum(1 for m in ladder_k if m["f_mhz"] < k1 - 0.5),
                       "ladder_k_modes_above_k1": sum(1 for m in ladder_k if m["f_mhz"] > k1 + 0.5),
                       "corotating_ring_family_mhz": [round(m["f_mhz"], 2) for m in corot],
                       "corotating_ring_family_k": [m.get("k") for m in corot]}
    (ROOT / "results" / f"{a.name}.json").write_text(json.dumps(summ, indent=1))
    print("wrote results/" + a.name + ".json", summ["k1_rank"])


if __name__ == "__main__":
    main()
