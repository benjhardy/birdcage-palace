"""Palace eigen solves at the capacitor values from Open BirdcageBuilder, our ladder model and the Palace tune.

Each case reuses the mesh of a validated reference case unchanged (same mesh, order, solver settings, mode count and
shift target); only the lumped capacitor values in config.json change.

    low-pass  : eighth-symmetry model at production density (runs/sym_e_m2o2_clean, 0.156 M unknowns)
    high-pass : full model, production mesh (runs/hp16_3T) - the symmetry planes are implemented for low-pass only
    band-pass : full model, production mesh (runs/bp16_3T)

    python -u compare/open-birdcagebuilder/palace_checks.py lp_obcb lp_tuned ...   (names from CASES below)

Results are appended to compare/open-birdcagebuilder/palace_checks.json (one entry per solve).
"""
from __future__ import annotations

import json
import shutil
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from birdcage_palace.design import CoilDesign  # noqa: E402
from birdcage_palace.modes import analyse  # noqa: E402
from birdcage_palace.palace import palace_busy, run, write_config  # noqa: E402

OUT = Path(__file__).with_name("palace_checks.json")
OBCB = json.loads(Path(__file__).with_name("obcb_raw.json").read_text())


def obcb_pf(name: str) -> float:
    return next(c for c in OBCB["cases"] if c["name"] == name)["result"]["cap_pf"]


def design_c(case: str, key: str) -> float:
    return json.loads((REPO / "results" / case).read_text())["design"][key] * 1e12


LP_LADDER = json.loads((REPO / "results/lp16_3T.json").read_text())["design"]["notes"]["c_leg_first_guess"] * 1e12
HP_LADDER = json.loads((REPO / "results/hp16_3T.json").read_text())["design"]["notes"]["c_ring_first_guess"] * 1e12
BP_LADDER = json.loads((REPO / "results/bp16_3T.json").read_text())["design"]["notes"]["c_leg_first_guess"] * 1e12

# name: (source case, coil type, C in pF, what)
CASES = {
    "lp_tuned": ("sym_e_m2o2_clean", "lowpass", design_c("lp16_3T.json", "c_leg"), "Palace-tuned C (baseline)"),
    "lp_obcb": ("sym_e_m2o2_clean", "lowpass", obcb_pf("lp_shield34_L25"), "Open BirdcageBuilder"),
    "lp_ladder": ("sym_e_m2o2_clean", "lowpass", LP_LADDER, "our ladder first guess"),
    "hp_obcb": ("hp16_3T", "highpass", obcb_pf("hp_shield34_L25"), "Open BirdcageBuilder"),
    "bp_obcb": ("bp16_3T", "bandpass", obcb_pf("bp_leg_shield34_L25"), "Open BirdcageBuilder"),
    "hp_ladder": ("hp16_3T", "highpass", HP_LADDER, "our ladder first guess"),
    "bp_ladder": ("bp16_3T", "bandpass", BP_LADDER, "our ladder first guess"),
}


def make_case(name: str) -> Path:
    src_name, ctype, c_pf, _ = CASES[name]
    src, out = REPO / "runs" / src_name, REPO / "runs" / f"obcb_{name}"
    out.mkdir(parents=True, exist_ok=True)
    info = json.loads((src / "build.json").read_text())
    for f in (info["mesh"], "build.json", "surfaces.json", "analytic.json"):
        if (src / f).exists():
            shutil.copy(src / f, out / f)
    des = CoilDesign.from_json(src / "design.json")
    assert des.coil_type == ctype
    c = c_pf * 1e-12
    if ctype in ("lowpass", "bandpass"):
        des.c_leg = c
    if ctype in ("highpass", "bandpass"):
        des.c_ring = c
    des.to_json(out / "design.json")
    cfg = json.loads((src / "config.json").read_text())
    ev = cfg["Solver"]["Eigenmode"]
    new = write_config(des, out, n_modes=ev["N"], target_hz=ev["Target"] * 1e9)
    # the solver block must be the reference one, apart from the capacitor values
    assert new["Solver"] == cfg["Solver"], (new["Solver"], cfg["Solver"])
    return out


def main(names: list[str]) -> None:
    results = json.loads(OUT.read_text()) if OUT.exists() else []
    for name in names:
        busy = palace_busy()
        if busy:
            raise SystemExit(f"another Palace job is running: {busy}")
        case = make_case(name)
        t0 = time.time()
        st = run(case, 4)
        res = analyse(case)
        hp = res.get("homogeneous_pair") or {}
        k1 = [m for m in res["modes"] if m["label"] == "k=1"]
        entry = {"name": name, "what": CASES[name][3], "coil_type": CASES[name][1], "source_case": CASES[name][0],
                 "c_pf": CASES[name][2], "case": f"runs/obcb_{name}", "f_k1_mhz": hp.get("f_mean_mhz"),
                 "f_k1_members_mhz": hp.get("f_mhz"), "k1_modes_found_mhz": [m["f_mhz"] for m in k1],
                 "roi_bt_cv": hp.get("roi_bt_cv"), "solve": st, "wall_total_s": round(time.time() - t0, 1),
                 "date": time.strftime("%Y-%m-%d %H:%M")}
        results = [r for r in results if r["name"] != name] + [entry]
        OUT.write_text(json.dumps(results, indent=1))
        print(json.dumps({k: entry[k] for k in ("name", "c_pf", "f_k1_mhz", "k1_modes_found_mhz")}), st, flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
