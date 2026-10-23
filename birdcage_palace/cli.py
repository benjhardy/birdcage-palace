"""birdcage: design, build, solve and identify the modes of an MRI birdcage coil with Gmsh + Palace.

    birdcage build --legs 16 --diameter 0.28 --length 0.25 --type lowpass --b0 3.0 --shield 0.34 --out runs/lp16
    birdcage solve --case runs/lp16 --np 4
    birdcage modes --case runs/lp16            # modes.json + figures
    birdcage tune  --case runs/lp16 --tol 0.005  # rescale C by (f_sim/f_target)^2 and re-solve until k=1 is on target

A JSON spec (for example the spec.json from the web builder, https://benjhardy.github.io/birdcage-builder/) can
replace the build flags; flags given on the command line override it:

    birdcage build --spec spec.json --out runs/my_coil

Palace is called as `run_palace {cfg} {np}` by default. Use --palace-cmd or the BIRDCAGE_PALACE_CMD environment
variable for another launcher, e.g. "palace -np {np} {cfg}" or "apptainer run palace.sif -np {np} {cfg}".

Symmetry (low-pass): build a quarter model with PEC on x=0 and PMC on y=0 (keeps the k=1 member with B1 along y),
then reuse its mesh for another PEC/PMC combination:

    birdcage build ... --symmetry quarter --sym-bc x=pec,y=pmc --out runs/q_ep
    birdcage variant --case runs/q_ep --sym-bc x=pmc,y=pec --out runs/q_me
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

from .design import CoilDesign, Phantom
from . import ladder


def _design_from_args(a) -> CoilDesign:
    d = CoilDesign(legs=a.legs, diameter=a.diameter, length=a.length, leg_width=a.leg_width, ring_width=a.ring_width,
                   gap=a.gap, coil_type=a.type, b0=a.b0, shield_diameter=None if a.shield in (0, None) else a.shield,
                   shield_length=a.shield_length, shield_ends=a.shield_ends, c_leg=a.c_leg, c_ring=a.c_ring,
                   bandpass_ratio=a.bandpass_ratio, h_coil=a.h_coil, h_gap=a.h_gap, h_max=a.h_max, order=a.order,
                   symmetry=a.symmetry, sym_cut=a.sym_cut, sym_bc=_parse_bc(a.sym_bc), rung_offset=a.rung_offset,
                   copper_sigma=a.copper_sigma or None, esr=a.esr or None,
                   air_radius=a.air_radius, air_length=a.air_length)
    if a.phantom:
        d.phantom = Phantom(radius=a.phantom_radius, length=a.phantom_length, eps_r=a.phantom_eps, sigma=a.phantom_sigma,
                            loss_model=a.phantom_loss_model)
    return d


def cmd_build(a) -> None:
    from .geometry import build
    from .palace import write_config
    des = _design_from_args(a)
    des.validate()
    an = ladder.analytic(des)
    if des.c_leg is None and an["c_leg"] is not None:
        des.c_leg = an["c_leg"]
    if des.c_ring is None and an["c_ring"] is not None:
        des.c_ring = an["c_ring"]
    des.notes = {"c_leg_first_guess": an["c_leg"], "c_ring_first_guess": an["c_ring"]}
    an = ladder.analytic(des)  # spectrum with the capacitors actually used
    out = Path(a.out)
    t0 = time.time()
    info = build(des, out)
    (out / "analytic.json").write_text(json.dumps(an, indent=1))
    write_config(des, out, n_modes=a.modes)
    print(f"built {out}: {info['tets']} tets, {len(info['ports'])} capacitors, {time.time() - t0:.1f} s")
    print(f"f_target {des.f_larmor / 1e6:.2f} MHz; first-guess C_leg={_pf(des.c_leg)} C_ring={_pf(des.c_ring)}")
    print("ladder spectrum (MHz):", ", ".join(f"k{k}:{f / 1e6:.1f}" for k, f in zip(an["spectrum"]["k"], an["spectrum"]["f_hz"])))


def _parse_bc(txt: str | None) -> dict:
    bc = {"x": "pec", "y": "pmc", "z": "pec"}
    for item in (txt or "").split(","):
        if item.strip():
            k, v = item.split("=")
            bc[k.strip().lower()] = v.strip().lower()
    return bc


def cmd_variant(a) -> None:
    """Copy a built case (same mesh) with other symmetry-plane conditions, capacitor value, FE order or mode count."""
    import shutil
    from .palace import write_config
    src, out = Path(a.case), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    info = json.loads((src / "build.json").read_text())
    for f in (info["mesh"], "build.json", "surfaces.json", "analytic.json"):
        if (src / f).exists() and not (out / f).exists():
            shutil.copy(src / f, out / f)
    des = CoilDesign.from_json(src / "design.json")
    if a.sym_bc:
        des.sym_bc = _parse_bc(a.sym_bc)
    if a.c_leg:
        des.c_leg = a.c_leg
    if a.order:
        des.order = a.order
    if a.c_rungs:
        des.c_leg_rungs = [float(x) for x in json.loads(Path(a.c_rungs).read_text())]
    _apply_material_args(des, a)
    des.validate()
    des.to_json(out / "design.json")
    cfg = json.loads((src / "config.json").read_text())
    write_config(des, out, n_modes=a.modes or cfg["Solver"]["Eigenmode"]["N"], target_hz=cfg["Solver"]["Eigenmode"]["Target"] * 1e9,
                 h_plane=a.probe_h)
    ph = des.phantom
    print(f"variant {out}: sym_bc={ {k: des.sym_bc[k] for k in des.planes} } C_leg={_pf(des.c_leg)} order={des.order}"
          + (f" phantom eps_r={ph.eps_r:g} sigma={ph.sigma:g} ({ph.loss_model})" if ph else "")
          + (f" copper={des.copper_sigma:g} S/m" if des.copper_sigma else "") + (f" ESR={des.esr:g} ohm" if des.esr else ""))


def _apply_material_args(des: CoilDesign, a) -> None:
    """Phantom material and coil-loss options shared by `variant` and `load` (the phantom geometry is in the mesh)."""
    if any(getattr(a, k, None) is not None for k in ("phantom_eps", "phantom_sigma", "phantom_loss_model")):
        if des.phantom is None:
            raise SystemExit("this case has no phantom region in its mesh: build it with --phantom (e.g. --phantom-eps 1 "
                             "--phantom-sigma 0 for an air-filled region), then change the material here")
        if a.phantom_eps is not None:
            des.phantom.eps_r = a.phantom_eps
        if a.phantom_sigma is not None:
            des.phantom.sigma = a.phantom_sigma
        if a.phantom_loss_model is not None:
            des.phantom.loss_model = a.phantom_loss_model
    if getattr(a, "copper_sigma", None) is not None:
        des.copper_sigma = a.copper_sigma or None
    if getattr(a, "esr", None) is not None:
        des.esr = a.esr or None


def cmd_portcase(a) -> None:
    """Copy a built case (same mesh) as a port-tuning case: every capacitor gap a resistive port, driven sweep."""
    import shutil
    from .palace import write_driven_config
    src, out = Path(a.case), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    info = json.loads((src / "build.json").read_text())
    for f in (info["mesh"], "build.json", "surfaces.json", "analytic.json"):
        if (src / f).exists() and not (out / f).exists():
            shutil.copy(src / f, out / f)
    des = CoilDesign.from_json(src / "design.json")
    if a.sym_bc:
        des.sym_bc = _parse_bc(a.sym_bc)
    des.validate()
    des.to_json(out / "design.json")
    # default sweep: 100-160 MHz in 2 MHz steps at 3 T, scaled with the Larmor frequency for other fields
    sc = des.f_larmor / (3.0 * 42.577478e6)
    fmin, fmax, fstep = (round(d * sc, 3) if v is None else v for v, d in ((a.fmin, 100.0), (a.fmax, 160.0), (a.fstep, 2.0)))
    write_driven_config(des, out, fmin * 1e6, fmax * 1e6, fstep * 1e6, a.r_ref, a.adaptive_tol, a.max_samples)
    print(f"port case {out}: {len(info['ports'])} ports, R_ref {a.r_ref:g} ohm, {fmin:g}-{fmax:g} MHz in {fstep:g} MHz steps, "
          f"planes {des.sym_bc}")


def _pf(c):
    return "-" if c is None else f"{c * 1e12:.3f} pF"


DEFAULT_PALACE_CMD = "run_palace {cfg} {np}"


def cmd_solve(a) -> dict:
    from .palace import palace_busy, run
    np_ = a.np
    busy = palace_busy(a.wsl_distro)
    if busy and np_ > a.np_shared:
        print(f"another Palace job is running ({len(busy)} processes): using {a.np_shared} ranks")
        np_ = a.np_shared
    st = run(Path(a.case), np_, palace_cmd=a.palace_cmd, distro=a.wsl_distro)
    print("solve:", st)
    return st


def cmd_modes(a) -> dict:
    from .modes import analyse, compare
    case = Path(a.case)
    res = analyse(case)
    an = json.loads((case / "analytic.json").read_text())
    res["comparison"] = compare(an, res)
    res["analytic"] = an
    st = case / "solve_stats.json"
    res["solve_stats"] = json.loads(st.read_text()) if st.exists() else None
    (case / "modes.json").write_text(json.dumps(res, indent=1))
    print(f"{'m':>3} {'f (MHz)':>9} {'Q':>9}  {'label':<28} {'purity':>6} {'ROI std/mean':>12} {'transverse':>10}")
    for m in res["modes"]:
        q = m["Q"]
        print(f"{m['m']:>3} {m['f_mhz']:9.2f} {q if q < 1e6 else float('inf'):9.1f}  {m['label'] + (' <- homogeneous' if m['homogeneous'] else ''):<28}"
              f" {m.get('k_purity', 0):6.3f} {m.get('roi_bt_cv', 0):12.3f} {m.get('roi_transverse_frac', 0):10.3f}")
    hp = res.get("homogeneous_pair")
    if hp:
        print("homogeneous k=1:", json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in hp.items()}))
    if not a.no_figs:
        from . import plots
        plots.all_figures(case, Path(a.figdir) if a.figdir else case / "figs")
    return res


def cmd_tune(a) -> None:
    from .design import CoilDesign
    from .palace import write_config
    case = Path(a.case)
    hist = []
    for it in range(a.max_iter):
        des = CoilDesign.from_json(case / "design.json")
        if not (case / "postpro" / "eig.csv").exists() or it > 0 or a.resolve_first:
            cmd_solve(a)
        a.no_figs = True
        res = cmd_modes(a)
        f = res["homogeneous_pair"]["f_mean_mhz"] * 1e6
        err = f / des.f_larmor - 1
        hist.append({"iter": it, "c_leg": des.c_leg, "c_ring": des.c_ring, "f_k1_mhz": f / 1e6, "rel_err": err,
                     "solve": json.loads((case / "solve_stats.json").read_text())})
        print(f"tune iter {it}: k=1 at {f / 1e6:.3f} MHz ({err * 100:+.2f} %), C_leg={_pf(des.c_leg)} C_ring={_pf(des.c_ring)}",
              flush=True)
        (case / "tune_history.json").write_text(json.dumps(hist, indent=1))
        if abs(err) < a.tol:
            break
        # 1/f^2 = L_eff (C + C_stray) is linear in C: secant on the last two points, plain (f/f0)^2 scaling first.
        #    A slope measured on another (e.g. coarser) mesh can seed the first step: --slope from `birdcage slope`.
        x = [h["c_leg"] or h["c_ring"] for h in hist]
        y = [1 / h["f_k1_mhz"] ** 2 for h in hist]
        y0 = 1 / (des.f_larmor / 1e6) ** 2
        slope = (y[-1] - y[-2]) / (x[-1] - x[-2]) if len(x) >= 2 and x[-1] != x[-2] else a.slope
        s = (x[-1] + (y0 - y[-1]) / slope) / x[-1] if slope else (f / des.f_larmor) ** 2
        hist[-1]["next_scale"] = s
        des.c_leg = des.c_leg * s if des.c_leg else None
        des.c_ring = des.c_ring * s if des.c_ring else None
        des.to_json(case / "design.json")
        an = ladder.analytic(des)
        (case / "analytic.json").write_text(json.dumps(an, indent=1))
        write_config(des, case, n_modes=json.loads((case / "config.json").read_text())["Solver"]["Eigenmode"]["N"])
    a.no_figs = False
    res = cmd_modes(a)
    summ = tune_summary(CoilDesign.from_json(case / "design.json"), res, hist)
    (case / "tune_summary.json").write_text(json.dumps(summ, indent=1))
    print("tune summary:", json.dumps(summ), flush=True)


def tune_summary(des: CoilDesign, res: dict, hist: list) -> dict:
    """Short record of a tune run: k=1 frequency, capacitors, B1 homogeneity, number of solves."""
    hp = res.get("homogeneous_pair") or {}
    f = hp.get("f_mean_mhz")
    first = (des.notes or {})
    return {"f_target_mhz": des.f_larmor / 1e6, "f_k1_mhz": f,
            "rel_err": None if f is None else f / (des.f_larmor / 1e6) - 1,
            "k1_split_mhz": hp.get("split_mhz"), "c_leg_pf": None if des.c_leg is None else des.c_leg * 1e12,
            "c_ring_pf": None if des.c_ring is None else des.c_ring * 1e12,
            "c_leg_first_guess_pf": None if first.get("c_leg_first_guess") is None else first["c_leg_first_guess"] * 1e12,
            "c_ring_first_guess_pf": None if first.get("c_ring_first_guess") is None else first["c_ring_first_guess"] * 1e12,
            "f_k1_at_first_guess_mhz": hist[0]["f_k1_mhz"] if hist else None,
            "roi_bt_cv": hp.get("roi_bt_cv"), "b1plus_roi_cv": hp.get("b1plus_roi_cv"),
            "b1minus_over_b1plus_roi": hp.get("b1minus_over_b1plus_roi"), "solves": len(hist),
            "palace_wall_s": [h["solve"].get("wall_s") for h in hist]}


def _sub(a, **kw):
    b = argparse.Namespace(**vars(a))
    for k, v in kw.items():
        setattr(b, k, v)
    return b


def _k1_q(res: dict | None) -> float | None:
    if not res:
        return None
    hp = res.get("homogeneous_pair") or {}
    qs = [m["Q"] for m in res["modes"] if m["m"] in hp.get("m", [])]
    return float(np.mean(qs)) if qs else None


def _single_mode_q(case: Path) -> tuple[float, float]:
    """f (MHz) and Q of the one mode of a coil-loss solve (first row of eig.csv)."""
    from .modes import _read_csv
    _, d = _read_csv(case / "postpro" / "eig.csv")
    return float(d[0, 1] * 1e3), float(d[0, 3])


def _loss_variant(src: Path, out: Path, des: CoilDesign, n_modes: int, target_hz: float | None) -> None:
    import shutil
    from .palace import write_config
    out.mkdir(parents=True, exist_ok=True)
    info = json.loads((src / "build.json").read_text())
    for f in (info["mesh"], "build.json", "surfaces.json", "analytic.json"):
        if (src / f).exists() and not (out / f).exists():
            shutil.copy(src / f, out / f)
    des.validate()
    des.to_json(out / "design.json")
    write_config(des, out, n_modes=n_modes, target_hz=target_hz)


def cmd_load(a) -> dict:
    """Fill the phantom region of a tuned case with a load and solve it on the same mesh (C not retuned).

    1. loaded, lossless coil: k = 1 frequency and shift, phantom-only Q, quadrature B1+ (modes.json in <out>/loaded)
    2. with --esr / --copper-sigma: the coil losses added, loaded and unloaded (one mode near f0 each): Q_L, Q_U, Q_U/Q_L
    Writes <out>/load_summary.json and adds it as "loaded" to the case's tune_summary.json.
    """
    case = Path(a.case)
    out = Path(a.out) if a.out else case.parent / (case.name + "_loaded")
    des0 = CoilDesign.from_json(case / "design.json")
    if des0.phantom is None:
        raise SystemExit("load needs a phantom region in the mesh: build with --phantom --phantom-eps 1 --phantom-sigma 0")
    cfg = json.loads((case / "config.json").read_text())
    n_modes = a.modes or min(cfg["Solver"]["Eigenmode"]["N"], 4)
    target = cfg["Solver"]["Eigenmode"]["Target"] * 1e9
    unl = json.loads((case / "modes.json").read_text()) if (case / "modes.json").exists() else None
    des = CoilDesign.from_json(case / "design.json")
    des.phantom.eps_r, des.phantom.sigma = a.phantom_eps, a.phantom_sigma
    des.copper_sigma = des.esr = None
    runs = {}
    _loss_variant(case, out / "loaded", des, n_modes, target)
    cmd_solve(_sub(a, case=str(out / "loaded")))
    res = cmd_modes(_sub(a, case=str(out / "loaded"), no_figs=True, figdir=None))
    hp = res.get("homogeneous_pair") or {}
    runs["loaded"] = json.loads((out / "loaded" / "solve_stats.json").read_text())
    f_u = ((unl or {}).get("homogeneous_pair") or {}).get("f_mean_mhz")
    f_l = hp.get("f_mean_mhz")
    summ = {"phantom": {"radius_m": des.phantom.radius, "length_m": des.phantom.length, "eps_r": a.phantom_eps,
                        "sigma_s_per_m": a.phantom_sigma},
            "c_leg_pf": None if des.c_leg is None else des.c_leg * 1e12,
            "c_ring_pf": None if des.c_ring is None else des.c_ring * 1e12,
            "f_k1_unloaded_mhz": f_u, "f_k1_loaded_mhz": f_l,
            "f_shift_rel": f_l / f_u - 1 if (f_u and f_l) else None,
            "Q_phantom": _k1_q(res), "Q_unloaded_lossless_coil": _k1_q(unl),
            "b1plus_roi_cv": hp.get("b1plus_roi_cv"), "b1minus_over_b1plus_roi": hp.get("b1minus_over_b1plus_roi"),
            "b1plus_roi_cv_unloaded": ((unl or {}).get("homogeneous_pair") or {}).get("b1plus_roi_cv")}
    if a.esr or a.copper_sigma:
        # one mode near f0: the nonlinear (copper) eigensolve can stall on strongly damped phantom modes
        for tag, eps, sig in (("loaded_coil_loss", a.phantom_eps, a.phantom_sigma), ("unloaded_coil_loss", 1.0, 0.0)):
            d = CoilDesign.from_json(case / "design.json")
            d.phantom.eps_r, d.phantom.sigma = eps, sig
            d.esr, d.copper_sigma = a.esr or None, a.copper_sigma or None
            _loss_variant(case, out / tag, d, 1, 0.9 * d.f_larmor)
            cmd_solve(_sub(a, case=str(out / tag)))
            runs[tag] = json.loads((out / tag / "solve_stats.json").read_text())
        fl, ql = _single_mode_q(out / "loaded_coil_loss")
        fu, qu = _single_mode_q(out / "unloaded_coil_loss")
        summ.update(coil_loss={"esr_ohm": a.esr, "copper_sigma_s_per_m": a.copper_sigma,
                               "note": "copper on strips and shield as a Palace Conductivity boundary (one face of each "
                                       "strip sheet: an upper bound on the strip loss); ESR as the parallel R at f0"},
                    Q_U=qu, Q_L=ql, QU_over_QL=qu / ql, f_k1_loaded_lossy_mhz=fl, f_k1_unloaded_lossy_mhz=fu)
    summ["solves"] = runs
    (out / "load_summary.json").write_text(json.dumps(summ, indent=1))
    ts = case / "tune_summary.json"
    if ts.exists():
        t = json.loads(ts.read_text())
        t["loaded"] = summ
        ts.write_text(json.dumps(t, indent=1))
    print("load summary:", json.dumps({k: v for k, v in summ.items() if k != "solves"}), flush=True)
    return summ


def _f_at_c(hist: list, c: float) -> float | None:
    """k=1 frequency at capacitor c from a tune history (1/f^2 linear in C; with one point, 1/f^2 proportional to C)."""
    x = [h["c_leg"] or h["c_ring"] for h in hist]
    y = [1 / h["f_k1_mhz"] ** 2 for h in hist]
    if not x:
        return None
    if len(set(x)) >= 2:
        k, b = np.polyfit(x, y, 1)
        return float((k * c + b) ** -0.5)
    return float((y[0] * c / x[0]) ** -0.5)


def cmd_compare(a) -> dict:
    """Compare tuned cases (e.g. with and without the shield): k=1 shift, tuned C, Q, B1 homogeneity."""
    rows = []
    for c in a.cases:
        case = Path(c)
        des = CoilDesign.from_json(case / "design.json")
        ts = json.loads((case / "tune_summary.json").read_text())
        hist = json.loads((case / "tune_history.json").read_text())
        res = json.loads((case / "modes.json").read_text())
        q = _k1_q(res)
        rows.append({"case": str(case), "shield": "none (absorbing boundary)" if not des.shield_diameter else
                     f"{des.shield_diameter * 100:g} cm, ends {des.shield_ends}",
                     "c_first_guess_pf": ts.get("c_leg_first_guess_pf") or ts.get("c_ring_first_guess_pf"),
                     "f_k1_at_first_guess_mhz": ts.get("f_k1_at_first_guess_mhz"),
                     "c_tuned_pf": ts.get("c_leg_pf") or ts.get("c_ring_pf"), "f_k1_tuned_mhz": ts.get("f_k1_mhz"),
                     "Q_k1": None if q is None or q > 1e9 else q,
                     "Q_note": "radiation Q (absorbing outer boundary)" if not des.shield_diameter else
                               "lossless conductors in a closed shield: no loss, Q infinite",
                     "b1plus_roi_cv": ts.get("b1plus_roi_cv"), "roi_bt_cv": ts.get("roi_bt_cv"),
                     "b1minus_over_b1plus_roi": ts.get("b1minus_over_b1plus_roi"), "solves": ts.get("solves"),
                     "palace_wall_s": ts.get("palace_wall_s"), "_hist": hist, "loaded": ts.get("loaded")})
    out = {"what": "birdcage compare (simulation only)", "cases": [{k: v for k, v in r.items() if k != "_hist"} for r in rows]}
    if len(rows) == 2:
        a0, b0 = rows
        cref = a0["c_tuned_pf"] * 1e-12
        fb = _f_at_c(b0["_hist"], cref)
        out["k1_shift"] = {
            "what": f"k=1 of case 2 at case 1's tuned C ({a0['c_tuned_pf']:.4f} pF), from case 2's 1/f^2 vs C fit",
            "f_case1_mhz": a0["f_k1_tuned_mhz"], "f_case2_at_same_c_mhz": fb,
            "shift_rel": None if fb is None else fb / a0["f_k1_tuned_mhz"] - 1}
        out["c_tuned_ratio_case2_over_case1"] = b0["c_tuned_pf"] / a0["c_tuned_pf"]
    def fmt(v, f):
        if v is None:
            return "-"
        return format(v, f) if f else str(v)
    lines = ["| | " + " | ".join(Path(r["case"]).name for r in rows) + " |", "|---|" + "---|" * len(rows)]
    for lab, key, f in (("shield", "shield", ""), ("C first guess (pF)", "c_first_guess_pf", ".3f"),
                        ("k=1 at first guess (MHz)", "f_k1_at_first_guess_mhz", ".2f"), ("tuned C (pF)", "c_tuned_pf", ".4f"),
                        ("k=1 tuned (MHz)", "f_k1_tuned_mhz", ".3f"), ("Q of k=1 (Palace)", "Q_k1", ".0f"),
                        ("quadrature B1+ ROI std/mean", "b1plus_roi_cv", ".2%"),
                        ("B1- / B1+ in the ROI", "b1minus_over_b1plus_roi", ".2%"), ("solves", "solves", "d")):
        lines.append(f"| {lab} | " + " | ".join(fmt(r[key], f) for r in rows) + " |")
    ks = out.get("k1_shift")
    if ks and ks["shift_rel"] is not None:
        lines.append(f"\nk=1 shift at the same C ({rows[0]['c_tuned_pf']:.4f} pF): {ks['f_case1_mhz']:.3f} -> "
                     f"{ks['f_case2_at_same_c_mhz']:.3f} MHz ({ks['shift_rel'] * 100:+.2f} %).")
    out["table_md"] = "\n".join(lines)
    o = Path(a.out)
    o.parent.mkdir(parents=True, exist_ok=True)
    o.write_text(json.dumps(out, indent=1))
    o.with_suffix(".md").write_text("# Comparison (simulation only)\n\n" + out["table_md"] + "\n")
    print(out["table_md"], flush=True)
    return out


def load_spec(path) -> dict:
    """Build options from a JSON spec: either flat, or the web builder's {"build": {...}} form. Keys are CLI option
    names with - or _ (leg-width or leg_width)."""
    d = json.loads(Path(path).read_text())
    if isinstance(d.get("build"), dict):
        d = d["build"]
    return {k.replace("-", "_"): v for k, v in d.items()}


def parse_args(argv=None) -> argparse.Namespace:
    ap = _parser()
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "build":
        pre = argparse.ArgumentParser(add_help=False)
        pre.add_argument("--spec")
        known, _ = pre.parse_known_args(argv[1:])
        if known.spec:
            spec = load_spec(known.spec)
            b = ap._build_parser
            dests = {act.dest for act in b._actions} - {"help", "spec", "fn"}
            bad = sorted(set(spec) - dests)
            if bad:
                ap.error(f"unknown key(s) in {known.spec}: {', '.join(bad)}")
            for act in b._actions:  # set_defaults skips choices= checks, so check them here
                if act.choices and act.dest in spec and spec[act.dest] not in act.choices:
                    ap.error(f"{known.spec}: {act.dest} must be one of {', '.join(map(str, act.choices))}")
            b.set_defaults(**spec)
    a = ap.parse_args(argv)
    if a.cmd == "build" and not a.out:
        ap.error("build needs --out (or \"out\" in the spec)")
    return a


def main(argv=None) -> None:
    a = parse_args(argv)
    a.fn(a)


def _material_options(p) -> None:
    p.add_argument("--phantom-eps", type=float, default=None, help="phantom eps_r (the case must have a phantom region)")
    p.add_argument("--phantom-sigma", type=float, default=None, help="phantom conductivity, S/m")
    p.add_argument("--phantom-loss-model", choices=["conductivity", "losstan"], default=None)
    p.add_argument("--copper-sigma", type=float, default=None, help="S/m on strips and shield; 0 = PEC")
    p.add_argument("--esr", type=float, default=None, help="ohm per capacitor; 0 = ideal")


def _parser() -> argparse.ArgumentParser:
    import os
    from . import __version__
    ap = argparse.ArgumentParser(prog="birdcage", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-V", "--version", action="version", version=f"birdcage-palace {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build", help="analytic first guess + Gmsh mesh + Palace config")
    ap._build_parser = b
    b.add_argument("--spec", default=None, help="JSON file with build options (keys = option names); flags override it")
    b.add_argument("--legs", type=int, default=16)
    b.add_argument("--diameter", type=float, default=0.28)
    b.add_argument("--length", type=float, default=0.25, help="rung length between end-ring centrelines, m")
    b.add_argument("--leg-width", type=float, default=0.012)
    b.add_argument("--ring-width", type=float, default=0.015)
    b.add_argument("--gap", type=float, default=0.005, help="capacitor gap, m")
    b.add_argument("--type", choices=["lowpass", "highpass", "bandpass"], default="lowpass")
    b.add_argument("--b0", type=float, default=3.0, help="field strength, T (target = 42.577 MHz/T x B0)")
    b.add_argument("--shield", type=float, default=0.34, help="shield diameter, m (0 = unshielded, absorbing boundary)")
    b.add_argument("--shield-length", type=float, default=None)
    b.add_argument("--shield-ends", choices=["pec", "absorbing"], default="pec")
    b.add_argument("--air-radius", type=float, default=None,
                   help="unshielded: radius of the air cylinder, m (default R + max(0.15 m, lambda0/4))")
    b.add_argument("--air-length", type=float, default=None,
                   help="unshielded: length of the air cylinder, m (default L + 2 max(0.15 m, lambda0/4))")
    b.add_argument("--c-leg", type=float, default=None, help="F; default from the ladder model")
    b.add_argument("--c-ring", type=float, default=None, help="F; default from the ladder model")
    b.add_argument("--bandpass-ratio", type=float, default=1.0)
    b.add_argument("--phantom", action="store_true")
    b.add_argument("--phantom-radius", type=float, default=0.08)
    b.add_argument("--phantom-length", type=float, default=0.20)
    b.add_argument("--phantom-eps", type=float, default=78.0)
    b.add_argument("--phantom-sigma", type=float, default=0.6)
    b.add_argument("--phantom-loss-model", choices=["conductivity", "losstan"], default="conductivity",
                   help="phantom loss as a conductivity, or as the loss tangent sigma/(w0 eps0 eps_r) (linear eigenproblem)")
    b.add_argument("--copper-sigma", type=float, default=None, help="S/m: copper surface loss on strips and shield (default PEC)")
    b.add_argument("--esr", type=float, default=None, help="ohm: series resistance of every capacitor (default ideal)")
    b.add_argument("--h-coil", type=float, default=0.005)
    b.add_argument("--h-gap", type=float, default=0.0025)
    b.add_argument("--h-max", type=float, default=0.035)
    b.add_argument("--order", type=int, default=2)
    b.add_argument("--modes", type=int, default=20, help="number of eigenmodes to request")
    b.add_argument("--symmetry", choices=["none", "half", "quarter", "eighth", "zhalf"], default="none",
                   help="model 1/2 (x>=0), 1/4 (x,y>=0), 1/8 (x,y,z>=0) or the z>=0 half (all rungs) of a low-pass coil")
    b.add_argument("--sym-cut", choices=["between", "through"], default="between",
                   help="symmetry planes between rungs (cage rotated by pi/N, default) or through rung centres")
    b.add_argument("--sym-bc", default="x=pec,y=pmc,z=pec",
                   help="PEC or PMC on each plane, e.g. x=pec,y=pmc (z=0 must be PEC)")
    b.add_argument("--rung-offset", type=float, default=None, help="angle of rung 0, rad (default 0, or pi/N with --sym-cut between)")
    b.add_argument("--out", default=None, help="case directory (required, here or in the spec)")
    b.set_defaults(fn=cmd_build)
    v = sub.add_parser("variant", help="copy a built case (same mesh) with other plane conditions / C / order / modes")
    v.add_argument("--case", required=True)
    v.add_argument("--out", required=True)
    v.add_argument("--sym-bc", default=None)
    v.add_argument("--c-leg", type=float, default=None)
    v.add_argument("--order", type=int, default=None)
    v.add_argument("--modes", type=int, default=None)
    v.add_argument("--c-rungs", default=None, help="JSON file with one capacitor value (F) per rung")
    _material_options(v)
    v.add_argument("--probe-h", type=float, default=0.006, help="probe spacing on the axial and coronal planes, m")
    v.set_defaults(fn=cmd_variant)
    pc = sub.add_parser("portcase", help="copy a built case as a driven port-tuning case (gaps -> resistive ports)")
    pc.add_argument("--case", required=True)
    pc.add_argument("--out", required=True)
    pc.add_argument("--sym-bc", default=None)
    pc.add_argument("--r-ref", type=float, default=500.0, help="port reference resistance, ohm")
    pc.add_argument("--fmin", type=float, default=None, help="MHz (default 100 MHz x B0 / 3 T)")
    pc.add_argument("--fmax", type=float, default=None, help="MHz (default 160 MHz x B0 / 3 T)")
    pc.add_argument("--fstep", type=float, default=None,
                    help="MHz (default 2 MHz x B0 / 3 T; the adaptive ROM is evaluated, and post-processed, at every step)")
    pc.add_argument("--adaptive-tol", type=float, default=1e-5)
    pc.add_argument("--max-samples", type=int, default=40)
    pc.set_defaults(fn=cmd_portcase)
    from . import porttune
    porttune.add_parser(sub)
    for name, fn, hlp in (("solve", cmd_solve, "run Palace"), ("modes", cmd_modes, "identify modes, write modes.json + figures"),
                          ("tune", cmd_tune, "retune C until the k=1 mode is on the Larmor frequency")):
        s = sub.add_parser(name, help=hlp)
        s.add_argument("--case", required=True)
        s.add_argument("--np", type=int, default=4)
        s.add_argument("--np-shared", type=int, default=4, help="ranks to use when another Palace job is running")
        s.add_argument("--figdir", default=None)
        s.add_argument("--no-figs", action="store_true")
        s.add_argument("--tol", type=float, default=0.005)
        s.add_argument("--max-iter", type=int, default=4)
        s.add_argument("--resolve-first", action="store_true")
        s.add_argument("--slope", type=float, default=None,
                       help="d(1/f_MHz^2)/dC (1/(MHz^2 F)) from an earlier pair of solves, seeds the first retune step")
        s.add_argument("--palace-cmd", default=os.environ.get("BIRDCAGE_PALACE_CMD", DEFAULT_PALACE_CMD),
                       help="Palace launcher with {cfg} and {np} placeholders (env BIRDCAGE_PALACE_CMD); "
                            f"default \"{DEFAULT_PALACE_CMD}\"")
        s.add_argument("--wsl-distro", default=os.environ.get("BIRDCAGE_WSL_DISTRO", "auto"),
                       help="Windows only: WSL distribution that runs Palace (env BIRDCAGE_WSL_DISTRO; default auto = the "
                            "WSL default distribution). On Linux Palace runs natively.")
        s.set_defaults(fn=fn)
        if name == "tune":
            tune_p = s
    ld = sub.add_parser("load", help="fill the phantom region of a tuned case with a load (same mesh, C not retuned)")
    for act in tune_p._actions:  # same solver options as tune
        if act.dest not in ("help", "fn") and act.option_strings:
            ld._add_action(act)
    ld.add_argument("--out", default=None, help="default <case>_loaded")
    ld.add_argument("--phantom-eps", type=float, default=78.0)
    ld.add_argument("--phantom-sigma", type=float, default=0.6)
    ld.add_argument("--esr", type=float, default=0.68, help="ohm per capacitor for Q_U and Q_L (0 = skip the coil-loss solves)")
    ld.add_argument("--copper-sigma", type=float, default=5.8e7, help="S/m on strips and shield for Q_U and Q_L (0 = PEC)")
    ld.add_argument("--modes", type=int, default=None)
    ld.set_defaults(fn=cmd_load)
    cp = sub.add_parser("compare", help="compare tuned cases (e.g. with and without the shield)")
    cp.add_argument("--cases", nargs="+", required=True)
    cp.add_argument("--out", default="compare.json")
    cp.set_defaults(fn=cmd_compare)
    return ap


if __name__ == "__main__":
    main(sys.argv[1:])
