"""birdcage: design, build, solve and identify the modes of an MRI birdcage coil with Gmsh + Palace.

    birdcage build --legs 16 --diameter 0.28 --length 0.25 --type lowpass --b0 3.0 --shield 0.34 --out runs/lp16
    birdcage solve --case runs/lp16 --np 4
    birdcage modes --case runs/lp16            # modes.json + figures
    birdcage tune  --case runs/lp16 --tol 0.005  # rescale C by (f_sim/f_target)^2 and re-solve until k=1 is on target

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

from .design import CoilDesign, Phantom
from . import ladder


def _design_from_args(a) -> CoilDesign:
    d = CoilDesign(legs=a.legs, diameter=a.diameter, length=a.length, leg_width=a.leg_width, ring_width=a.ring_width,
                   gap=a.gap, coil_type=a.type, b0=a.b0, shield_diameter=None if a.shield in (0, None) else a.shield,
                   shield_length=a.shield_length, shield_ends=a.shield_ends, c_leg=a.c_leg, c_ring=a.c_ring,
                   bandpass_ratio=a.bandpass_ratio, h_coil=a.h_coil, h_gap=a.h_gap, h_max=a.h_max, order=a.order,
                   symmetry=a.symmetry, sym_cut=a.sym_cut, sym_bc=_parse_bc(a.sym_bc), rung_offset=a.rung_offset)
    if a.phantom:
        d.phantom = Phantom(radius=a.phantom_radius, length=a.phantom_length, eps_r=a.phantom_eps, sigma=a.phantom_sigma)
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
    des.validate()
    des.to_json(out / "design.json")
    cfg = json.loads((src / "config.json").read_text())
    write_config(des, out, n_modes=a.modes or cfg["Solver"]["Eigenmode"]["N"], target_hz=cfg["Solver"]["Eigenmode"]["Target"] * 1e9)
    print(f"variant {out}: sym_bc={ {k: des.sym_bc[k] for k in des.planes} } C_leg={_pf(des.c_leg)} order={des.order}")


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
    write_driven_config(des, out, a.fmin * 1e6, a.fmax * 1e6, a.fstep * 1e6, a.r_ref, a.adaptive_tol, a.max_samples)
    print(f"port case {out}: {len(info['ports'])} ports, R_ref {a.r_ref:g} ohm, {a.fmin:g}-{a.fmax:g} MHz")


def _pf(c):
    return "-" if c is None else f"{c * 1e12:.3f} pF"


def cmd_solve(a) -> dict:
    from .palace import palace_busy, run
    np_ = a.np
    busy = palace_busy()
    if busy and np_ > a.np_shared:
        print(f"another Palace job is running ({len(busy)} processes): using {a.np_shared} ranks")
        np_ = a.np_shared
    st = run(Path(a.case), np_)
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
        print(f"tune iter {it}: k=1 at {f / 1e6:.3f} MHz ({err * 100:+.2f} %), C_leg={_pf(des.c_leg)} C_ring={_pf(des.c_ring)}")
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
    cmd_modes(a)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(prog="birdcage", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build", help="analytic first guess + Gmsh mesh + Palace config")
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
    b.add_argument("--c-leg", type=float, default=None, help="F; default from the ladder model")
    b.add_argument("--c-ring", type=float, default=None, help="F; default from the ladder model")
    b.add_argument("--bandpass-ratio", type=float, default=1.0)
    b.add_argument("--phantom", action="store_true")
    b.add_argument("--phantom-radius", type=float, default=0.08)
    b.add_argument("--phantom-length", type=float, default=0.20)
    b.add_argument("--phantom-eps", type=float, default=78.0)
    b.add_argument("--phantom-sigma", type=float, default=0.6)
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
    b.add_argument("--out", required=True)
    b.set_defaults(fn=cmd_build)
    v = sub.add_parser("variant", help="copy a built case (same mesh) with other plane conditions / C / order / modes")
    v.add_argument("--case", required=True)
    v.add_argument("--out", required=True)
    v.add_argument("--sym-bc", default=None)
    v.add_argument("--c-leg", type=float, default=None)
    v.add_argument("--order", type=int, default=None)
    v.add_argument("--modes", type=int, default=None)
    v.add_argument("--c-rungs", default=None, help="JSON file with one capacitor value (F) per rung")
    v.set_defaults(fn=cmd_variant)
    pc = sub.add_parser("portcase", help="copy a built case as a driven port-tuning case (gaps -> resistive ports)")
    pc.add_argument("--case", required=True)
    pc.add_argument("--out", required=True)
    pc.add_argument("--sym-bc", default=None)
    pc.add_argument("--r-ref", type=float, default=500.0, help="port reference resistance, ohm")
    pc.add_argument("--fmin", type=float, default=100.0, help="MHz")
    pc.add_argument("--fmax", type=float, default=160.0, help="MHz")
    pc.add_argument("--fstep", type=float, default=2.0, help="MHz (the adaptive ROM is evaluated, and post-processed, at every step)")
    pc.add_argument("--adaptive-tol", type=float, default=1e-5)
    pc.add_argument("--max-samples", type=int, default=40)
    pc.set_defaults(fn=cmd_portcase)
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
        s.set_defaults(fn=fn)
    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main(sys.argv[1:])
