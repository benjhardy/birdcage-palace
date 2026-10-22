"""Phantom loading study: the tuned 16-rung low-pass coil (C = 1.828 pF, NOT retuned) with a centred cylindrical
phantom, swept in conductivity and in permittivity. Plan and cost: docs/loading-study.md.

Every case is three Palace eigen solves on the eighth model (x = 0, y = 0, z = 0 planes; production mesh density,
order 2), all on one mesh per phantom radius:
  em   x=0 PEC, y=0 PMC, z=0 PEC: the k = 1 member with B1 along y; lossless coil, so its Q is the phantom Q
  me   x=0 PMC, y=0 PEC, z=0 PEC: the partner (B1 along x); with em it gives the quadrature B1+ map
  full em with the coil losses added (copper strips and shield, Palace Conductivity boundary; capacitor ESR 0.68 ohm
       as the equivalent parallel R): its Q is the loaded Q, Q_L. The case eps_r = 1, sigma = 0 gives Q_U.

    python tools/loading_study.py setup          # meshes + case folders (runs/loading/...)
    python tools/loading_study.py run [--only r8] # solve every unsolved case in order (Palace in WSL, foreground)
    python tools/loading_study.py analyse        # results/loading_sweep.json

Simulation only.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from birdcage_palace import modes as M  # noqa: E402
from birdcage_palace.design import CoilDesign  # noqa: E402

RUNS = ROOT / "runs" / "loading"
C_TUNED = 1.828e-12
ESR, CU = 0.68, 5.8e7
PROBE_H = 0.004
SIGMAS = [0, 0.05, 0.1, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0]   # S/m at eps_r 78
EPSS = [1, 5, 10, 20, 30, 40, 50, 60, 70, 78]                           # at sigma 0
RADII = [0.05, 0.065, 0.08, 0.095, 0.11]                                # m, at (78, 0.6) and (1, 0)
LENGTH = 0.20


def cases() -> list[dict]:
    """In run order: the sigma sweep, the eps_r sweep, conductivity only, then the other radii."""
    out, seen = [], set()
    def add(r, eps, sig, sweep):
        key = (r, eps, sig)
        if key in seen:
            for c in out:
                if (c["r"], c["eps"], c["sigma"]) == key:
                    c["sweeps"].append(sweep)
            return
        seen.add(key)
        out.append({"name": f"r{round(r * 1000):03d}_e{eps:g}_s{sig:g}".replace(".", "p"), "r": r, "eps": eps,
                    "sigma": sig, "sweeps": [sweep]})
    add(0.08, 1, 0, "unloaded")
    for s in SIGMAS:
        add(0.08, 78, s, "sigma")
    for e in EPSS:
        add(0.08, e, 0, "eps")
    add(0.08, 1, 0.6, "sigma_only")
    for r in RADII:
        if r != 0.08:
            add(r, 1, 0, "radius_unloaded")
            add(r, 78, 0.6, "radius")
    return out


def py(*args, check=True):
    cmd = [sys.executable, "-m", "birdcage_palace", *map(str, args)]
    r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    if check and r.returncode:
        raise RuntimeError(f"{' '.join(cmd)}\n{r.stdout}\n{r.stderr}")
    return r.stdout


def base_dir(r: float) -> Path:
    return RUNS / f"base_r{round(r * 1000):03d}"


def setup(a) -> None:
    RUNS.mkdir(parents=True, exist_ok=True)
    for r in RADII:
        b = base_dir(r)
        if not (b / "build.json").exists():
            print(py("build", "--symmetry", "eighth", "--sym-bc", "x=pec,y=pmc,z=pec", "--c-leg", C_TUNED, "--phantom",
                     "--phantom-radius", r, "--phantom-length", LENGTH, "--phantom-eps", 78, "--phantom-sigma", 0.6,
                     "--modes", 4, "--out", b).strip().splitlines()[0])
    for c in cases():
        b = base_dir(c["r"])
        mat = ["--phantom-eps", c["eps"], "--phantom-sigma", c["sigma"], "--phantom-loss-model", "conductivity",
               "--probe-h", PROBE_H]
        for role, bc, extra in (("em", "x=pec,y=pmc,z=pec", ["--copper-sigma", 0, "--esr", 0, "--modes", 4]),
                                ("me", "x=pmc,y=pec,z=pec", ["--copper-sigma", 0, "--esr", 0, "--modes", 4]),
                                ("full", "x=pec,y=pmc,z=pec", ["--copper-sigma", CU, "--esr", ESR, "--modes", 1])):
            d = RUNS / c["name"] / role
            if (d / "config.json").exists():
                continue
            py("variant", "--case", b, "--out", d, "--sym-bc", bc, *mat, *extra)
            if role == "full":
                # one mode near f0: the nonlinear (copper) solve stalls on the strongly damped phantom modes
                cfg = json.loads((d / "config.json").read_text())
                cfg["Solver"]["Eigenmode"].update(N=1, Target=0.9 * CoilDesign().f_larmor / 1e9)
                (d / "config.json").write_text(json.dumps(cfg, indent=1))
    (RUNS / "cases.json").write_text(json.dumps(cases(), indent=1))
    print(f"{len(cases())} cases, {3 * len(cases())} solves in {RUNS}")


def run(a) -> None:
    log = RUNS / "timings.jsonl"
    for c in cases():
        if a.only and not c["name"].startswith(a.only):
            continue
        for role in ("em", "me", "full"):
            d = RUNS / c["name"] / role
            if (d / "postpro" / "eig.csv").exists() and json.loads((d / "solve_stats.json").read_text()).get("exit") == 0:
                continue
            t0 = time.time()
            out = py("solve", "--case", d, "--np", a.np, check=False)
            st = json.loads((d / "solve_stats.json").read_text()) if (d / "solve_stats.json").exists() else {}
            rec = {"case": c["name"], "role": role, "wall_s": round(time.time() - t0, 1), **st,
                   "t_end": time.strftime("%Y-%m-%d %H:%M:%S")}
            with open(log, "a") as fh:
                fh.write(json.dumps(rec) + "\n")
            print(json.dumps(rec), flush=True)


# ---------------------------------------------------------------- analysis
def _unfold_E(case: Path, E: np.ndarray) -> np.ndarray:
    """Probe E from the representative points to the full point set (PEC plane: E(Mr) = -M E(r); PMC: +M E(r))."""
    pr = json.loads((case / "probes.json").read_text())
    bc = pr["sym_bc"]
    sig = {i: (-1.0 if bc[a] == "pec" else 1.0) for i, a in enumerate("xyz") if a in bc}
    out = E[:, np.array(pr["rep_index"]), :].copy()
    for j, fl in enumerate(pr["flips"]):
        for i in fl:
            f = np.full(3, sig[i])
            f[i] *= -1
            out[:, j, :] *= f
    return out


def _probe_E(case: Path) -> np.ndarray:
    h, d = M._read_csv(case / "postpro" / "probe-E.csv")
    cols = M._complex_columns(h, d, r"E_([xyz])\[(\d+)\]")
    n = max(int(k[1]) for k in cols)
    E = np.zeros((len(d), n, 3), complex)
    for (ax, i), v in cols.items():
        E[:, int(i) - 1, "xyz".index(ax)] = v
    return _unfold_E(case, E)


def _energies(case: Path) -> np.ndarray:
    h, d = M._read_csv(case / "postpro" / "domain-E.csv")
    return d


def _k1_index(case: Path) -> int:
    res = M.analyse(case)
    hp = res.get("homogeneous_pair")
    if not hp:
        raise RuntimeError(f"{case}: no k=1 mode")
    return hp["m"][0] - 1, res


def analyse_case(c: dict, dirs: tuple | None = None) -> dict:
    root = RUNS / c["name"]
    em, me, full = (Path(d) for d in dirs) if dirs else (root / "em", root / "me", root / "full")
    ia, ra = _k1_index(em)
    ib, rb = _k1_index(me)
    La, Lb = M.load(em), M.load(me)
    pr = json.loads((em / "probes.json").read_text())
    pts = np.array(pr["points"])
    n_ax, n_cor = pr["axial"], pr["coronal"]
    sl = {"axial": slice(0, n_ax), "coronal": slice(n_ax, n_ax + n_cor), "roi": slice(n_ax + n_cor, None)}
    Ba, Bb = La["B"][ia], Lb["B"][ib]
    Ea, Eb = _probe_E(em)[ia], _probe_E(me)[ib]
    # quadrature: same alignment / normalisation as modes.quadrature, applied to E as well
    roi = sl["roi"]
    def align(B):
        v = B[roi].sum(0)
        j = np.argmax(np.abs(v))
        return np.exp(-1j * np.angle(v[j]))
    pa, pb = align(Ba), align(Bb)
    na, nb = np.linalg.norm((Ba * pa)[roi]), np.linalg.norm((Bb * pb)[roi])
    ca, cb = pa / na, pb / nb
    best = None
    for s in (1, -1):
        B = ca * Ba + s * 1j * cb * Bb
        bp = (B[:, 0] + 1j * B[:, 1]) / 2
        if best is None or np.abs(bp).sum() > np.abs(best[1]).sum():
            bm = np.conj(B[:, 0] - 1j * B[:, 1]) / 2
            best = (s, bp, bm, ca * Ea + s * 1j * cb * Eb, (abs(ca), abs(cb)))
    s, bp, bm, Eq, (aa, ab) = best
    scale = 1e-6 / np.abs(bp[roi]).mean()   # mean |B1+| over the ROI = 1 uT
    bp, bm, Eq = bp * scale, bm * scale, Eq * scale
    r_xy = np.hypot(pts[:, 0], pts[:, 1])
    # 3 mm inside the surface: the mesh facets the cylinder, so a probe just inside can fall in an air element
    inph = (r_xy < c["r"] - 0.003) & (np.abs(pts[:, 2]) < LENGTH / 2 - 0.003)
    ax = np.zeros(len(pts), bool); ax[sl["axial"]] = True
    cor = np.zeros(len(pts), bool); cor[sl["coronal"]] = True
    Emag = np.sqrt((np.abs(Eq) ** 2).sum(1))
    sar = 0.5 * c["sigma"] * Emag ** 2          # W/m^3 at 1 uT mean B1+ (sigma |E|^2 / 2, no density)
    ab1 = np.abs(bp)
    sl_ph = ax & inph
    # stored energy and power, for B1+ per sqrt(W): the eigenvectors' energies scale with |c|^2
    ea, eb = _energies(em)[ia], _energies(me)[ib]
    rec = {k: c[k] for k in ("name", "r", "eps", "sigma", "sweeps")}
    fa, fb = ra["modes"][ia]["f_mhz"], rb["modes"][ib]["f_mhz"]
    qa, qb = ra["modes"][ia]["Q"], rb["modes"][ib]["Q"]
    rec.update(f_mhz=[fa, fb], f_mean_mhz=(fa + fb) / 2, Q_phantom=[qa, qb],
               b1plus_roi_cv=float(ab1[roi].std() / ab1[roi].mean()),
               b1minus_over_b1plus_roi=float(np.abs(bm[roi]).mean() / ab1[roi].mean()),
               b1plus_phantom_slice_cv=float(ab1[sl_ph].std() / ab1[sl_ph].mean()),
               b1plus_centre_over_mean=float(ab1[ax][np.argmin(r_xy[ax])] / ab1[roi].mean()),
               E_peak_phantom_vpm=float(Emag[inph].max()), E_mean_phantom_axial_vpm=float(Emag[sl_ph].mean()),
               sar_like_peak_wpm3=float(sar[inph].max()), sar_like_mean_axial_wpm3=float(sar[sl_ph].mean()),
               k_purity=[ra["modes"][ia].get("k_purity"), rb["modes"][ib].get("k_purity")],
               energy_cols="E_elec, E_mag, E_cap, E_ind (J), one eighth of the coil, Palace eigenvector norm",
               energy_em=ea[1:].tolist(), energy_me=eb[1:].tolist(), quad_coeff=[float(aa * scale), float(ab * scale)],
               currents_em=np.abs(La["I"][ia]).tolist() if "I" in La else None)
    if (full / "solve_stats.json").exists() and (full / "postpro" / "eig.csv").exists():
        h, d = M._read_csv(full / "postpro" / "eig.csv")
        if d.ndim == 2 and len(d):
            rec["f_full_mhz"] = float(d[0, 1] * 1e3)
            rec["Q_L"] = float(d[0, 3])
    for role, d in (("em", em), ("me", me), ("full", full)):
        p = d / "solve_stats.json"
        if p.exists():
            st = json.loads(p.read_text())
            rec[f"solve_{role}"] = {k: st.get(k) for k in ("wall_s", "unknowns", "ranks", "exit")}
    # maps for the figures (axial z = 0 plane and coronal y = 0 plane), at 1 uT mean B1+
    rec["_maps"] = {"axial_xy": pts[sl["axial"], :2].tolist(), "coronal_xz": pts[sl["coronal"]][:, [0, 2]].tolist(),
                    "b1p_axial": ab1[sl["axial"]].tolist(), "b1p_coronal": ab1[sl["coronal"]].tolist(),
                    "E_axial": Emag[sl["axial"]].tolist(), "E_coronal": Emag[sl["coronal"]].tolist(),
                    # complex quadrature transverse B (for the time animation) on the axial plane
                    "Bq_axial_re": np.real((ca * Ba + s * 1j * cb * Bb)[sl["axial"], :2] * scale).tolist(),
                    "Bq_axial_im": np.imag((ca * Ba + s * 1j * cb * Bb)[sl["axial"], :2] * scale).tolist()}
    return rec


def analyse(a) -> None:
    recs, maps = [], {}
    for c in cases():
        root = RUNS / c["name"]
        if not all((root / r / "solve_stats.json").exists() for r in ("em", "me")):
            continue
        try:
            rec = analyse_case(c)
        except Exception as e:  # noqa: BLE001
            print(f"{c['name']}: {e}")
            continue
        maps[c["name"]] = rec.pop("_maps")
        recs.append(rec)
    unl = {r["r"]: r for r in recs if r["eps"] == 1 and r["sigma"] == 0}
    for r in recs:
        u = unl.get(r["r"])
        if u:
            r["f_shift_mhz"] = r["f_mean_mhz"] - u["f_mean_mhz"]
            r["f_shift_rel"] = r["f_mean_mhz"] / u["f_mean_mhz"] - 1
            if "Q_L" in u and "Q_L" in r:
                r["Q_U"] = u["Q_L"]
                r["QU_over_QL"] = u["Q_L"] / r["Q_L"]
    (RUNS / "maps.json").write_text(json.dumps(maps))
    out = {"what": "Loading study: tuned 16-rung low-pass birdcage (D 28 cm, L 25 cm, 34 cm closed PEC/copper shield, "
                   "C 1.828 pF per rung, not retuned) with a centred homogeneous cylindrical phantom (length 20 cm). "
                   "Palace eigenmode, eighth model, production mesh density (5 / 2.5 / 35 mm, 20 mm in the phantom), "
                   "order 2. Simulation only.",
           "method": {"Q_phantom": "Palace Q of the k=1 member with lossless coil (phantom loss only)",
                      "Q_L": "Palace Q with phantom + copper (5.8e7 S/m on strips and shield, Conductivity boundary) "
                             "+ capacitor ESR 0.68 ohm (parallel R at f0)",
                      "Q_U": "Q_L of the eps_r = 1, sigma = 0 case on the same mesh",
                      "B1": "quadrature B1+ of the two k=1 members (two eighth runs), normalised to mean |B1+| = 1 uT "
                            "over the ROI r < 8.4 cm, |z| < 4.17 cm; phantom slice = z = 0 plane inside the phantom",
                      "E": "|E| of the same quadrature field (1 uT mean B1+); SAR-like = sigma |E|^2 / 2 in W/m^3 "
                           "(no tissue density), at the probes (z = 0 and y = 0 planes, 3D ROI grid) more than 3 mm inside the phantom"},
           "cases": recs}
    (ROOT / "results" / "loading_sweep.json").write_text(json.dumps(out, indent=1))
    print(f"{len(recs)} cases -> results/loading_sweep.json")
    for r in recs:
        print(f"{r['name']:<16} f {r['f_mean_mhz']:8.3f}  dF {r.get('f_shift_rel', 0) * 100:+6.2f} %  Qph {r['Q_phantom'][0]:9.1f}"
              f"  QL {r.get('Q_L', float('nan')):7.1f}  QU/QL {r.get('QU_over_QL', float('nan')):6.2f}  B1+cv {r['b1plus_roi_cv'] * 100:5.1f} %"
              f"  slice {r['b1plus_phantom_slice_cv'] * 100:5.1f} %  Epk {r['E_peak_phantom_vpm']:6.1f}  SARpk {r['sar_like_peak_wpm3']:7.2f}")


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("setup").set_defaults(fn=setup)
    r = sub.add_parser("run")
    r.add_argument("--np", type=int, default=4)
    r.add_argument("--only", default=None)
    r.set_defaults(fn=run)
    sub.add_parser("analyse").set_defaults(fn=analyse)
    a = ap.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
