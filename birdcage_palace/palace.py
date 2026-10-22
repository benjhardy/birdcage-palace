"""Palace eigenmode configuration and runner."""
from __future__ import annotations

import json
import os
import platform
import re
import shlex
import subprocess
import time
from pathlib import Path

import numpy as np

from .design import CoilDesign


def probe_points(des: CoilDesign, h_plane: float = 0.006, h_roi: float = 0.012) -> dict:
    """Probe sets: axial plane z=0, coronal plane y=0, and the 3D ROI (r < 0.6 R, |z| < length/6)."""
    R, l = des.radius, des.length
    rmax = 0.92 * R
    g = np.arange(-rmax, rmax + 1e-9, h_plane)
    ax = [(x, y, 0.0) for x in g for y in g if x * x + y * y <= rmax ** 2]
    zc = np.arange(-0.75 * l, 0.75 * l + 1e-9, h_plane)
    if des.shield_len:
        zc = zc[np.abs(zc) < des.shield_len / 2 - 0.004]
    cor = [(x, 0.0, z) for z in zc for x in g]
    gr = np.arange(-0.6 * R, 0.6 * R + 1e-9, h_roi)
    gz = np.arange(-l / 6, l / 6 + 1e-9, h_roi)
    roi = [(x, y, z) for z in gz for x in gr for y in gr if x * x + y * y <= (0.6 * R) ** 2]
    # nudge points off exact symmetry planes so no probe lands on an interior mesh face of the cap patches
    def _f(p):
        return [round(float(c) + 1.1e-5, 7) for c in p]
    return {"axial": [_f(p) for p in ax], "coronal": [_f(p) for p in cor], "roi": [_f(p) for p in roi]}


def symmetry_reps(points: list, planes: list[str]) -> tuple[list, list[int], list[list[int]]]:
    """Map every full-domain probe to its mirror image in the modelled part (|x|, |y|, |z| on the cut axes).

    Returns the unique representative points, the index of each full point's representative, and which axes were
    flipped to reach it (the field there is rebuilt from the representative by the mirror rule of each plane).
    """
    reps, key_idx, rep_index, flips = [], {}, [], []
    for p in points:
        q, fl = list(p), []
        for i, a in enumerate("xyz"):
            if a in planes and q[i] < 0:
                q[i], fl = -q[i], fl + [i]
        key = tuple(round(c, 7) for c in q)
        if key not in key_idx:
            key_idx[key] = len(reps)
            reps.append(q)
        rep_index.append(key_idx[key])
        flips.append(fl)
    return reps, rep_index, flips


EPS0 = 8.8541878128e-12


def materials(des: CoilDesign, attrs: dict) -> list:
    """Palace materials: air, plus the phantom with its loss as a conductivity or as the equivalent loss tangent."""
    mats = [{"Attributes": [attrs["air"]], "Permittivity": 1.0, "Permeability": 1.0}]
    if attrs.get("phantom"):
        p = des.phantom
        m = {"Attributes": [attrs["phantom"]], "Permittivity": p.eps_r, "Permeability": 1.0}
        if p.sigma > 0:
            if p.loss_model == "losstan":
                f = p.f_ref or des.f_larmor
                m["LossTan"] = p.sigma / (2 * np.pi * f * EPS0 * p.eps_r)
            else:
                m["Conductivity"] = p.sigma
        mats.append(m)
    return mats


def metal_boundaries(des: CoilDesign, attrs: dict) -> tuple[list, list]:
    """PEC attributes and (with copper_sigma) Palace Conductivity boundaries for the strips and the shield.

    The strips are internal sheets: Palace applies one surface impedance (1+i)/(sigma delta) to the sheet, i.e. the
    loss of the whole strip current on one face. A real strip shares its current between both faces, so this is an
    upper bound on the strip loss (at most 2x). The shield is an outer boundary (External: true).
    """
    pec, cond = [], []
    if des.copper_sigma:
        cond.append({"Attributes": [attrs["coil"]], "Conductivity": des.copper_sigma})
        if attrs.get("shield"):
            cond.append({"Attributes": [attrs["shield"]], "Conductivity": des.copper_sigma, "External": True})
    else:
        pec = [attrs["coil"]] + ([attrs["shield"]] if attrs.get("shield") else [])
    return pec, cond


def write_config(des: CoilDesign, case: Path, n_modes: int = 20, target_hz: float | None = None,
                 tol: float = 1e-8, save: int = 0, probes: bool = True, h_plane: float = 0.006) -> dict:
    case = Path(case)
    info = json.loads((case / "build.json").read_text())
    a = info["attrs"]
    mats = materials(des, a)
    pec, cond = metal_boundaries(des, a)
    # Symmetry planes: PEC is set explicitly; PMC (tangential H = 0) is Palace's natural boundary, so a PMC plane is
    # simply left out of the boundary list.
    sym = a.get("sym") or {}
    pec += [sym[ax] for ax in sorted(sym) if des.sym_bc.get(ax) == "pec"]
    bnd = {"PEC": {"Attributes": pec}}
    if cond:
        bnd["Conductivity"] = cond
    if a.get("absorbing"):
        bnd["Absorbing"] = {"Attributes": [a["absorbing"]], "Order": 2}  # 2nd order: see geometry.air_clearance
    lp = []
    for i, p in enumerate(info["ports"], start=1):
        c = des.c_leg if p["kind"] == "leg" else des.c_ring
        if p["kind"] == "leg" and des.c_leg_rungs:
            c = des.c_leg_rungs[p["n"]]
        el = {"Index": i, "Attributes": [p["attr"]], "Direction": p["dir"], "C": c * p.get("c_scale", 1.0)}
        if des.esr:
            # series ESR -> parallel R at w0 (exact at w0); a cut element keeps the per-square value: R / c_scale
            w0 = 2 * np.pi * des.f_larmor
            qc = 1 / (w0 * c * des.esr)
            el["R"] = des.esr * (1 + qc * qc) / p.get("c_scale", 1.0)
        lp.append(el)
    bnd["LumpedPort"] = lp
    if target_hz is None:  # search from well below the lowest mode we care about
        target_hz = 0.55 * des.f_larmor if des.coil_type == "lowpass" else 0.25 * des.f_larmor
        if a.get("absorbing"):
            # with an absorbing boundary a low target sits among strongly damped box modes and Krylov-Schur stalls
            # (unshielded LP16 at 0.55 f0: no mode converged in 10 restarts; at 0.9 f0: 254 s). 0.75 f0 still catches
            # a k = 1 that the first guess puts 25 % low.
            target_hz = 0.75 * des.f_larmor
    domains = {"Materials": mats}
    pts = None
    if probes:
        pts = probe_points(des, h_plane=h_plane)
        allp = pts["axial"] + pts["coronal"] + pts["roi"]
        meta = {k: len(v) for k, v in pts.items()} | {"points": allp}
        solve_pts = allp
        if des.planes:
            solve_pts, rep_index, flips = symmetry_reps(allp, des.planes)
            meta |= {"rep_points": solve_pts, "rep_index": rep_index, "flips": flips, "planes": des.planes,
                     "sym_bc": {ax: des.sym_bc[ax] for ax in des.planes}}
        domains["Postprocessing"] = {"Probe": [{"Index": i + 1, "Center": p} for i, p in enumerate(solve_pts)]}
        (case / "probes.json").write_text(json.dumps(meta))
    if a.get("absorbing"):
        # An absorbing boundary adds strongly damped "box" modes (Q < 1) that the Krylov-Schur solver converges
        # slowly: in the unshielded package test a 6-mode solve found the k = 1 pair in 2 restarts and then spent
        # 20+ min (and once 48 min) on the next modes. Ask for the k = 1 pair only, at 1e-6 (plenty for f and Q).
        tol = max(tol, 1e-6)
        n_modes = min(n_modes, 2)
    cfg = {
        "Problem": {"Type": "Eigenmode", "Verbose": 1, "Output": "postpro"},
        "Model": {"Mesh": info["mesh"], "L0": 1.0},
        "Domains": domains,
        "Boundaries": bnd,
        "Solver": {"Order": des.order, "Device": "CPU",
                   "Eigenmode": {"Target": target_hz / 1e9, "N": n_modes, "Tol": tol, "Save": save},
                   "Linear": {"Type": "MUMPS", "MGMaxLevels": 1, "KSPType": "GMRES", "Tol": 1e-10, "MaxIts": 200,
                              "EstimatorTol": 1e-3, "EstimatorMaxIts": 200, "EstimatorMG": True}},
    }
    if save == 0:
        cfg["Problem"]["OutputFormats"] = {"Paraview": False, "GridFunction": False}
    (case / "config.json").write_text(json.dumps(cfg, indent=1))
    return cfg


def _wslpath(p: Path) -> str:
    p = Path(p).resolve()
    drive = p.drive.rstrip(":").lower()
    return f"/mnt/{drive}" + p.as_posix()[2:]


def palace_busy(distro: str = "Ubuntu-24.04") -> list[str]:
    cmd = ["pgrep", "-af", "palace-x86_64.bin"]
    if platform.system() == "Windows":
        cmd = ["wsl", "-d", distro, "--"] + cmd
    r = subprocess.run(cmd, capture_output=True, text=True)
    return [ln for ln in r.stdout.splitlines() if "palace-x86_64.bin" in ln and "pgrep" not in ln]


def run(case: Path, np_: int = 4, palace_cmd: str = "run_palace {cfg} {np}", distro: str = "Ubuntu-24.04") -> dict:
    """Run Palace on case/config.json. On Windows the command runs inside WSL (Palace via Apptainer)."""
    case = Path(case)
    pcmd = palace_cmd.format(cfg="config.json", np=np_)
    # GNU time records the peak memory per rank when it is installed; without it Palace still runs
    inner = (f"if [ -x /usr/bin/time ]; then /usr/bin/time -v {pcmd} > palace.log 2> time.log; "
             f"else {pcmd} > palace.log 2> time.log; fi")
    t0 = time.time()
    if platform.system() == "Windows":
        cmd = ["wsl", "-d", distro, "--", "bash", "-lc", f"cd {shlex.quote(_wslpath(case))} && {inner}"]
    else:
        cmd = ["bash", "-lc", f"cd {shlex.quote(str(case.resolve()))} && {inner}"]
    r = subprocess.run(cmd)
    wall = time.time() - t0
    log = (case / "palace.log").read_text(errors="replace") if (case / "palace.log").exists() else ""
    tlog = (case / "time.log").read_text(errors="replace") if (case / "time.log").exists() else ""
    m = re.findall(r"global unknowns = (\d+)", log) or re.findall(r"(\d+) global unknowns", log)
    rss = re.findall(r"Maximum resident set size \(kbytes\): (\d+)", tlog)
    peak = re.findall(r"[Pp]eak.*memory.*?([\d.]+ ?[GMK]i?B)", log)
    stats = {"ranks": np_, "wall_s": round(wall, 1), "exit": r.returncode,
             "unknowns": int(m[-1]) if m else None, "peak_mem_reported": peak[-1] if peak else None,
             "max_rss_one_rank_gb": round(int(rss[-1]) / 1048576, 2) if rss else None}
    (case / "solve_stats.json").write_text(json.dumps(stats, indent=1))
    return stats


def write_driven_config(des: CoilDesign, case: Path, f_min_hz: float = 100e6, f_max_hz: float = 160e6,
                        f_step_hz: float = 2e6, r_ref: float = 500.0, adaptive_tol: float = 1e-5,
                        max_samples: int = 40) -> dict:
    """Port-tuning model: every capacitor gap becomes a resistive lumped port (reference r_ref), each its own excitation.

    The driven sweep returns the multi-port S-matrix of the bare coil (no lumped C), which a circuit simulator then
    terminates with the tuning capacitors and the drive network. A gap cut by a symmetry plane keeps the same
    per-square resistance (Rs = R w / l): R is divided by the same factor the capacitor is multiplied by.
    """
    case = Path(case)
    info = json.loads((case / "build.json").read_text())
    a = info["attrs"]
    mats = materials(des, a)
    pec, cond = metal_boundaries(des, a)
    sym = a.get("sym") or {}
    pec += [sym[ax] for ax in sorted(sym) if des.sym_bc.get(ax) == "pec"]
    bnd = {"PEC": {"Attributes": pec}}
    if cond:
        bnd["Conductivity"] = cond
    if a.get("absorbing"):
        bnd["Absorbing"] = {"Attributes": [a["absorbing"]], "Order": 2}  # 2nd order: see geometry.air_clearance
    bnd["LumpedPort"] = [{"Index": i, "Attributes": [p["attr"]], "Direction": p["dir"],
                          "R": r_ref / p.get("c_scale", 1.0), "Excitation": i}
                         for i, p in enumerate(info["ports"], start=1)]
    cfg = {
        "Problem": {"Type": "Driven", "Verbose": 1, "Output": "postpro",
                    "OutputFormats": {"Paraview": False, "GridFunction": False}},
        "Model": {"Mesh": info["mesh"], "L0": 1.0},
        "Domains": {"Materials": mats},
        "Boundaries": bnd,
        "Solver": {"Order": des.order, "Device": "CPU",
                   "Driven": {"Samples": [{"Type": "Linear", "MinFreq": f_min_hz / 1e9, "MaxFreq": f_max_hz / 1e9,
                                           "FreqStep": f_step_hz / 1e9, "SaveStep": 0}],
                              "AdaptiveTol": adaptive_tol, "AdaptiveMaxSamples": max_samples},
                   "Linear": {"Type": "MUMPS", "MGMaxLevels": 1, "KSPType": "GMRES", "Tol": 1e-10, "MaxIts": 200,
                              "EstimatorTol": 1e-3, "EstimatorMaxIts": 200}},
    }
    (case / "config.json").write_text(json.dumps(cfg, indent=1))
    (case / "driven.json").write_text(json.dumps({"r_ref": r_ref, "f_min_hz": f_min_hz, "f_max_hz": f_max_hz,
                                                  "f_step_hz": f_step_hz, "adaptive_tol": adaptive_tol}, indent=1))
    return cfg
