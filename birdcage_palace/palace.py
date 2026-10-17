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


def write_config(des: CoilDesign, case: Path, n_modes: int = 20, target_hz: float | None = None,
                 tol: float = 1e-8, save: int = 0, probes: bool = True) -> dict:
    case = Path(case)
    info = json.loads((case / "build.json").read_text())
    a = info["attrs"]
    mats = [{"Attributes": [a["air"]], "Permittivity": 1.0, "Permeability": 1.0}]
    if a.get("phantom"):
        p = des.phantom
        mats.append({"Attributes": [a["phantom"]], "Permittivity": p.eps_r, "Permeability": 1.0, "Conductivity": p.sigma})
    pec = [a["coil"]] + ([a["shield"]] if a.get("shield") else [])
    bnd = {"PEC": {"Attributes": pec}}
    if a.get("absorbing"):
        bnd["Absorbing"] = {"Attributes": [a["absorbing"]], "Order": 1}
    lp = []
    for i, p in enumerate(info["ports"], start=1):
        c = des.c_leg if p["kind"] == "leg" else des.c_ring
        lp.append({"Index": i, "Attributes": [p["attr"]], "Direction": p["dir"], "C": c})
    bnd["LumpedPort"] = lp
    if target_hz is None:  # search from well below the lowest mode we care about
        target_hz = 0.55 * des.f_larmor if des.coil_type == "lowpass" else 0.25 * des.f_larmor
    domains = {"Materials": mats}
    pts = None
    if probes:
        pts = probe_points(des)
        allp = pts["axial"] + pts["coronal"] + pts["roi"]
        domains["Postprocessing"] = {"Probe": [{"Index": i + 1, "Center": p} for i, p in enumerate(allp)]}
        (case / "probes.json").write_text(json.dumps({k: len(v) for k, v in pts.items()} | {"points": allp}))
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
    inner = f"/usr/bin/time -v {palace_cmd.format(cfg='config.json', np=np_)} > palace.log 2> time.log"
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
