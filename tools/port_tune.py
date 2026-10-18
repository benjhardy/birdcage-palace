"""Port tuning of the birdcage: the 16-port S-matrix of the bare coil (every capacitor gap a port) from Palace, then
the capacitors and the quadrature drive network tuned in circuitRF (open-source circuit simulator, headless CLI in WSL).

    python tools/port_tune.py assemble --cases runs/pt_e --out runs/pt_e/birdcage16.s16p
    python tools/port_tune.py resonance --snp runs/pt_e/birdcage16.s16p [--c 1.828e-12]
    python tools/port_tune.py tune --snp runs/pt_e/birdcage16.s16p --name unloaded --esr 0.68
    python tools/port_tune.py check --cases runs/ptc_e --full runs/ptc_full

The four cases <prefix>_em/_me/_mm/_ee are eighth (or quarter) port cases with planes between rungs; the full matrix
is rebuilt from them exactly (birdcage_palace/sparams.py).
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from birdcage_palace import sparams  # noqa: E402

GAMMA_BAR = 42.577478e6
CRF = "$HOME/.local/share/circuitRF/current/circuitRF"


def wslpath(p: Path) -> str:
    p = Path(p).resolve()
    return f"/mnt/{p.drive.rstrip(':').lower()}" + p.as_posix()[2:]


def wsl(cmd: str, cwd: Path, timeout=3600) -> subprocess.CompletedProcess:
    return subprocess.run(["wsl.exe", "-d", "Ubuntu-24.04", "--", "bash", "-lc", f"cd {wslpath(cwd)} && {cmd}"],
                          capture_output=True, text=True, timeout=timeout)


def read_snp(path: Path):
    import skrf as rf
    n = rf.Network(str(path))
    return n.f, n.s, float(np.real(n.z0[0, 0]))


# ------------------------------------------------------------------------------------------- network algebra (python)
def y_matrix(S, r):
    I = np.eye(S.shape[-1])
    return np.linalg.solve(r * (I + S), (I - S))


def resonances(f, S, r, c, f_lo=110e6, f_hi=150e6, step=2e3):
    """Resonances of the lossless coil with every gap closed by c (scalar or per-gap array): frequencies where an
    eigenvalue of the gap susceptance matrix Im(Y + jwC) crosses zero. S is cubic-spline interpolated (the bare-coil
    S is smooth), crossings are located by linear interpolation on a fine grid."""
    from scipy.interpolate import CubicSpline
    fi = np.arange(max(f_lo, f[0]), min(f_hi, f[-1]), step)
    Si = CubicSpline(f, S.real, axis=0)(fi) + 1j * CubicSpline(f, S.imag, axis=0)(fi)
    Y = y_matrix(Si, r)
    cc = np.broadcast_to(np.atleast_1d(c), (S.shape[1],))
    Yt = Y + 1j * 2 * np.pi * fi[:, None, None] * np.einsum("i,ij->ij", cc, np.eye(S.shape[1]))[None]
    B = (Yt - np.conj(np.transpose(Yt, (0, 2, 1)))) / 2j
    ev = np.linalg.eigvalsh(B)
    out = []
    for j in range(ev.shape[1]):
        e = ev[:, j]
        idx = np.where(np.sign(e[:-1]) != np.sign(e[1:]))[0]
        for i in idx:
            if abs(e[i]) + abs(e[i + 1]) < 1.0:  # a zero crossing, not a pole (poles jump through +-inf)
                out.append(float(fi[i] - e[i] * (fi[i + 1] - fi[i]) / (e[i + 1] - e[i])))
    return sorted(out)


def k1_frequency(f, S, r, c, f_lo=110e6, f_hi=150e6):
    """Lowest resonance (the k = 1 pair: mean of the two lowest crossings when both are found)."""
    rs = resonances(f, S, r, c, f_lo, f_hi)
    if not rs:
        raise ValueError("no resonance in the band")
    pair = rs[:2] if len(rs) > 1 and rs[1] / rs[0] - 1 < 0.01 else rs[:1]
    return float(np.mean(pair)), rs


def c_for_frequency(f, S, r, f_target, c0=1.8e-12):
    """Common gap capacitor that puts the k=1 resonance on f_target (secant on 1/f^2, linear in C)."""
    c = [c0, c0 * 1.02]
    fr = [k1_frequency(f, S, r, cc)[0] for cc in c]
    for _ in range(8):
        y = [1 / x ** 2 for x in fr]
        slope = (y[-1] - y[-2]) / (c[-1] - c[-2])
        cn = c[-1] + (1 / f_target ** 2 - y[-1]) / slope
        c.append(cn)
        fr.append(k1_frequency(f, S, r, cn)[0])
        if abs(fr[-1] / f_target - 1) < 1e-7:
            break
    return c[-1], fr[-1]


# ------------------------------------------------------------------------------------------- circuitRF netlist
ISO = (2, 10)  # rungs on the bisector of the two drive rungs (45 deg from each): a trim here couples the two k=1
#                members, which cancels the residual coupling of the drive gaps through the other modes
RANGES = {"Ct": (1.0, 3.0, "lin"), "Ci": (0.5, 3.5, "lin"), "Cd1": (0.5, 5.0, "lin"), "Cd2": (0.5, 5.0, "lin"),
          "Cm1": (0.05, 300.0, "log"), "Cm2": (0.05, 300.0, "log"), "Cp1": (2.0, 200.0, "log"), "Cp2": (2.0, 200.0, "log")}
FEED = "divider"  # "series": series match capacitor straight across the gap; "divider": the gap capacitor is split
#                   into Cd (top) and a tap capacitor Cp to the feed, then a series Cm to the 50 ohm port


def rung_var(i: int, drive=(0, 4), iso: bool = True) -> str:
    if i == drive[0]:
        return "Cd1"
    if i == drive[1]:
        return "Cd2"
    return "Ci" if (iso and i in ISO) else "Ct"


def write_cnl(work: Path, snp: Path, name: str, f0: float, esr: float, vals: dict, tuned: list[str] | None = None,
              drive=(0, 4), band_khz: float = 100.0, n_ports=16, iso: bool = True, span_mhz: float = 1.0,
              step_khz: float = 2.0, optimize: str = "algorithm=minimax maxiter=300 seed=1",
              tie: dict | None = None, feed: str | None = None) -> Path:
    """Netlist: the EM 16-port, each gap closed by its capacitor (with ESR), the two drive gaps also fed through a
    series match capacitor from a 50 ohm port. vals: Ct (most rungs), Ci (isolation rungs 2, 10), Cd1/Cd2 (drive gaps),
    Cm1/Cm2 (series match), in F. tie = {"Cd2": "Cd1", "Cm2": "Cm1"} makes the two drive ports share one value.
    Goals at f0 +- band: |S11|, |S22| <= -20 dB and |S21| <= -20 dB."""
    local = work / Path(snp).name
    if local.resolve() != Path(snp).resolve():
        shutil.copy(snp, local)
    tie = tie or {}
    feed = feed or FEED
    nodes = " ".join(f"g{i}" for i in range(n_ports))
    L = [f"; birdcage port tuning ({name}): Palace {n_ports}-port of the bare coil (every capacitor gap a port);",
         "; the capacitors (with ESR) and the quadrature drive network are added here",
         f"Resr = {esr:g} Ohm"]
    names = ("Ct", "Ci", "Cd1", "Cd2", "Cm1", "Cm2") + (("Cp1", "Cp2") if feed == "divider" else ())
    for k in names:
        L.append(f"{k} = {tie[k]}" if k in tie else f"{k} = {vals[k] * 1e12:.6f} pF")
    L.append(f'SnP:COIL {nodes} NumPorts={n_ports} File="{local.name}"')
    for i in range(n_ports):
        cv = rung_var(i, drive, iso)
        if feed == "divider" and i in drive:
            k = drive.index(i) + 1  # gap -> Cd -> tap a_k -> Cp -> ground; ESR in series with each capacitor
            L += [f"C:C{i} g{i} e{i} C={cv}", f"R:RS{i} e{i} a{k} R=Resr",
                  f"C:CP{k} a{k} t{k} C=Cp{k}", f"R:RP{k} t{k} 0 R=Resr"]
        elif esr > 0:
            L += [f"C:C{i} g{i} e{i} C={cv}", f"R:RS{i} e{i} 0 R=Resr"]
        else:
            L.append(f"C:C{i} g{i} 0 C={cv}")
    for k, i in enumerate(drive, start=1):
        node = f"a{k}" if feed == "divider" else f"g{i}"
        L += [f"C:CM{k} {node} p{k} C=Cm{k}", f"Port:P{k} p{k} 0 Num={k} Z=50 Ohm"]
    fa, fb = f0 - span_mhz * 1e6, f0 + span_mhz * 1e6
    L.append(f"analysis SP1 type=sparam start={fa / 1e6:.4f} stop={fb / 1e6:.4f} step={step_khz / 1e3:g} Unit=MHz")
    lo, hi = (f0 - band_khz * 1e3) / 1e6, (f0 + band_khz * 1e3) / 1e6
    if tuned:
        for k in tuned:
            mn, mx, sc = RANGES[k]
            L.append(f"tune {k} min={mn:g} pF max={mx:g} pF scale={sc} opt=1")
        L += [f"goal match1 = dB(SP1.S(1,1)) analysis=SP1 over=freq lo={lo:.4f} MHz hi={hi:.4f} MHz le -20",
              f"goal match2 = dB(SP1.S(2,2)) analysis=SP1 over=freq lo={lo:.4f} MHz hi={hi:.4f} MHz le -20",
              f"goal isol = dB(SP1.S(2,1)) analysis=SP1 over=freq lo={lo:.4f} MHz hi={hi:.4f} MHz le -20",
              f"optimize {optimize}"]
    p = work / f"tune_{name}.cnl"
    p.write_text("\n".join(L) + "\n")
    return p


def crf_opt(work: Path, cnl: Path) -> dict:
    r = wsl(f"{CRF} opt {cnl.name} --json", work)
    (work / f"{cnl.stem}.opt.json").write_text(r.stdout)
    if r.stderr:
        (work / f"{cnl.stem}.opt.err").write_text(r.stderr)
    return json.loads(r.stdout)


def crf_sparam(work: Path, cnl: Path, sets: dict, out_name: str) -> Path:
    s = " ".join(f"--set {k}={v * 1e12:.6f}pF" for k, v in sets.items())
    r = wsl(f"{CRF} sparam {cnl.name} {s} -o {out_name}", work)
    if r.returncode:
        raise RuntimeError(r.stdout + r.stderr)
    return work / out_name


def _val_pf(txt: str) -> float:
    m = re.match(r"\s*([-+\d.eE]+)\s*([pfnu]?)F?", txt)
    v = float(m.group(1))
    return v * {"p": 1e-12, "f": 1e-15, "n": 1e-9, "u": 1e-6, "": 1.0}[m.group(2)]


def rung_caps(vals: dict, n=16, drive=(0, 4), iso=True, feed: str | None = None) -> np.ndarray:
    """Capacitance across each gap with the feed disconnected (divider: Cd in series with Cp)."""
    feed = feed or FEED
    C = np.array([vals[rung_var(i, drive, iso)] for i in range(n)])
    if feed == "divider":
        for k, i in enumerate(drive, start=1):
            C[i] = 1 / (1 / vals[f"Cd{k}"] + 1 / vals[f"Cp{k}"])
    return C


def terminated(Y, w, vals, esr, drive=(0, 4), iso=True, feed: str | None = None):
    """Nodal admittance of the coil (bare-coil gap admittance Y) with the capacitors and the drive network.
    Nodes: 16 gaps, then for each drive port its tap (divider) and its port node. Returns (Yn, port node indices)."""
    feed = feed or FEED
    n = Y.shape[0]
    nd = len(drive)
    extra = (2 * nd) if feed == "divider" else nd
    Yn = np.zeros((n + extra, n + extra), complex)
    Yn[:n, :n] = Y

    def stamp(i, j, y):
        Yn[i, i] += y
        if j is not None:
            Yn[j, j] += y
            Yn[i, j] -= y
            Yn[j, i] -= y
    yc = lambda c: 1 / (esr + 1 / (1j * w * c))  # noqa: E731
    ports = []
    for i in range(n):
        if not (feed == "divider" and i in drive):
            stamp(i, None, yc(vals[rung_var(i, drive, iso)]))
    for k, i in enumerate(drive, start=1):
        if feed == "divider":
            a, pn = n + 2 * (k - 1), n + 2 * (k - 1) + 1
            stamp(i, a, yc(vals[f"Cd{k}"]))
            stamp(a, None, yc(vals[f"Cp{k}"]))
            stamp(a, pn, 1j * w * vals[f"Cm{k}"])
        else:
            pn = n + (k - 1)
            stamp(i, pn, 1j * w * vals[f"Cm{k}"])
        ports.append(pn)
    return Yn, ports


def python_two_port(f, S, r, vals: dict, esr, drive=(0, 4), iso=True, feed: str | None = None):
    """Cross-check: terminate the EM 16-port in python and return the 2-port S (50 ohm) at the drive ports."""
    from scipy.interpolate import CubicSpline
    fe, Se = S
    Si = CubicSpline(fe, Se.real, axis=0)(f) + 1j * CubicSpline(fe, Se.imag, axis=0)(f)
    Y = y_matrix(Si, r)
    out = np.zeros((len(f), 2, 2), complex)
    for k in range(len(f)):
        Yn, pn = terminated(Y[k], 2 * np.pi * f[k], vals, esr, drive, iso, feed)
        Z = np.linalg.inv(Yn)[np.ix_(pn, pn)]
        out[k] = np.linalg.solve(Z + 50 * np.eye(2), Z - 50 * np.eye(2))
    return out


def rung_currents(f0, S, r, vals, esr, e=(1, 1j), drive=(0, 4), iso=True, feed: str | None = None):
    """Gap (= rung) currents with the two ports driven by EMFs e behind 50 ohm (quadrature by default)."""
    from scipy.interpolate import CubicSpline
    fe, Se = S
    Si = CubicSpline(fe, Se.real, axis=0)([f0]) + 1j * CubicSpline(fe, Se.imag, axis=0)([f0])
    Y = y_matrix(Si, r)[0]
    Yn, pn = terminated(Y, 2 * np.pi * f0, vals, esr, drive, iso, feed)
    for p in pn:
        Yn[p, p] += 1 / 50
    J = np.zeros(Yn.shape[0], complex)
    for p, ek in zip(pn, e):
        J[p] = ek / 50
    V = np.linalg.solve(Yn, J)
    n = Y.shape[0]
    return Y @ V[:n]


def k1_content(I):
    P = np.abs(np.fft.fft(I)) ** 2
    return float(max(P[1], P[-1]) / P.sum())


def best_vals(opt: dict, start: dict) -> dict:
    v = dict(start)
    for x in opt["result"]["optimize"]["variables"]:
        v[x["key"]] = _val_pf(x["best"])
    return v


def read_crf_s2p(path: Path):
    import skrf as rf
    n = rf.Network(str(path))
    return n.f, n.s


def tune(snp: Path, work: Path, name: str, f0: float, esr: float, band_khz: float, c_eigen: float,
         iters: int = 0) -> dict:
    """Port tuning in three steps:
      1. resonance only: every gap the same C; the S-matrix gives the k=1 frequencies at the eigenmode C and the C that
         puts the pair on f0 (unloaded, lossless network; compared with the eigen solves);
      2. circuitRF on the rotation-averaged matrix (ideal 16-fold symmetry): Ct, Ci (isolation rungs), Cd and Cm shared
         by both drive ports;
      3. circuitRF on the matrix as solved (the mesh splits the k=1 pair slightly): the drive ports get their own Cd, Cm
         (per-port trims), starting from step 2.
    Returns everything for the figures and the README."""
    work.mkdir(parents=True, exist_ok=True)
    f, S, r = read_snp(snp)
    sym = snp.with_name(snp.stem + "_sym" + snp.suffix)
    fs, Ss, _ = read_snp(sym)
    out = {"snp": str(snp), "f0_hz": f0, "esr_ohm": esr, "band_khz": band_khz, "c_eigen": c_eigen}
    try:
        out["k1_at_c_eigen_hz"] = resonances(f, S, r, c_eigen, f0 * 0.9, f0 * 1.1)[:2]
        out["k1_at_c_eigen_sym_hz"] = resonances(fs, Ss, r, c_eigen, f0 * 0.9, f0 * 1.1)[:2]
        c_res, f_res = c_for_frequency(f, S, r, f0, c_eigen)
        out["c_resonance_only"] = c_res
        out["k1_at_c_resonance_only_hz"] = resonances(f, S, r, c_res, f0 * 0.95, f0 * 1.05)[:2]
    except ValueError as e:  # lossy (loaded) matrix: no clean susceptance zero crossing
        out["resonance_note"] = str(e)
        c_res = c_eigen
    if FEED == "divider":  # gap = Cd in series with a large tap Cp (series value ~ c_res); Cm nearly a short
        cp0 = 20e-12
        cd0 = 1 / (1 / c_res - 1 / cp0)
        v0 = {"Ct": c_res, "Ci": c_res, "Cd1": cd0, "Cd2": cd0, "Cp1": cp0, "Cp2": cp0, "Cm1": 30e-12, "Cm2": 30e-12}
        tie = {"Cd2": "Cd1", "Cp2": "Cp1", "Cm2": "Cm1"}
        free1 = ["Ct", "Ci", "Cd1", "Cp1", "Cm1"]
    else:
        v0 = {"Ct": c_res, "Ci": c_res, "Cd1": 0.85 * c_res, "Cd2": 0.85 * c_res, "Cm1": 0.6e-12, "Cm2": 0.6e-12}
        tie = {"Cd2": "Cd1", "Cm2": "Cm1"}
        free1 = ["Ct", "Ci", "Cd1", "Cm1"]
    out["feed"] = FEED
    out["start"] = v0
    cnl = write_cnl(work, sym, f"{name}_sym", f0, esr, v0, free1, band_khz=band_khz, tie=tie)
    o1 = crf_opt(work, cnl)
    v1 = best_vals(o1, v0)
    for k2, k1 in tie.items():
        v1[k2] = v1[k1]
    out["sym"] = {"vals": v1, "opt": {k: o1["result"]["optimize"][k] for k in ("outcome", "finishReason", "bestCost", "iterations", "evaluations")},
                  "goals": o1["result"]["optimize"]["goals"]}
    tuned = free1 + list(tie)
    cnl2 = write_cnl(work, snp, name, f0, esr, v1, tuned, band_khz=band_khz)
    o2 = crf_opt(work, cnl2)
    v2 = best_vals(o2, v1)
    out["raw"] = {"vals": v2, "opt": {k: o2["result"]["optimize"][k] for k in ("outcome", "finishReason", "bestCost", "iterations", "evaluations")},
                  "goals": o2["result"]["optimize"]["goals"]}
    # wide sweeps (circuitRF) before and after; before = every gap at the eigenmode C, same match capacitors
    vb = dict(v2, Ct=c_eigen, Ci=c_eigen)  # every rung at the eigenmode C, drive network as tuned
    if FEED == "divider":
        for k in (1, 2):
            vb[f"Cd{k}"] = 1 / (1 / c_eigen - 1 / v2[f"Cp{k}"])
    else:
        vb["Cd1"] = vb["Cd2"] = c_eigen
    for tag, vals in (("start", v0), ("sym", v1), ("raw", v2)):
        Iq = rung_currents(f0, (f, S), r, vals, esr)
        out.setdefault("k1_content_quadrature", {})[tag] = k1_content(Iq)
        out.setdefault("rung_current_max_over_min", {})[tag] = float(np.abs(Iq).max() / np.abs(Iq).min())
    for tag, vals in (("before", vb), ("after", v2)):
        c = write_cnl(work, snp, f"{name}_{tag}", f0, esr, vals, None, span_mhz=3.0, step_khz=5.0)
        rr = wsl(f"{CRF} sparam {c.name} -o out_{name}_{tag}.s2p", work)
        if rr.returncode:
            raise RuntimeError(rr.stdout + rr.stderr)
        ff, s2 = read_crf_s2p(work / f"out_{name}_{tag}.s2p")
        out[f"sweep_{tag}"] = {"f_hz": ff.tolist(), "S11_db": (20 * np.log10(np.abs(s2[:, 0, 0]))).tolist(),
                               "S22_db": (20 * np.log10(np.abs(s2[:, 1, 1]))).tolist(),
                               "S21_db": (20 * np.log10(np.abs(s2[:, 1, 0]) + 1e-30)).tolist()}
        if tag == "after":
            py = python_two_port(ff, (f, S), r, v2, esr)
            out["crosscheck_max_abs_dS_vs_python"] = float(np.abs(py - s2).max())
    # optimizer trajectory of step 2 (from the naive start): the same seeded run stopped after 1, 2, ... iterations
    traj = []
    for k in range(1, iters + 1):
        c = write_cnl(work, sym, f"{name}_it", f0, esr, v0, free1, band_khz=band_khz, tie=tie,
                      optimize=f"algorithm=minimax maxiter={k} seed=1")
        ok = crf_opt(work, c)
        vk = best_vals(ok, v0)
        for k2, k1 in tie.items():
            vk[k2] = vk[k1]
        traj.append({"iter": k, "vals": vk, "cost": ok["result"]["optimize"]["bestCost"],
                     "goals": {g["name"]: g["value"] for g in ok["result"]["optimize"]["goals"]}})
        if ok["result"]["optimize"]["outcome"] == "goalsMet":
            break
    out["trajectory"] = traj
    (work / f"tune_{name}.json").write_text(json.dumps(out, indent=1))
    return out


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    a1 = sub.add_parser("assemble")
    a1.add_argument("--cases", required=True)
    a1.add_argument("--out", required=True)
    a1.add_argument("--rotation-average", action="store_true", help="also write <out>_sym with the C16-averaged matrix")
    a2 = sub.add_parser("check")
    a2.add_argument("--cases", required=True)
    a2.add_argument("--full", required=True)
    a3 = sub.add_parser("tune")
    a3.add_argument("--snp", required=True)
    a3.add_argument("--name", required=True)
    a3.add_argument("--esr", type=float, default=0.68)
    a3.add_argument("--band-khz", type=float, default=2.0)
    a3.add_argument("--c-eigen", type=float, default=1.8278810808108165e-12)
    a3.add_argument("--b0", type=float, default=3.0)
    a3.add_argument("--iters", type=int, default=0)
    a = ap.parse_args()
    if a.cmd == "tune":
        snp = Path(a.snp)
        res = tune(snp, snp.parent / "circuitrf", a.name, a.b0 * GAMMA_BAR, a.esr, a.band_khz, a.c_eigen, a.iters)
        print(json.dumps({k: res[k] for k in res if not k.startswith("sweep") and k != "trajectory"}, indent=1, default=str))
    elif a.cmd == "assemble":
        cases = {c: Path(f"{a.cases}_{c}") for c in ("em", "me", "mm", "ee")}
        f, S, meta = sparams.assemble(cases)
        r = json.loads((cases["em"] / "driven.json").read_text())["r_ref"]
        sparams.write_touchstone(Path(a.out), f, S, r, [f"birdcage 16-port (gaps), rebuilt from {a.cases}_{{em,me,mm,ee}}",
                                                       f"rung 0 at {np.degrees(meta['phi0']):.3f} deg, rung n at phi0 + n*22.5 deg"])
        if a.rotation_average:
            Ss = sparams.rotation_average(S)
            q = Path(a.out)
            sparams.write_touchstone(q.with_name(q.stem + "_sym" + q.suffix), f, Ss, r,
                                     [f"birdcage 16-port (gaps), rebuilt from {a.cases}_{{em,me,mm,ee}} and averaged over the 16 rotations"])
            print("rotation-averaged: max |S - S_sym| =", f"{np.abs(S - Ss).max():.2e}")
        U = np.abs(np.einsum("fij,fkj->fik", S, S.conj()) - np.eye(S.shape[1])).max()
        print(f"wrote {a.out}: {len(f)} frequencies, {S.shape[1]} ports, R {r:g} ohm; max |S S^H - I| = {U:.2e}")
    elif a.cmd == "check":
        cases = {c: Path(f"{a.cases}_{c}") for c in ("em", "me", "mm", "ee")}
        f, S, _ = sparams.assemble(cases)
        ff, Sf = sparams.read_port_s(Path(a.full))
        col = Sf[:, :, 0]
        d = np.abs(S[:, :, 0] - col)
        print(f"column 1, full model vs rebuilt from the eighth models: max |dS| = {d.max():.2e}, "
              f"|S11| full {np.abs(col[:, 0]).min():.4f}..{np.abs(col[:, 0]).max():.4f}")
        out = {"f_hz": f.tolist(), "max_abs_dS_col1": float(d.max()), "max_abs_dS_per_port": d.max(0).tolist(),
               "S_col1_full_re": col.real.tolist(), "S_col1_full_im": col.imag.tolist(),
               "S_col1_rebuilt_re": S[:, :, 0].real.tolist(), "S_col1_rebuilt_im": S[:, :, 0].imag.tolist()}
        (Path(a.full) / "assembly_check.json").write_text(json.dumps(out))


if __name__ == "__main__":
    main()
