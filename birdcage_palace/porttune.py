"""Port tuning: the multi-port S-matrix of the bare coil (every capacitor gap a port) from Palace, then the capacitors
and the quadrature drive network tuned in circuitRF (open-source RF circuit simulator, headless CLI).

    birdcage porttune assemble  --cases runs/pt_e --out runs/pt_e/birdcage16.s16p --rotation-average
    birdcage porttune resonance --snp runs/pt_e/birdcage16.s16p [--c 1.828e-12]
    birdcage porttune tune      --snp runs/pt_e/birdcage16.s16p --name unloaded --esr 0.68 --band-khz 2
    birdcage porttune verify    --snp runs/pt_e/birdcage16.s16p --ref-case runs/z_ref --tuned-case runs/z_tuned
    birdcage porttune check     --cases runs/ptc_e --full runs/ptc_full

The four cases <prefix>_em/_me/_mm/_ee are eighth (or quarter) port cases with the planes between rungs; the full
matrix is rebuilt from them exactly (sparams.py). The step-by-step guide is docs/circuitrf-tuning.md.

circuitRF is found from --crf, then the BIRDCAGE_CIRCUITRF environment variable, then `circuitRF` / `circuitrf` on
the PATH, then ~/.local/share/circuitRF/current/circuitRF (where its install.sh puts it). On Linux it runs natively.
On Windows it runs inside WSL: --wsl-distro or BIRDCAGE_WSL_DISTRO names the distribution ("auto" = the WSL default).
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import re
import shlex
import shutil
import subprocess
from pathlib import Path

import numpy as np

from . import sparams

GAMMA_BAR = 42.577478e6
F_REF = 3.0 * GAMMA_BAR  # the 3 T reference coil; frequency bands below are scaled from it
DEFAULT_CRF = "~/.local/share/circuitRF/current/circuitRF"
CRF_NAMES = ("circuitRF", "circuitrf")


# ------------------------------------------------------------------------------------------- running circuitRF
def wslpath(p: Path) -> str:
    p = Path(p).resolve()
    return f"/mnt/{p.drive.rstrip(':').lower()}" + p.as_posix()[2:]


class CircuitRF:
    """Where circuitRF is and how to call it.

    crf: path of the executable (None = BIRDCAGE_CIRCUITRF, else auto). distro: WSL distribution used on Windows
    (None = BIRDCAGE_WSL_DISTRO, else "auto" = the WSL default distribution). On any other system circuitRF runs
    natively and wsl.exe is never called."""

    def __init__(self, crf: str | None = None, distro: str | None = None, system: str | None = None,
                 environ: dict | None = None):
        env = os.environ if environ is None else environ
        self.crf = crf or env.get("BIRDCAGE_CIRCUITRF") or None
        self.distro = distro or env.get("BIRDCAGE_WSL_DISTRO") or "auto"
        self.system = system or platform.system()

    @property
    def via_wsl(self) -> bool:
        return self.system == "Windows"

    def native_exe(self) -> str:
        """Linux/macOS: the executable path (raises FileNotFoundError with the install hint if there is none)."""
        if self.crf:
            p = os.path.expanduser(self.crf)
            if shutil.which(p) or Path(p).is_file():
                return p
            raise FileNotFoundError(f"circuitRF not found at {self.crf} (--crf / BIRDCAGE_CIRCUITRF)")
        for n in CRF_NAMES:
            w = shutil.which(n)
            if w:
                return w
        p = os.path.expanduser(DEFAULT_CRF)
        if Path(p).is_file():
            return p
        raise FileNotFoundError("circuitRF not found: not on the PATH and not at " + DEFAULT_CRF +
                                ". Install it (docs/circuitrf-tuning.md, step 0) or pass --crf / set BIRDCAGE_CIRCUITRF.")

    def _wsl_exe_word(self) -> str:
        """The executable as a bash word, resolved inside WSL."""
        if self.crf:
            c = self.crf
            if c.startswith("~/"):
                return '"$HOME"/' + shlex.quote(c[2:])
            return shlex.quote(c)
        alts = " || ".join(f"command -v {n}" for n in CRF_NAMES)
        return f'"$({alts} || echo "$HOME/.local/share/circuitRF/current/circuitRF")"'

    def argv(self, args: list[str], cwd: Path) -> list[str]:
        if not self.via_wsl:
            return [self.native_exe()] + [str(a) for a in args]
        script = f"cd {shlex.quote(wslpath(cwd))} && exec {self._wsl_exe_word()} " + " ".join(shlex.quote(str(a)) for a in args)
        pre = ["wsl.exe"] + ([] if self.distro.lower() in ("auto", "default", "") else ["-d", self.distro])
        return pre + ["--exec", "bash", "-lc", script]

    def run(self, args: list[str], cwd: Path, timeout: float = 3600) -> subprocess.CompletedProcess:
        return subprocess.run(self.argv(args, cwd), cwd=None if self.via_wsl else str(cwd), capture_output=True,
                              text=True, encoding="utf-8", errors="replace", timeout=timeout, stdin=subprocess.DEVNULL)

    def version(self, cwd: Path = Path(".")) -> str:
        """`circuitRF --version` (with a timeout: in 1.0.9, --help and a bare `help` start the GUI and wait)."""
        r = self.run(["--version"], cwd, timeout=60)
        return r.stdout.strip() or r.stderr.strip()


def add_crf_options(p: argparse.ArgumentParser) -> None:
    p.add_argument("--crf", default=None,
                   help="circuitRF executable (env BIRDCAGE_CIRCUITRF; default: circuitRF/circuitrf on the PATH, "
                        f"else {DEFAULT_CRF})")
    p.add_argument("--wsl-distro", default=None,
                   help="Windows only: WSL distribution that runs circuitRF (env BIRDCAGE_WSL_DISTRO; default auto = the "
                        "WSL default). Ignored on Linux, where circuitRF runs natively.")


# ------------------------------------------------------------------------------------------- Touchstone
_UNIT = {"HZ": 1.0, "KHZ": 1e3, "MHZ": 1e6, "GHZ": 1e9, "THZ": 1e12}


def read_touchstone(path: Path) -> tuple[np.ndarray, np.ndarray, float]:
    """Touchstone v1 S-parameter file (any port count; RI, MA or DB) -> f (Hz), S (nf, P, P), reference R."""
    path = Path(path)
    m = re.search(r"\.s(\d+)p$", path.name, re.I)
    if not m:
        raise ValueError(f"{path}: not a .sNp file")
    P = int(m.group(1))
    unit, fmt, z0 = 1e9, "MA", 50.0
    nums: list[float] = []
    for line in path.read_text(errors="replace").splitlines():
        line = line.split("!", 1)[0].strip()
        if not line:
            continue
        if line.startswith("#"):
            tok = line[1:].upper().split()
            for i, t in enumerate(tok):
                if t in _UNIT:
                    unit = _UNIT[t]
                elif t in ("RI", "MA", "DB"):
                    fmt = t
                elif t == "R" and i + 1 < len(tok):
                    z0 = float(tok[i + 1])
            continue
        if line.startswith("["):
            raise ValueError(f"{path}: Touchstone v2 keywords are not supported")
        nums += [float(x) for x in line.split()]
    per = 1 + 2 * P * P
    if len(nums) % per:
        raise ValueError(f"{path}: {len(nums)} numbers is not a multiple of {per}")
    a = np.array(nums).reshape(-1, per)
    f = a[:, 0] * unit
    x, y = a[:, 1::2], a[:, 2::2]
    if fmt == "RI":
        v = x + 1j * y
    elif fmt == "MA":
        v = x * np.exp(1j * np.radians(y))
    else:
        v = 10 ** (x / 20) * np.exp(1j * np.radians(y))
    S = v.reshape(-1, P, P)
    if P == 2:  # two-port files list S11 S21 S12 S22
        S = S.transpose(0, 2, 1)
    return f, S, z0


def read_snp(path: Path):
    return read_touchstone(path)


def meta_path(snp: Path) -> Path:
    return Path(snp).with_suffix(".json")


def read_meta(snp: Path) -> dict:
    p = meta_path(snp)
    return json.loads(p.read_text()) if p.exists() else {}


# ------------------------------------------------------------------------------------------- network algebra
def y_matrix(S, r):
    I = np.eye(S.shape[-1])
    return np.linalg.solve(r * (I + S), (I - S))


def _spline(f, S, fi):
    from scipy.interpolate import CubicSpline
    return CubicSpline(f, S.real, axis=0)(fi) + 1j * CubicSpline(f, S.imag, axis=0)(fi)


def resonances(f, S, r, c, f_lo=None, f_hi=None, step=2e3):
    """Resonances of the lossless coil with every gap closed by c (scalar or per-gap array): frequencies where an
    eigenvalue of the gap susceptance matrix Im(Y + jwC) crosses zero. S is cubic-spline interpolated (the bare-coil
    S is smooth), crossings are located by linear interpolation on a fine grid. Default band: the whole file."""
    f_lo = f[0] if f_lo is None else max(f_lo, f[0])
    f_hi = f[-1] if f_hi is None else min(f_hi, f[-1])
    fi = np.arange(f_lo, f_hi, step)
    Y = y_matrix(_spline(f, S, fi), r)
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


def k1_frequency(f, S, r, c, f_lo=None, f_hi=None):
    """Lowest resonance (the k = 1 pair: mean of the two lowest crossings when both are found)."""
    rs = resonances(f, S, r, c, f_lo, f_hi)
    if not rs:
        raise ValueError("no resonance in the band")
    pair = rs[:2] if len(rs) > 1 and rs[1] / rs[0] - 1 < 0.01 else rs[:1]
    return float(np.mean(pair)), rs


def c_k1_estimate(f, S, r, f0) -> float:
    """Common gap C that makes the k = 1 harmonic resonant at f0: the bare-coil gap admittance projected on the
    rung-current pattern cos(phi_n), B_k1 = Im(v^H Y v) / v^H v, then C = -B_k1 / w0 (exact for the rotation-averaged
    matrix, whose eigenvectors are the harmonics)."""
    N = S.shape[1]
    Y = y_matrix(_spline(f, S, np.array([f0]))[0], r)
    v = np.cos(2 * np.pi * np.arange(N) / N)
    b = float(np.imag(v @ Y @ v) / (v @ v))
    return -b / (2 * np.pi * f0)


def c_for_frequency(f, S, r, f_target, c0=None):
    """Common gap capacitor that puts the k=1 resonance on f_target (secant on 1/f^2, linear in C)."""
    if c0 is None:
        c0 = c_k1_estimate(f, S, r, f_target)
    lo, hi = f_target * 110e6 / F_REF, f_target * 150e6 / F_REF  # 110-150 MHz at 3 T
    c = [c0, c0 * 1.02]
    fr = [k1_frequency(f, S, r, cc, lo, hi)[0] for cc in c]
    for _ in range(8):
        y = [1 / x ** 2 for x in fr]
        slope = (y[-1] - y[-2]) / (c[-1] - c[-2])
        cn = c[-1] + (1 / f_target ** 2 - y[-1]) / slope
        c.append(cn)
        fr.append(k1_frequency(f, S, r, cn, lo, hi)[0])
        if abs(fr[-1] / f_target - 1) < 1e-7:
            break
    return c[-1], fr[-1]


# ------------------------------------------------------------------------------------------- circuitRF netlist
def layout(N: int) -> tuple[tuple[int, int], tuple[int, ...]]:
    """Drive gaps and isolation-trim gaps for an N-rung coil (N a multiple of 4).

    The two drive gaps are rungs 0 and N/4 (90 degrees apart). The isolation trim sits on the rungs nearest the
    bisector of the two drives (45 degrees from both) and the opposite ones: a common change there couples the two
    k = 1 members, which cancels the residual coupling of the drive gaps through the other modes.
    16 rungs: rungs 2 and 10. 12 rungs: 45 degrees falls between rungs 1 and 2, so rungs 1, 2, 7 and 8."""
    if N % 4:
        raise ValueError("port tuning needs a multiple of 4 rungs")
    drive = (0, N // 4)
    b = N / 8
    near = (int(b),) if b == int(b) else (int(np.floor(b)), int(np.ceil(b)))
    iso = tuple(sorted(set(near) | {(i + N // 2) % N for i in near}))
    return drive, iso


# tuned-value ranges (pF) for a coil whose k = 1 resonance needs about 1.825 pF per gap; scaled by c_res / 1.825 pF
RANGES = {"Ct": (1.0, 3.0, "lin"), "Ci": (0.5, 3.5, "lin"), "Cd1": (0.5, 5.0, "lin"), "Cd2": (0.5, 5.0, "lin"),
          "Cm1": (0.05, 300.0, "log"), "Cm2": (0.05, 300.0, "log"), "Cp1": (2.0, 200.0, "log"), "Cp2": (2.0, 200.0, "log")}
C_REF = 1.825e-12
SCALED = ("Ct", "Ci", "Cd1", "Cd2", "Cp1", "Cp2")  # Cm (series to the 50 ohm port) keeps its absolute range
FEED = "divider"  # "series": series match capacitor straight across the gap; "divider": the gap capacitor is split
#                   into Cd (top) and a tap capacitor Cp to the feed, then a series Cm to the 50 ohm port


def range_scale(c_res: float) -> float:
    """c_res / 1.825 pF, rounded to 2 decimals (1.00 for the 16-rung reference coil)."""
    return round(c_res / C_REF, 2) or 0.01


def rung_var(i: int, drive=(0, 4), iso=(2, 10)) -> str:
    if i == drive[0]:
        return "Cd1"
    if i == drive[1]:
        return "Cd2"
    return "Ci" if (iso and i in iso) else "Ct"


def write_cnl(work: Path, snp: Path, name: str, f0: float, esr: float, vals: dict, tuned: list[str] | None = None,
              drive=(0, 4), band_khz: float = 100.0, n_ports=16, iso=(2, 10), span_mhz: float = 1.0,
              step_khz: float = 2.0, optimize: str = "algorithm=minimax maxiter=300 seed=1",
              tie: dict | None = None, feed: str | None = None, scale: float = 1.0) -> Path:
    """Netlist: the EM N-port, each gap closed by its capacitor (with ESR), the two drive gaps also fed from a 50 ohm
    port (divider: Cd on top, tap Cp, series Cm). vals: Ct (most rungs), Ci (isolation rungs), Cd1/Cd2 (drive gaps),
    Cp1/Cp2 (taps), Cm1/Cm2 (series to the ports), in F. tie = {"Cd2": "Cd1", ...} makes the two drive ports share
    one value. Goals at f0 +- band: |S11|, |S22| <= -20 dB and |S21| <= -20 dB. scale multiplies the C ranges."""
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
            s = scale if k in SCALED else 1.0
            L.append(f"tune {k} min={mn * s:g} pF max={mx * s:g} pF scale={sc} opt=1")
        L += [f"goal match1 = dB(SP1.S(1,1)) analysis=SP1 over=freq lo={lo:.4f} MHz hi={hi:.4f} MHz le -20",
              f"goal match2 = dB(SP1.S(2,2)) analysis=SP1 over=freq lo={lo:.4f} MHz hi={hi:.4f} MHz le -20",
              f"goal isol = dB(SP1.S(2,1)) analysis=SP1 over=freq lo={lo:.4f} MHz hi={hi:.4f} MHz le -20",
              f"optimize {optimize}"]
    p = work / f"tune_{name}.cnl"
    p.write_text("\n".join(L) + "\n", newline="\n")  # LF everywhere: the netlist is edited and run in Linux
    return p


def crf_opt(work: Path, cnl: Path, crf: CircuitRF | None = None) -> dict:
    crf = crf or CircuitRF()
    r = crf.run(["opt", cnl.name, "--json"], work)
    (work / f"{cnl.stem}.opt.json").write_text(r.stdout)
    if r.stderr:
        (work / f"{cnl.stem}.opt.err").write_text(r.stderr)
    try:
        return json.loads(r.stdout)
    except json.JSONDecodeError:
        raise RuntimeError(f"circuitRF opt {cnl.name} gave no JSON (exit {r.returncode}):\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")


def crf_sparam(work: Path, cnl: Path, out_name: str, sets: dict | None = None, crf: CircuitRF | None = None) -> Path:
    crf = crf or CircuitRF()
    args = ["sparam", cnl.name]
    for k, v in (sets or {}).items():
        args += ["--set", f"{k}={v * 1e12:.6f}pF"]
    r = crf.run(args + ["-o", out_name], work)
    if r.returncode:
        raise RuntimeError(r.stdout + r.stderr)
    return work / out_name


def _val_pf(txt: str) -> float:
    m = re.match(r"\s*([-+\d.eE]+)\s*([pfnu]?)F?", txt)
    v = float(m.group(1))
    return v * {"p": 1e-12, "f": 1e-15, "n": 1e-9, "u": 1e-6, "": 1.0}[m.group(2)]


def rung_caps(vals: dict, n=16, drive=(0, 4), iso=(2, 10), feed: str | None = None) -> np.ndarray:
    """Capacitance across each gap with the feed disconnected (divider: Cd in series with Cp)."""
    feed = feed or FEED
    C = np.array([vals[rung_var(i, drive, iso)] for i in range(n)])
    if feed == "divider":
        for k, i in enumerate(drive, start=1):
            C[i] = 1 / (1 / vals[f"Cd{k}"] + 1 / vals[f"Cp{k}"])
    return C


def terminated(Y, w, vals, esr, drive=(0, 4), iso=(2, 10), feed: str | None = None):
    """Nodal admittance of the coil (bare-coil gap admittance Y) with the capacitors and the drive network.
    Nodes: N gaps, then for each drive port its tap (divider) and its port node. Returns (Yn, port node indices)."""
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


def python_two_port(f, S, r, vals: dict, esr, drive=(0, 4), iso=(2, 10), feed: str | None = None):
    """Cross-check: terminate the EM N-port in python and return the 2-port S (50 ohm) at the drive ports."""
    fe, Se = S
    Y = y_matrix(_spline(fe, Se, f), r)
    out = np.zeros((len(f), 2, 2), complex)
    for k in range(len(f)):
        Yn, pn = terminated(Y[k], 2 * np.pi * f[k], vals, esr, drive, iso, feed)
        Z = np.linalg.inv(Yn)[np.ix_(pn, pn)]
        out[k] = np.linalg.solve(Z + 50 * np.eye(2), Z - 50 * np.eye(2))
    return out


def rung_currents(f0, S, r, vals, esr, e=(1, 1j), drive=(0, 4), iso=(2, 10), feed: str | None = None):
    """Gap (= rung) currents with the two ports driven by EMFs e behind 50 ohm (quadrature by default)."""
    fe, Se = S
    Y = y_matrix(_spline(fe, Se, np.array([f0])), r)[0]
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


def _opt_summary(o: dict) -> dict:
    return {k: o["result"]["optimize"][k] for k in ("outcome", "finishReason", "bestCost", "iterations", "evaluations")}


def tune(snp: Path, work: Path, name: str, f0: float, esr: float, band_khz: float, c_eigen: float,
         iters: int = 0, crf: CircuitRF | None = None) -> dict:
    """Port tuning in three steps:
      1. resonance only: every gap the same C; the S-matrix gives the k=1 frequencies at the reference C (eigenmode
         loop result, or the ladder first guess the coil was built with) and the C that puts the pair on f0;
      2. circuitRF on the rotation-averaged matrix (ideal N-fold symmetry): Ct, Ci (isolation rungs), Cd, Cp and Cm
         shared by both drive ports;
      3. circuitRF on the matrix as solved (the mesh splits the k=1 pair slightly): the drive ports get their own Cd,
         Cp, Cm (per-port trims), starting from step 2.
    Writes <work>/tune_<name>.json and <snp dir>/c_rungs_<name>.json (C across each gap, feeds disconnected)."""
    crf = crf or CircuitRF()
    snp, work = Path(snp), Path(work)
    work.mkdir(parents=True, exist_ok=True)
    f, S, r = read_snp(snp)
    N = S.shape[1]
    drive, iso = layout(N)
    sym = snp.with_name(snp.stem + "_sym" + snp.suffix)
    if not sym.exists():  # rotation average (the exact N-fold symmetry)
        sparams.write_touchstone(sym, f, sparams.rotation_average(S), r,
                                 [f"{N}-port averaged over the {N} rotations (written by porttune tune)"])
    fs, Ss, _ = read_snp(sym)
    out = {"snp": str(snp), "n_ports": N, "drive_gaps": list(drive), "isolation_gaps": list(iso), "f0_hz": f0,
           "esr_ohm": esr, "band_khz": band_khz, "c_eigen": c_eigen}
    lo, hi = f0 * 0.85, f0 * 1.15
    try:  # (the bands are those of the 3 T runs in the README, scaled with f0)
        out["k1_at_c_eigen_hz"] = resonances(f, S, r, c_eigen)[:2]
        out["k1_at_c_eigen_sym_hz"] = resonances(fs, Ss, r, c_eigen)[:2]
        c_res, f_res = c_for_frequency(f, S, r, f0, c_eigen if abs(c_eigen / c_k1_estimate(f, S, r, f0) - 1) < 0.02 else None)
        out["c_resonance_only"] = c_res
        out["k1_at_c_resonance_only_hz"] = resonances(f, S, r, c_res, lo, hi)[:2]
    except ValueError as e:  # lossy (loaded) matrix: no clean susceptance zero crossing
        out["resonance_note"] = str(e)
        c_res = c_k1_estimate(f, S, r, f0)
    sc = range_scale(c_res)
    out["range_scale"] = sc
    if FEED == "divider":  # gap = Cd in series with a large tap Cp (series value ~ c_res); Cm nearly a short
        cp0 = 20e-12 * sc
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
    kw = dict(drive=drive, iso=iso, n_ports=N, scale=sc)
    cnl = write_cnl(work, sym, f"{name}_sym", f0, esr, v0, free1, band_khz=band_khz, tie=tie, **kw)
    o1 = crf_opt(work, cnl, crf)
    v1 = best_vals(o1, v0)
    for k2, k1 in tie.items():
        v1[k2] = v1[k1]
    out["sym"] = {"vals": v1, "opt": _opt_summary(o1), "goals": o1["result"]["optimize"]["goals"]}
    tuned = free1 + list(tie)
    cnl2 = write_cnl(work, snp, name, f0, esr, v1, tuned, band_khz=band_khz, **kw)
    o2 = crf_opt(work, cnl2, crf)
    v2 = best_vals(o2, v1)
    out["raw"] = {"vals": v2, "opt": _opt_summary(o2), "goals": o2["result"]["optimize"]["goals"]}
    # wide sweeps (circuitRF) before and after; before = every gap at the reference C, drive network as tuned
    vb = dict(v2, Ct=c_eigen, Ci=c_eigen)
    if FEED == "divider":
        for k in (1, 2):
            vb[f"Cd{k}"] = 1 / (1 / c_eigen - 1 / v2[f"Cp{k}"])
    else:
        vb["Cd1"] = vb["Cd2"] = c_eigen
    tk = dict(drive=drive, iso=iso)
    for tag, vals in (("start", v0), ("sym", v1), ("raw", v2)):
        Iq = rung_currents(f0, (f, S), r, vals, esr, **tk)
        out.setdefault("k1_content_quadrature", {})[tag] = k1_content(Iq)
        out.setdefault("rung_current_max_over_min", {})[tag] = float(np.abs(Iq).max() / np.abs(Iq).min())
    span = max(3.0, 1.5 * abs(f0 - np.mean(out.get("k1_at_c_eigen_hz") or [f0])) / 1e6)
    for tag, vals in (("before", vb), ("after", v2)):
        c = write_cnl(work, snp, f"{name}_{tag}", f0, esr, vals, None, span_mhz=span if tag == "before" else 3.0,
                      step_khz=5.0 if tag == "after" or span <= 3.0 else 20.0, **kw)
        crf_sparam(work, c, f"out_{name}_{tag}.s2p", crf=crf)
        ff, s2, _ = read_touchstone(work / f"out_{name}_{tag}.s2p")
        out[f"sweep_{tag}"] = {"f_hz": ff.tolist(), "S11_db": (20 * np.log10(np.abs(s2[:, 0, 0]))).tolist(),
                               "S22_db": (20 * np.log10(np.abs(s2[:, 1, 1]))).tolist(),
                               "S21_db": (20 * np.log10(np.abs(s2[:, 1, 0]) + 1e-30)).tolist()}
        if tag == "after":
            py = python_two_port(ff, (f, S), r, v2, esr, **tk)
            out["crosscheck_max_abs_dS_vs_python"] = float(np.abs(py - s2).max())
            i0 = int(np.argmin(np.abs(ff - f0)))
            out["at_f0"] = {"f_hz": float(ff[i0]), **{k: float(out["sweep_after"][k][i0]) for k in ("S11_db", "S22_db", "S21_db")}}
    # optimizer trajectory of step 2 (from the naive start): the same seeded run stopped after 1, 2, ... iterations
    traj = []
    for k in range(1, iters + 1):
        c = write_cnl(work, sym, f"{name}_it", f0, esr, v0, free1, band_khz=band_khz, tie=tie,
                      optimize=f"algorithm=minimax maxiter={k} seed=1", **kw)
        ok = crf_opt(work, c, crf)
        vk = best_vals(ok, v0)
        for k2, k1 in tie.items():
            vk[k2] = vk[k1]
        traj.append({"iter": k, "vals": vk, "cost": ok["result"]["optimize"]["bestCost"],
                     "goals": {g["name"]: g["value"] for g in ok["result"]["optimize"]["goals"]}})
        if ok["result"]["optimize"]["outcome"] == "goalsMet":
            break
    out["trajectory"] = traj
    cr = rung_caps(v2, N, drive, iso)
    out["c_rungs"] = cr.tolist()
    try:
        out["k1_predicted_feeds_disconnected_hz"] = resonances(f, S, r, cr, lo, hi)[:2]
    except ValueError:
        pass
    (snp.parent / f"c_rungs_{name}.json").write_text(json.dumps(cr.tolist()))
    (work / f"tune_{name}.json").write_text(json.dumps(out, indent=1))
    return out


def summary(res: dict) -> str:
    """Plain-text result of a tune."""
    v = res["raw"]["vals"]
    pf = lambda x: f"{x * 1e12:.4f} pF"  # noqa: E731
    N, (d1, d2), iso = res["n_ports"], res["drive_gaps"], res["isolation_gaps"]
    L = [f"port tuning, {N} gaps, f0 {res['f0_hz'] / 1e6:.4f} MHz, ESR {res['esr_ohm']:g} ohm, band +-{res['band_khz']:g} kHz"]
    if "k1_at_c_eigen_hz" in res:
        L.append(f"  k=1 with every gap at {pf(res['c_eigen'])}: " + " / ".join(f"{x / 1e6:.3f}" for x in res["k1_at_c_eigen_hz"]) + " MHz")
    if "c_resonance_only" in res:
        L.append(f"  every gap equal, C for f0: {pf(res['c_resonance_only'])}")
    L += [f"  C on the other {N - 2 - len(iso)} rungs (Ct): {pf(v['Ct'])}",
          f"  C on rungs {', '.join(map(str, iso))} (Ci, isolation): {pf(v['Ci'])}",
          f"  drive gaps {d1} and {d2}: Cd {pf(v['Cd1'])} / {pf(v['Cd2'])}, Cp {pf(v['Cp1'])} / {pf(v['Cp2'])}, "
          f"Cm {pf(v['Cm1'])} / {pf(v['Cm2'])}"]
    a = res.get("at_f0")
    if a:
        L.append(f"  at {a['f_hz'] / 1e6:.4f} MHz: |S11| {a['S11_db']:.1f} dB, |S22| {a['S22_db']:.1f} dB, |S21| {a['S21_db']:.1f} dB")
    g = {x["name"]: x["value"] for x in res["raw"]["goals"]}
    L.append(f"  worst over the band: |S11| {g['match1']:.1f}, |S22| {g['match2']:.1f}, |S21| {g['isol']:.1f} dB; "
             f"circuitRF: {res['sym']['opt']['outcome']} / {res['raw']['opt']['outcome']}, evaluations "
             f"{res['sym']['opt']['evaluations']} + {res['raw']['opt']['evaluations']}")
    L.append(f"  k=1 share of the rung-current power, quadrature drive: {res['k1_content_quadrature']['raw']:.3f}; "
             f"circuitRF vs numpy termination max |dS| {res['crosscheck_max_abs_dS_vs_python']:.1e}")
    return "\n".join(L)


# ------------------------------------------------------------------------------------------- commands
def cmd_assemble(a) -> dict:
    from .design import CoilDesign
    cases = {c: Path(f"{a.cases}_{c}") for c in ("em", "me", "mm", "ee")}
    missing = [str(p) for p in cases.values() if not (p / "postpro" / "port-S.csv").exists()]
    if missing:
        raise SystemExit("not solved yet (no postpro/port-S.csv): " + ", ".join(missing))
    f, S, meta = sparams.assemble(cases)
    r = json.loads((cases["em"] / "driven.json").read_text())["r_ref"]
    N = S.shape[1]
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    sparams.write_touchstone(out, f, S, r, [f"birdcage {N}-port (gaps), rebuilt from {a.cases}_{{em,me,mm,ee}}",
                                           f"rung 0 at {np.degrees(meta['phi0']):.3f} deg, rung n at phi0 + n*{360 / N:g} deg"])
    msg = []
    if a.rotation_average:
        Ss = sparams.rotation_average(S)
        sparams.write_touchstone(out.with_name(out.stem + "_sym" + out.suffix), f, Ss, r,
                                 [f"birdcage {N}-port (gaps), rebuilt from {a.cases}_{{em,me,mm,ee}} and averaged over the {N} rotations"])
        msg.append(f"rotation-averaged: max |S - S_sym| = {np.abs(S - Ss).max():.2e}")
    U = float(np.abs(np.einsum("fij,fkj->fik", S, S.conj()) - np.eye(N)).max())
    des = CoilDesign.from_json(cases["em"] / "design.json")
    stats = {}
    for c, p in cases.items():
        st = p / "solve_stats.json"
        if st.exists():
            stats[c] = json.loads(st.read_text())
    m = {"n_ports": N, "phi0": meta["phi0"], "r_ref": r, "f0_hz": des.f_larmor, "b0": des.b0,
         "c_build": des.c_leg, "coil_type": des.coil_type, "cases": {c: str(p) for c, p in cases.items()},
         "mesh": {"h_coil": des.h_coil, "h_gap": des.h_gap, "h_max": des.h_max, "order": des.order},
         "f_min_hz": float(f[0]), "f_max_hz": float(f[-1]), "n_freq": len(f), "lossless_err": U, "solve_stats": stats}
    meta_path(out).write_text(json.dumps(m, indent=1))
    msg.append(f"wrote {out}: {len(f)} frequencies, {N} ports, R {r:g} ohm; max |S S^H - I| = {U:.2e}")
    print("\n".join(msg), flush=True)
    return m


def _f0_c(a, snp) -> tuple[float, float]:
    meta = read_meta(snp)
    f0 = a.b0 * GAMMA_BAR if a.b0 else meta.get("f0_hz", 3.0 * GAMMA_BAR)
    c = a.c_eigen if a.c_eigen else meta.get("c_build") or 1.8278810808108165e-12
    return f0, c


def cmd_resonance(a) -> dict:
    snp = Path(a.snp)
    f, S, r = read_snp(snp)
    f0, cb = _f0_c(argparse.Namespace(b0=a.b0, c_eigen=a.c), snp)
    out = {"snp": str(snp), "f0_hz": f0, "c_ref": cb}
    try:
        out["k1_at_c_ref_hz"] = resonances(f, S, r, cb)[:2]
    except ValueError as e:
        out["k1_at_c_ref_note"] = str(e)
    c_res, f_res = c_for_frequency(f, S, r, f0)
    out["c_for_f0"] = c_res
    out["k1_at_c_for_f0_hz"] = resonances(f, S, r, c_res, 0.95 * f0, 1.05 * f0)[:2]
    out["all_resonances_at_c_for_f0_mhz"] = [x / 1e6 for x in resonances(f, S, r, c_res)]
    U = float(np.abs(np.einsum("fij,fkj->fik", S, S.conj()) - np.eye(S.shape[1])).max())
    out["lossless_err"] = U
    k = out.get("k1_at_c_ref_hz")
    print(f"{snp.name}: {S.shape[1]} ports, {len(f)} frequencies {f[0] / 1e6:.2f}-{f[-1] / 1e6:.2f} MHz, max |S S^H - I| {U:.1e}")
    if k:
        print(f"  k=1 with every gap at {cb * 1e12:.4f} pF: " + " / ".join(f"{x / 1e6:.3f}" for x in k) +
              f" MHz ({(np.mean(k) / f0 - 1) * 100:+.2f} % from f0)")
    print(f"  every gap equal, C for k=1 at {f0 / 1e6:.4f} MHz: {c_res * 1e12:.4f} pF -> " +
          " / ".join(f"{x / 1e6:.4f}" for x in out["k1_at_c_for_f0_hz"]) + " MHz")
    print("  resonances in the sweep band at that C (MHz): " + ", ".join(f"{x:.2f}" for x in out["all_resonances_at_c_for_f0_mhz"]))
    o = snp.with_name(snp.stem + "_resonance.json")
    o.write_text(json.dumps(out, indent=1))
    return out


def cmd_tune(a) -> dict:
    snp = Path(a.snp)
    f0, c = _f0_c(a, snp)
    crf = CircuitRF(a.crf, a.wsl_distro)
    work = Path(a.work) if a.work else snp.parent / "circuitrf"
    res = tune(snp, work, a.name, f0, a.esr, a.band_khz, c, a.iters, crf)
    print(summary(res), flush=True)
    print(f"wrote {work / f'tune_{a.name}.json'}, {work / f'tune_{a.name}.cnl'} (the tuned netlist) and "
          f"{snp.parent / f'c_rungs_{a.name}.json'}")
    return res


def _pair(case: Path) -> list[float]:
    m = json.loads((Path(case) / "modes.json").read_text())
    hp = m.get("homogeneous_pair")
    if not hp:
        raise SystemExit(f"{case}: no k=1 pair in modes.json")
    return sorted(hp["f_mhz"]) if isinstance(hp["f_mhz"], list) else [hp["f_mhz"]]


def cmd_verify(a) -> dict:
    """Palace eigen solves of the z-half model (all rungs) at the reference C and at the tuned set, against the
    resonances the N-port predicts for the same capacitors. The meshes differ, so compare the shifts."""
    from .design import CoilDesign
    snp = Path(a.snp)
    f, S, r = read_snp(snp)
    ref, tun = Path(a.ref_case), Path(a.tuned_case)
    dref, dtun = CoilDesign.from_json(ref / "design.json"), CoilDesign.from_json(tun / "design.json")
    cref = np.array(dref.c_leg_rungs) if dref.c_leg_rungs else np.full(S.shape[1], dref.c_leg)
    ctun = np.array(dtun.c_leg_rungs) if dtun.c_leg_rungs else np.full(S.shape[1], dtun.c_leg)
    pr = [x / 1e6 for x in resonances(f, S, r, cref)[:2]]
    pt = [x / 1e6 for x in resonances(f, S, r, ctun)[:2]]
    er, et = _pair(ref), _pair(tun)
    out = {"palace_ref_mhz": er, "palace_tuned_mhz": et, "nport_ref_mhz": pr, "nport_tuned_mhz": pt,
           "shift_palace_mhz": float(np.mean(et) - np.mean(er)), "shift_nport_mhz": float(np.mean(pt) - np.mean(pr)),
           "mesh_offset_rel": float(np.mean(er) / np.mean(pr) - 1)}
    out["shift_diff_khz"] = (out["shift_palace_mhz"] - out["shift_nport_mhz"]) * 1e3
    fmt = lambda v: " / ".join(f"{x:.3f}" for x in v)  # noqa: E731
    print(f"{'':28} {'Palace (z-half)':>20} {'N-port prediction':>20}")
    print(f"{'reference C':28} {fmt(er):>20} {fmt(pr):>20}")
    print(f"{'tuned set (feeds open)':28} {fmt(et):>20} {fmt(pt):>20}")
    print(f"shift: Palace {out['shift_palace_mhz']:+.4f} MHz, N-port {out['shift_nport_mhz']:+.4f} MHz "
          f"(difference {out['shift_diff_khz']:+.1f} kHz); mesh offset at the reference C {out['mesh_offset_rel'] * 100:+.3f} %")
    if a.out:
        Path(a.out).write_text(json.dumps(out, indent=1))
    return out


def cmd_check(a) -> dict:
    cases = {c: Path(f"{a.cases}_{c}") for c in ("em", "me", "mm", "ee")}
    f, S, _ = sparams.assemble(cases)
    ff, Sf = sparams.read_port_s(Path(a.full))
    col = Sf[:, :, 0]
    d = np.abs(S[:, :, 0] - col)
    print(f"column 1, full model vs rebuilt from the reduced models: max |dS| = {d.max():.2e}, "
          f"|S11| full {np.abs(col[:, 0]).min():.4f}..{np.abs(col[:, 0]).max():.4f}")
    out = {"f_hz": f.tolist(), "max_abs_dS_col1": float(d.max()), "max_abs_dS_per_port": d.max(0).tolist(),
           "S_col1_full_re": col.real.tolist(), "S_col1_full_im": col.imag.tolist(),
           "S_col1_rebuilt_re": S[:, :, 0].real.tolist(), "S_col1_rebuilt_im": S[:, :, 0].imag.tolist()}
    (Path(a.full) / "assembly_check.json").write_text(json.dumps(out))
    return out


def add_parser(sub) -> argparse.ArgumentParser:
    """`porttune` and its sub-commands, on an argparse sub-parser collection."""
    pt = sub.add_parser("porttune", help="port tuning with circuitRF: assemble the N-port, check the resonance, tune, verify",
                        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ps = pt.add_subparsers(dest="pt_cmd", required=True)
    a1 = ps.add_parser("assemble", help="rebuild the full N-port from the four PEC/PMC port cases")
    a1.add_argument("--cases", required=True, help="prefix: <prefix>_em, _me, _mm, _ee")
    a1.add_argument("--out", required=True, help="Touchstone file, e.g. runs/x/birdcage12.s12p")
    a1.add_argument("--rotation-average", action="store_true", help="also write <out>_sym with the rotation-averaged matrix")
    a1.set_defaults(fn=cmd_assemble)
    a2 = ps.add_parser("resonance", help="k=1 frequency with every gap closed by one C, and the C that puts it on f0")
    a2.add_argument("--snp", required=True)
    a2.add_argument("--c", type=float, default=None, help="F; default the C the coil was built with (from <snp>.json)")
    a2.add_argument("--b0", type=float, default=None, help="T; default from <snp>.json, else 3")
    a2.set_defaults(fn=cmd_resonance)
    a3 = ps.add_parser("tune", help="circuitRF: capacitors + capacitive-divider feed against the match/isolation goals")
    a3.add_argument("--snp", required=True)
    a3.add_argument("--name", required=True, help="tag for the output files (tune_<name>.json, tune_<name>.cnl, ...)")
    a3.add_argument("--esr", type=float, default=0.68, help="ohm in series with every capacitor (0.68: Q ~1000 at 128 MHz)")
    a3.add_argument("--band-khz", type=float, default=2.0, help="goals hold over f0 +- this band")
    a3.add_argument("--c-eigen", "--c-ref", dest="c_eigen", type=float, default=None,
                    help="F, reference C for the 'before' sweep (eigenmode-loop result or ladder guess); default from <snp>.json")
    a3.add_argument("--b0", type=float, default=None, help="T; default from <snp>.json, else 3")
    a3.add_argument("--iters", type=int, default=0, help="also record the optimiser trajectory, 1..iters iterations")
    a3.add_argument("--work", default=None, help="circuitRF working folder (default <snp dir>/circuitrf)")
    add_crf_options(a3)
    a3.set_defaults(fn=cmd_tune)
    a4 = ps.add_parser("verify", help="Palace z-half eigen solves vs the N-port prediction for the same capacitors")
    a4.add_argument("--snp", required=True)
    a4.add_argument("--ref-case", required=True, help="z-half case solved at one C on every gap (with modes.json)")
    a4.add_argument("--tuned-case", required=True, help="the same mesh with the tuned set (variant --c-rungs)")
    a4.add_argument("--out", default=None)
    a4.set_defaults(fn=cmd_verify)
    a5 = ps.add_parser("check", help="rebuilt first column vs a full-model driven solve")
    a5.add_argument("--cases", required=True)
    a5.add_argument("--full", required=True)
    a5.set_defaults(fn=cmd_check)
    return pt


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(prog="port_tune")
    sub = ap.add_subparsers(dest="cmd", required=True)
    pt = add_parser(sub)
    # tools/port_tune.py <sub-command> ...: the porttune sub-commands at top level
    a = pt.parse_args(argv)
    a.fn(a)
