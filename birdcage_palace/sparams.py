"""Multi-port S-matrix of the bare coil (every capacitor gap a port) from Palace driven sweeps.

A full model gives the N x N matrix directly. With symmetry planes it is rebuilt exactly from the reduced models:
the coil is invariant under the mirrors x -> -x and y -> -y, so S commutes with the port permutations they induce.
A quarter (or eighth) model with plane conditions c = (PEC/PMC on x = 0, on y = 0) computes the response to the
symmetry-adapted excitation a = sum_g chi_c(g) e_{g(j)} (chi = -1 for each mirror across a PEC plane, +1 across a PMC
plane, the parity of the rung current), i.e. S^c_ij = sum_g chi_c(g) S_{i, g(j)}. Character orthogonality inverts it:

    S_{g_p(i), g_q(j)} = 1/4 sum_c chi_c(g_p) chi_c(g_q) S^c_ij

so the four PEC/PMC combinations together give every entry. The z = 0 plane (eighth model) is PEC in all four: the
z-even (PMC) family has no voltage across the gaps and does not couple to the ports.
"""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import numpy as np


def read_port_s(case: Path) -> tuple[np.ndarray, np.ndarray]:
    """f (Hz) and S (nf, P, P) from postpro/port-S.csv (|S| in dB, arg in degrees)."""
    p = Path(case) / "postpro" / "port-S.csv"
    with open(p, newline="") as fh:
        rows = [r for r in csv.reader(fh) if r]
    head = [h.strip() for h in rows[0]]
    data = np.array([[float(x) for x in r] for r in rows[1:]])
    f = data[:, 0] * 1e9
    idx = {}
    for j, h in enumerate(head):
        m = re.match(r"(\|S|arg\(S)\[(\d+)\]\[(\d+)\]", h)
        if m:
            idx[(m.group(1)[0], int(m.group(2)), int(m.group(3)))] = j
    P = max(k[1] for k in idx)
    S = np.zeros((len(f), P, P), complex)
    for o in range(1, P + 1):
        for e in range(1, P + 1):
            if ("|", o, e) in idx:
                mag = 10 ** (data[:, idx[("|", o, e)]] / 20)
                ph = np.radians(data[:, idx[("a", o, e)]])
                S[:, o - 1, e - 1] = mag * np.exp(1j * ph)
    # sort / unique in frequency (adaptive sweeps print in sample order)
    order = np.argsort(f)
    f, S = f[order], S[order]
    keep = np.r_[True, np.diff(f) > 0]
    return f[keep], S[keep]


def _images(N: int, phi0: float, port_phi: list[float], planes: list[str]):
    """For each full rung n: (model port index, flips applied (sx, sy) as 0/1)."""
    out = []
    pphi = np.array(port_phi)
    for n in range(N):
        ph = phi0 + 2 * np.pi * n / N
        c = [np.cos(ph), np.sin(ph)]
        fl = [0, 0]
        for i, a in enumerate("xy"):
            if a in planes and c[i] < 0:
                c[i], fl[i] = -c[i], 1
            elif a in planes and abs(c[i]) < 1e-9:
                raise ValueError("assembly needs the planes between rungs (--sym-cut between)")
        rp = np.arctan2(c[1], c[0])
        d = np.abs(np.angle(np.exp(1j * (pphi - rp))))
        j = int(np.argmin(d))
        assert d[j] < 1e-6
        out.append((j, fl))
    return out


def assemble(cases: dict[str, Path]) -> tuple[np.ndarray, np.ndarray, dict]:
    """cases: {'em': path, 'me': ..., 'mm': ..., 'ee': ...} (first letter plane x = 0, second y = 0; e = PEC, m = PMC).
    Returns f, the full N x N S-matrix in rung order n = 0..N-1, and meta."""
    first = Path(next(iter(cases.values())))
    des = json.loads((first / "design.json").read_text())
    info = json.loads((first / "build.json").read_text())
    N, phi0 = des["legs"], info["phi0"]
    planes = [p for p in info["planes"] if p in "xy"]
    if planes != ["x", "y"]:
        raise ValueError("assembly expects a quarter or eighth model (planes x = 0 and y = 0)")
    im = _images(N, phi0, [p["phi"] for p in info["ports"]], planes)
    data = {}
    f0 = None
    for c, case in cases.items():
        f, S = read_port_s(case)
        if f0 is None:
            f0 = f
        elif len(f) != len(f0) or np.abs(f - f0).max() > 1:
            raise ValueError("frequency grids differ between the combinations")
        data[c] = S
    chi = {c: (-1.0 if c[0] == "e" else 1.0, -1.0 if c[1] == "e" else 1.0) for c in cases}
    if len(chi) != 4:
        raise ValueError("need all four PEC/PMC combinations")
    S = np.zeros((len(f0), N, N), complex)
    for p in range(N):
        i, fp = im[p]
        for q in range(N):
            j, fq = im[q]
            acc = 0
            for c, (sx, sy) in chi.items():
                acc = acc + (sx ** (fp[0] + fq[0])) * (sy ** (fp[1] + fq[1])) * data[c][:, i, j]
            S[:, p, q] = acc / 4
    return f0, S, {"N": N, "phi0": phi0, "rung_phi": [phi0 + 2 * np.pi * n / N for n in range(N)]}


def write_touchstone(path: Path, f: np.ndarray, S: np.ndarray, r_ref: float, comments: list[str] = ()) -> Path:
    """Touchstone v1 (RI, Hz), any port count; four complex values per line."""
    P = S.shape[1]
    lines = [f"! {c}" for c in comments] + [f"# HZ S RI R {r_ref:g}"]
    for k, fk in enumerate(f):
        vals = []
        for i in range(P):
            row = S[k, i, :]
            chunks = [row[m:m + 4] for m in range(0, P, 4)]
            for ci, ch in enumerate(chunks):
                txt = " ".join(f"{v.real:+.10e} {v.imag:+.10e}" for v in ch)
                prefix = f"{fk:.6f} " if (i == 0 and ci == 0) else "  "
                vals.append(prefix + txt)
        lines += vals
    Path(path).write_text("\n".join(lines) + "\n")
    return Path(path)


def s_to_z(S: np.ndarray, r_ref: float) -> np.ndarray:
    I = np.eye(S.shape[-1])
    return r_ref * np.linalg.solve(I - S, I + S)


def terminate_caps_resonance(f: np.ndarray, S: np.ndarray, r_ref: float, c: np.ndarray, esr: float = 0.0):
    """Frequencies where the coil with every gap closed by capacitor c (F, per gap) resonates: the gap admittance
    matrix Y(f) + Y_C(f) is singular. Returns the smallest singular value of (Y + Y_C) normalised, per frequency."""
    I = np.eye(S.shape[-1])
    Y = np.linalg.solve(r_ref * (I + S), (I - S))  # Y = (1/R)(I+S)^-1 (I-S)
    w = 2 * np.pi * f
    zc = esr + 1 / (1j * w[:, None] * np.atleast_1d(c)[None, :])
    YT = Y + np.einsum("fi,ij->fij", 1 / zc, I)
    sv = np.linalg.svd(YT, compute_uv=False)
    return sv[:, -1] / sv[:, 0]


def rotation_average(S: np.ndarray) -> np.ndarray:
    """Project S onto matrices that commute with the cyclic rotation of the N rungs (S[p, q] = S[p+1, q+1]).

    The exact coil has the 16-fold rotation symmetry; a finite-element mesh does not, which splits each degenerate
    pair slightly. Averaging over the N rotations removes that discretisation asymmetry (the mirror symmetries are
    already exact in an assembled matrix)."""
    N = S.shape[-1]
    out = np.zeros_like(S)
    idx = np.arange(N)
    for r in range(N):
        P = (idx + r) % N
        out += S[:, P][:, :, P]
    return out / N
