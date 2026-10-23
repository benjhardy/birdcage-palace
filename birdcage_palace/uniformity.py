"""Quasi-static field uniformity of each ladder mode (Biot-Savart), the model behind the web builder's spectrum colours.

The mode currents follow from Kirchhoff on the ladder: mode k has mesh currents M_n = sin(k (phi_n + pi/N)) /
(2 sin(k pi/N)) in loop n (between rung n and n+1). The top end-ring segment n carries +M_n, the bottom one -M_n
(counter-rotating), and rung n carries M_n - M_(n-1) = cos(k phi_n). k = 0 is the uniform counter-rotating ring
mode (no rung current); the co-rotating ring mode has +M_n in both rings. Rungs are filaments on the strip
centrelines from ring centreline to ring centreline; ring segments are 4 straight chords each. With a shield, every
element gets the image the ladder model uses (radius Rs^2/R, opposite current), and with a closed PEC can both are
mirrored once in each end cap (z mirrored, opposite current: tangential components reverse, axial ones stay).

Metric, as birdcage_palace.modes on Palace's probes: over the ROI (r < 0.6 R, |z| < length/6; the Palace grid of
12 mm, coarsened to 0.06 R for coils larger than R = 20 cm) the coefficient of variation std/mean of the transverse
|B| = sqrt(Bx^2 + By^2), and the transverse fraction sum |B_t|^2 / sum |B|^2.

The web page (web/birdcage-builder/index.html, `Uniformity` in <script id="core">) is a line-by-line port;
tests/test_uniformity.py checks the two against each other and against the Palace results in results/*.json.
"""
from __future__ import annotations

import math

import numpy as np

MU0_4PI = 1e-7
N_CHORD = 4      # straight chords per end-ring segment
H_ROI = 0.012    # m, probe spacing of birdcage_palace.palace.probe_points
NUDGE = 1.1e-5   # m, the probe offset palace.probe_points applies


def _get(des, name, default=None):
    return des.get(name, default) if isinstance(des, dict) else getattr(des, name, default)


def roi_points(des) -> np.ndarray:
    R, l = _get(des, "diameter") / 2, _get(des, "length")
    h = max(H_ROI, 0.06 * R)
    # the grid of palace.probe_points, written as start + i h (numpy's arange rounds differently) and with a
    # 1e-9 relative margin on the circle, so the JS port picks exactly the same points
    gr = [-0.6 * R + i * h for i in range(math.ceil((1.2 * R + 1e-9) / h))]
    gz = [-l / 6 + i * h for i in range(math.ceil((l / 3 + 1e-9) / h))]
    r2 = (0.6 * R) ** 2 * (1 + 1e-9)
    pts = [(x, y, z) for z in gz for x in gr for y in gr if x * x + y * y <= r2]
    return np.array(pts) + NUDGE


def unit_sources(des) -> list[list[tuple]]:
    """3N unit sources: rung n, top ring segment n, bottom ring segment n (unit current, +z / +phi), each a list of
    straight segments (a, b, sign) including the shield and end-cap images."""
    N, R, l = _get(des, "legs"), _get(des, "diameter") / 2, _get(des, "length")
    Rs = _get(des, "shield_diameter")
    ends = _get(des, "shield_ends", "pec") or "pec"
    Ls = (_get(des, "shield_length") or l + 0.20) if Rs else None
    ph = [2 * math.pi * n / N for n in range(N)]
    base = []
    for n in range(N):
        c, s = math.cos(ph[n]), math.sin(ph[n])
        base.append([((R * c, R * s, -l / 2), (R * c, R * s, l / 2), 1.0)])
    for zc in (l / 2, -l / 2):
        for n in range(N):
            seg = []
            for j in range(N_CHORD):
                t0 = ph[n] + 2 * math.pi / N * j / N_CHORD
                t1 = ph[n] + 2 * math.pi / N * (j + 1) / N_CHORD
                seg.append(((R * math.cos(t0), R * math.sin(t0), zc), (R * math.cos(t1), R * math.sin(t1), zc), 1.0))
            base.append(seg)
    if not Rs:
        return base
    q = (Rs / 2) ** 2 / R / R  # image radius / R
    out = []
    for src in base:
        img = [((a[0] * q, a[1] * q, a[2]), (b[0] * q, b[1] * q, b[2]), -w) for a, b, w in src]
        both = src + img
        if ends == "pec":
            for zc in (Ls / 2, -Ls / 2):
                both = both + [((a[0], a[1], 2 * zc - a[2]), (b[0], b[1], 2 * zc - b[2]), -w) for a, b, w in src + img]
        out.append(both)
    return out


def segment_field(P: np.ndarray, a, b) -> np.ndarray:
    """B of a unit current on the straight segment a -> b at points P (n, 3)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = b - a
    t = d / np.linalg.norm(d)
    r1, r2 = P - a, P - b
    c = np.cross(t, r1)
    rho2 = (c * c).sum(1)
    f = (r1 @ t) / np.linalg.norm(r1, axis=1) - (r2 @ t) / np.linalg.norm(r2, axis=1)
    return MU0_4PI * c * (f / rho2)[:, None]


def unit_fields(des, P: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    P = roi_points(des) if P is None else P
    F = np.zeros((3 * _get(des, "legs"), len(P), 3))
    for i, src in enumerate(unit_sources(des)):
        for a, b, w in src:
            F[i] += w * segment_field(P, a, b)
    return P, F


def mode_weights(N: int, k: int, family: str = "ladder") -> np.ndarray:
    """Unit-source weights [rungs, top ring, bottom ring] of ladder mode k, or of the co-rotating ring mode."""
    if family == "corot":
        return np.concatenate([np.zeros(N), np.ones(N), np.ones(N)])
    ph = 2 * math.pi * np.arange(N) / N
    M = np.ones(N) if k == 0 else np.sin(k * (ph + math.pi / N)) / (2 * math.sin(k * math.pi / N))
    return np.concatenate([M - np.roll(M, 1), M, -M])


def metric(B: np.ndarray) -> dict:
    bt = np.hypot(B[:, 0], B[:, 1])
    return {"bt_cv": float(bt.std() / bt.mean()), "transverse_frac": float((bt ** 2).sum() / (B ** 2).sum())}


def spectrum_uniformity(des) -> dict:
    """Per-mode metric for k = 0..N/2 (k = 0: counter-rotating ring mode) and the co-rotating ring mode."""
    N = _get(des, "legs")
    _, F = unit_fields(des)
    modes = [dict(k=k, **metric(np.tensordot(mode_weights(N, k), F, 1))) for k in range(N // 2 + 1)]
    return {"modes": modes, "corot": metric(np.tensordot(mode_weights(N, 0, "corot"), F, 1)), "n_points": F.shape[1]}
