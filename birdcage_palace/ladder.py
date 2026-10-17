"""Analytic ladder-network model of a birdcage coil.

The coil is treated as N meshes (one between each pair of neighbouring rungs). Mesh current J_n flows up rung n+1,
back along top ring segment n, down rung n and along bottom ring segment n. Because the network is periodic, its
eigenvectors are the discrete Fourier modes J_n = exp(i 2 pi k n / N), k = 0..N/2, and the mutual-inductance matrices
reduce to one number per mode (Leifer 1997):

    Lleg_k  = sum_m M_leg(0, m) cos(theta_k m)
    Lring_k = sum_m [M_tt(0, m) - M_tb(0, m)] cos(theta_k m)      (top and bottom ring segments carry opposite current)

with theta_k = 2 pi k / N and s_k = sin(theta_k / 2). Kirchhoff around one mesh then gives

    low-pass  (C in rungs):        w_k^2 = 1 / (Cl (Lleg_k + Lring_k / (2 s_k^2)))
    high-pass (C in ring segments): w_k^2 = 1 / (Cr (Lring_k + 2 s_k^2 Lleg_k))
    band-pass (both):               w_k^2 = (1/Cr + 2 s_k^2/Cl) / (Lring_k + 2 s_k^2 Lleg_k)

Rung current I_n = J_n - J_{n-1} is proportional to cos(2 pi k n / N + phi0), so k = 1 gives the sinusoidal current
distribution and the homogeneous transverse B1 (Hayes et al. 1985). k = 1..N/2-1 are degenerate pairs (cos and sin);
k = 0 (no rung current, counter-rotating ring currents) sits at DC for low-pass and is the highest mode for high-pass;
k = N/2 is single. The mesh model misses the co-rotating end-ring mode; for high/band-pass it is added separately as
w^2 = 1 / (Cr sum_m [M_tt + M_tb](0, m)).

Inductances come from a Neumann double integral over the conductor centre lines, regularised with the geometric mean
distance of a thin strip (0.2235 w) so the same kernel gives strip self-inductance (Grover) and mutuals. An RF shield of
radius Rs is represented by image conductors at Rs^2/R carrying the opposite current (exact for rungs in a long shield,
approximate for the rings). This is the same family of model as BirdcageBuilder (Chin et al. 2002); it is a
quasi-static first guess, the full-wave answer comes from Palace.
"""
from __future__ import annotations

import numpy as np

from .design import CoilDesign

MU0 = 4e-7 * np.pi
GMD_STRIP = 0.2235  # geometric mean distance of a thin strip, in units of its width


def _leg_points(R, phi, z0, z1, n=400):
    t = (np.arange(n) + 0.5) / n
    z = z0 + (z1 - z0) * t
    p = np.stack([np.full(n, R * np.cos(phi)), np.full(n, R * np.sin(phi)), z], 1)
    d = np.tile([0.0, 0.0, (z1 - z0) / n], (n, 1))
    return p, d


def _arc_points(R, phi0, phi1, z, n=160):
    t = (np.arange(n) + 0.5) / n
    ph = phi0 + (phi1 - phi0) * t
    p = np.stack([R * np.cos(ph), R * np.sin(ph), np.full(n, z)], 1)
    ds = R * (phi1 - phi0) / n
    d = np.stack([-np.sin(ph) * ds, np.cos(ph) * ds, np.zeros(n)], 1)
    return p, d


def _neumann(pa, da, pb, db, delta):
    r2 = ((pa[:, None, :] - pb[None, :, :]) ** 2).sum(-1)
    dot = da @ db.T
    return MU0 / (4 * np.pi) * float((dot / np.sqrt(r2 + delta ** 2)).sum())


def inductance_spectrum(des: CoilDesign) -> dict:
    """Return Lleg_k, Lring_k (k = 0..N/2) and the co-rotating ring inductance, in H."""
    N, R, l = des.legs, des.radius, des.length
    phis = 2 * np.pi * np.arange(N) / N
    dl, dr = GMD_STRIP * des.leg_width, GMD_STRIP * des.ring_width
    Rs = des.shield_diameter / 2 if des.shield_diameter else None
    Ri = Rs ** 2 / R if Rs else None

    p0, d0 = _leg_points(R, 0.0, -l / 2, l / 2)
    m_leg = np.zeros(N)
    for m in range(N):
        pm, dm = _leg_points(R, phis[m], -l / 2, l / 2)
        m_leg[m] = _neumann(p0, d0, pm, dm, dl if m == 0 else 0.0)
        if Ri:
            pi_, di_ = _leg_points(Ri, phis[m], -l / 2, l / 2)
            m_leg[m] -= _neumann(p0, d0, pi_, di_, 0.0)

    dphi = 2 * np.pi / N
    q0, e0 = _arc_points(R, 0.0, dphi, l / 2)
    m_tt, m_tb = np.zeros(N), np.zeros(N)
    for m in range(N):
        qm, em = _arc_points(R, phis[m], phis[m] + dphi, l / 2)
        qb, eb = _arc_points(R, phis[m], phis[m] + dphi, -l / 2)
        m_tt[m] = _neumann(q0, e0, qm, em, dr)
        m_tb[m] = _neumann(q0, e0, qb, eb, dr)
        if Ri:
            qi, ei = _arc_points(Ri, phis[m], phis[m] + dphi, l / 2)
            qj, ej = _arc_points(Ri, phis[m], phis[m] + dphi, -l / 2)
            m_tt[m] -= _neumann(q0, e0, qi, ei, 0.0)
            m_tb[m] -= _neumann(q0, e0, qj, ej, 0.0)

    ks = np.arange(N // 2 + 1)
    cosm = np.cos(2 * np.pi * np.outer(ks, np.arange(N)) / N)
    return {"k": ks, "Lleg": cosm @ m_leg, "Lring": cosm @ (m_tt - m_tb),
            "Lring_corot": float((m_tt + m_tb).sum()), "M_leg": m_leg, "M_tt": m_tt, "M_tb": m_tb}


def mode_frequencies(des: CoilDesign, ind: dict, c_leg: float | None, c_ring: float | None) -> dict:
    """Ladder-model resonances (Hz) for k = 0..N/2, plus the co-rotating ring mode where it exists."""
    N = des.legs
    k = ind["k"]
    s2 = np.sin(np.pi * k / N) ** 2
    Ll, Lr = ind["Lleg"], ind["Lring"]
    with np.errstate(divide="ignore", invalid="ignore"):
        if des.coil_type == "lowpass":
            w2 = 1.0 / (c_leg * (Ll + Lr / (2 * s2)))
            w2 = np.where(k == 0, 0.0, w2)
        elif des.coil_type == "highpass":
            w2 = 1.0 / (c_ring * (Lr + 2 * s2 * Ll))
        else:
            w2 = (1 / c_ring + 2 * s2 / c_leg) / (Lr + 2 * s2 * Ll)
    f = np.sqrt(w2) / (2 * np.pi)
    out = {"k": k.tolist(), "f_hz": f.tolist(), "degenerate": [bool(0 < kk < N // 2) for kk in k]}
    if des.coil_type != "lowpass":
        out["f_corot_hz"] = float(1 / (2 * np.pi * np.sqrt(c_ring * ind["Lring_corot"])))
    return out


def capacitor_for_target(des: CoilDesign, ind: dict, f_target: float, k: int = 1) -> tuple[float | None, float | None]:
    """Capacitance that puts mode k at f_target. Band-pass keeps c_ring = bandpass_ratio * c_leg."""
    N = des.legs
    w2 = (2 * np.pi * f_target) ** 2
    s2 = np.sin(np.pi * k / N) ** 2
    Ll, Lr = ind["Lleg"][k], ind["Lring"][k]
    if des.coil_type == "lowpass":
        return 1.0 / (w2 * (Ll + Lr / (2 * s2))), None
    if des.coil_type == "highpass":
        return None, 1.0 / (w2 * (Lr + 2 * s2 * Ll))
    a = des.bandpass_ratio  # Cr = a Cl  ->  w2 (Lr + 2 s2 Ll) = (1/a + 2 s2) / Cl
    cl = (1 / a + 2 * s2) / (w2 * (Lr + 2 * s2 * Ll))
    return cl, a * cl


def leg_currents(N: int, k: int, phase: float = 0.0) -> np.ndarray:
    return np.cos(2 * np.pi * k * np.arange(N) / N + phase)


def analytic(des: CoilDesign, f_target: float | None = None) -> dict:
    """Inductances, first-guess capacitors and the full ladder spectrum for a design."""
    f_target = f_target or des.f_larmor
    ind = inductance_spectrum(des)
    cl, cr = des.c_leg, des.c_ring
    if (des.coil_type in ("lowpass", "bandpass") and cl is None) or (des.coil_type in ("highpass", "bandpass") and cr is None):
        cl, cr = capacitor_for_target(des, ind, f_target)
    spec = mode_frequencies(des, ind, cl, cr)
    return {"f_target_hz": f_target, "c_leg": cl, "c_ring": cr, "spectrum": spec,
            "Lleg_k_nH": (ind["Lleg"] * 1e9).tolist(), "Lring_k_nH": (ind["Lring"] * 1e9).tolist(),
            "Lleg_self_nH": float(ind["M_leg"][0] * 1e9), "Lring_seg_self_nH": float(ind["M_tt"][0] * 1e9),
            "model": "ladder network, Neumann/GMD inductances" + (", shield by images" if des.shield_diameter else "")}
