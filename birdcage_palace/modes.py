"""Read Palace eigenmode output and identify birdcage modes.

For every eigenmode:
  * capacitor currents (port-I.csv) -> rung (or ring-segment) current pattern -> spatial DFT -> mode index k and purity
  * probe B field in the ROI (r < 0.6 R, |z| < length/6) -> transverse B1 homogeneity (std/mean of |B_t|) and the
    transverse fraction sum|B_t|^2 / sum|B|^2
The homogeneous mode is the k = 1 mode with a transverse field (fraction > 0.9) and the lowest ROI std/mean; its
near-degenerate partner is the other k = 1 mode. Combining the pair with a 90 degree phase gives the circularly
polarised B1+ = (Bx + i By)/2 of an ideal quadrature drive.
"""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import numpy as np


def _read_csv(path: Path) -> tuple[list[str], np.ndarray]:
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh))
    head = [h.strip() for h in rows[0]]
    data = np.array([[float(x) for x in r] for r in rows[1:] if r], dtype=float)
    return head, data


def _complex_columns(head, data, pattern: str) -> dict:
    """Map 'key' -> complex column for headers like Re{V[3]} (V) / Im{V[3]} (V), or real-only headers."""
    out: dict = {}
    rx = re.compile(r"^(Re|Im)\{" + pattern + r"\}")
    rr = re.compile(r"^" + pattern + r" \(")
    for j, h in enumerate(head):
        m = rx.match(h)
        if m:
            key = m.groups()[1:]
            out.setdefault(key, np.zeros(len(data), complex))
            out[key] += data[:, j] * (1 if m.group(1) == "Re" else 1j)
            continue
        m = rr.match(h)
        if m:
            out.setdefault(m.groups(), np.zeros(len(data), complex))
            out[m.groups()] += data[:, j]
    return out


def load(case: Path) -> dict:
    case = Path(case)
    post = case / "postpro"
    head, eig = _read_csv(post / "eig.csv")
    col = {h.split(" ")[0]: j for j, h in enumerate(head)}
    f_re = eig[:, col.get("Re{f}", 1)]
    f_im = eig[:, col["Im{f}"]] if "Im{f}" in col else np.zeros_like(f_re)
    q = eig[:, col["Q"]] if "Q" in col else np.full_like(f_re, np.inf)
    err = eig[:, col["Error"]] if "Error" in col else None
    res = {"f_ghz": f_re, "f_im_ghz": f_im, "Q": q, "err": err}
    for name, sym in (("port-I.csv", "I"), ("port-V.csv", "V")):
        p = post / name
        if p.exists():
            h, d = _read_csv(p)
            cols = _complex_columns(h, d, sym + r"\[(\d+)\]")
            idx = sorted(int(k[0]) for k in cols)
            res[sym] = np.stack([cols[(str(i),)] for i in idx], 1)
    pb = post / "probe-B.csv"
    if pb.exists():
        h, d = _read_csv(pb)
        cols = _complex_columns(h, d, r"B_([xyz])\[(\d+)\]")
        n = max(int(k[1]) for k in cols)
        B = np.zeros((len(d), n, 3), complex)
        for (ax, i), v in cols.items():
            B[:, int(i) - 1, "xyz".index(ax)] = v
        res["B"] = B
    return res


def _dft_power(x: np.ndarray) -> np.ndarray:
    """|DFT|^2 folded to k = 0..N/2 (k and N-k summed), normalised to 1."""
    N = len(x)
    c = np.abs(np.fft.fft(x)) ** 2
    p = np.zeros(N // 2 + 1)
    for k in range(N):
        p[min(k, N - k)] += c[k]
    return p / p.sum() if p.sum() > 0 else p


def _align(v: np.ndarray) -> np.ndarray:
    """Remove the arbitrary complex phase of an eigenvector (make its largest entry real positive)."""
    j = np.argmax(np.abs(v))
    return v * np.exp(-1j * np.angle(v[j])) if np.abs(v[j]) > 0 else v


def analyse(case: Path) -> dict:
    case = Path(case)
    des = json.loads((case / "design.json").read_text())
    info = json.loads((case / "build.json").read_text())
    pr = json.loads((case / "probes.json").read_text())
    r = load(case)
    N = des["legs"]
    ports = info["ports"]
    n_ax, n_cor = pr["axial"], pr["coronal"]
    sl = {"axial": slice(0, n_ax), "coronal": slice(n_ax, n_ax + n_cor), "roi": slice(n_ax + n_cor, None)}
    leg_idx = [i for i, p in enumerate(ports) if p["kind"] == "leg"]
    top_idx = [i for i, p in enumerate(ports) if p["kind"] == "ring_top"]
    bot_idx = [i for i, p in enumerate(ports) if p["kind"] == "ring_bot"]
    modes = []
    for m in range(len(r["f_ghz"])):
        I = r["I"][m] if "I" in r else None
        rec = {"m": m + 1, "f_mhz": float(r["f_ghz"][m] * 1e3), "Q": float(r["Q"][m])}
        if I is not None:
            if leg_idx:
                cur = _align(I[leg_idx])
                rec["pattern_source"] = "rung capacitor currents"
            else:  # high-pass: ring-segment (mesh) currents carry the same harmonic
                cur = _align(I[top_idx])
                rec["pattern_source"] = "top end-ring capacitor currents"
            p = _dft_power(cur)
            k = int(np.argmax(p))
            rec["cap_current_sq"] = float(np.sum(np.abs(I) ** 2))
            rec.update(k=k, k_purity=float(p[k]), dft_power=p.tolist(),
                       currents_re=np.real(cur).tolist(), currents_im=np.imag(cur).tolist())
            if top_idx and bot_idx:
                t, b = _align(I[top_idx]), I[bot_idx]
                b = b * np.exp(-1j * np.angle(I[top_idx][np.argmax(np.abs(I[top_idx]))]))
                ratio = np.vdot(t, b) / max(np.vdot(t, t).real, 1e-300)
                rec["ring_bot_over_top"] = float(np.real(ratio))
                rec["ring_current_frac"] = float(np.sum(np.abs(I[top_idx + bot_idx]) ** 2) / np.sum(np.abs(I) ** 2))
        if "B" in r:
            B = r["B"][m]
            roi = B[sl["roi"]]
            bt = np.sqrt(np.abs(roi[:, 0]) ** 2 + np.abs(roi[:, 1]) ** 2)
            rec["roi_bt_cv"] = float(bt.std() / bt.mean())
            rec["roi_transverse_frac"] = float((bt ** 2).sum() / (np.abs(roi) ** 2).sum())
            c = B[sl["axial"]][np.argmin(np.sum(np.array(pr["points"][:n_ax])[:, :2] ** 2, 1))]
            rec["center_bt_over_roi_mean"] = float(np.sqrt(abs(c[0]) ** 2 + abs(c[1]) ** 2) / bt.mean())
        modes.append(rec)

    # Modes whose capacitors carry (almost) no current are not ladder modes: end-ring / rung-half resonances set by
    # stray capacitance. Flag them so they are not mistaken for a k mode.
    isq = [md.get("cap_current_sq", 0) for md in modes]
    ref = max(isq) if isq else 0
    for md in modes:
        md["cap_free"] = bool(ref and md.get("cap_current_sq", ref) < 1e-2 * ref)
        if md["cap_free"]:
            md["k_dft"] = md.pop("k", None)
            md["k_purity"] = 0.0
    # Label modes. k=1 candidates: dominant k=1 with transverse field.
    k1 = [md for md in modes if md.get("k") == 1 and md.get("roi_transverse_frac", 1) > 0.9]
    k1 = sorted(k1, key=lambda md: md.get("roi_bt_cv", 0))[:2]
    for md in modes:
        kk = md.get("k")
        if md.get("cap_free"):
            md["label"] = "no capacitor current (ring/stray-C)"
        elif kk is None:
            md["label"] = "?"
        elif kk == 0:
            if des["coil_type"] == "lowpass":
                md["label"] = "k=0 (rungs in phase)"
            else:
                rt = md.get("ring_bot_over_top", 0)
                md["label"] = "end-ring k=0 " + ("co-rotating" if rt < 0 else "counter-rotating")
        else:
            md["label"] = f"k={kk}"
        md["homogeneous"] = md in k1
    out = {"design": des, "modes": modes, "N": N}
    if len(k1) == 2:
        fa, fb = k1[0]["f_mhz"], k1[1]["f_mhz"]
        out["homogeneous_pair"] = {"m": [k1[0]["m"], k1[1]["m"]], "f_mhz": [fa, fb], "f_mean_mhz": (fa + fb) / 2,
                                   "split_mhz": abs(fa - fb), "split_rel": abs(fa - fb) / ((fa + fb) / 2),
                                   "roi_bt_cv": [k1[0]["roi_bt_cv"], k1[1]["roi_bt_cv"]]}
        if "B" in r:
            q = quadrature(r["B"][k1[0]["m"] - 1], r["B"][k1[1]["m"] - 1])
            roi = q["b1p"][sl["roi"]]
            out["homogeneous_pair"]["b1plus_roi_cv"] = float(np.abs(roi).std() / np.abs(roi).mean())
            out["homogeneous_pair"]["b1minus_over_b1plus_roi"] = float(np.abs(q["b1m"][sl["roi"]]).mean() / np.abs(roi).mean())
    elif len(k1) == 1:
        out["homogeneous_pair"] = {"m": [k1[0]["m"]], "f_mhz": [k1[0]["f_mhz"]], "f_mean_mhz": k1[0]["f_mhz"],
                                   "roi_bt_cv": [k1[0]["roi_bt_cv"]]}
    return out


def quadrature(Ba: np.ndarray, Bb: np.ndarray) -> dict:
    """B1+ / B1- of the k=1 pair driven 90 degrees apart, Ba + i Bb or Ba - i Bb, whichever gives the larger B1+.

    Each eigenvector is first made real (it is, up to a phase, for a lossless mode) and the pair is orthogonalised
    and equal-amplitude normalised over all probes so the combination is a pure rotating field.
    """
    def realify(B):
        v = B.reshape(-1)
        j = np.argmax(np.abs(v))
        return np.real(B * np.exp(-1j * np.angle(v[j])))
    a, b = realify(Ba), realify(Bb)
    a = a / np.linalg.norm(a)
    b = b - np.sum(a * b) * a
    b = b / np.linalg.norm(b)
    best = None
    for s in (1, -1):
        B = a + s * 1j * b
        bp = (B[:, 0] + 1j * B[:, 1]) / 2
        bm = np.conj(B[:, 0] - 1j * B[:, 1]) / 2
        if best is None or np.abs(bp).sum() > np.abs(best["b1p"]).sum():
            best = {"b1p": bp, "b1m": bm, "sign": s}
    return best


def compare(analytic: dict, res: dict) -> list[dict]:
    """Pair Palace modes with ladder-model modes by k (degenerate pairs averaged)."""
    rows = []
    spec = analytic["spectrum"]
    for k, fa in zip(spec["k"], spec["f_hz"]):
        ms = [m for m in res["modes"] if m.get("k") == k and m.get("k_purity", 0) > 0.5]
        if k == 0 and res["design"]["coil_type"] == "lowpass":
            continue
        if res["design"]["coil_type"] != "lowpass" and k == 0:
            ms = [m for m in ms if "counter" in m.get("label", "")]
        ms = sorted(ms, key=lambda m: m["f_mhz"])[: (2 if 0 < k < res["N"] // 2 else 1)]
        fp = [m["f_mhz"] for m in ms]
        rows.append({"k": k, "f_ladder_mhz": fa / 1e6, "f_palace_mhz": fp,
                     "ratio": (float(np.mean(fp)) / (fa / 1e6)) if fp and fa > 0 else None})
    return rows
