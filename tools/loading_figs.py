"""Figures for the phantom loading study (results/loading_sweep.json + runs/loading/maps.json).

    python tools/loading_figs.py --out <dir>        # cards, sweep plots, maps, GIF/MP4, index.html
    python tools/loading_figs.py --out docs/img --docs   # the page / README subset

Simulation only. Every number on a figure comes from results/loading_sweep.json.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402
from matplotlib.patches import Circle, Rectangle  # noqa: E402
from PIL import Image  # noqa: E402

INK, INK2 = "#0b0b0b", "#52514e"
B_CMAP, E_CMAP = "magma", "viridis"
HONEST = "Palace FEM · homogeneous cylindrical phantom · simulation only"
DPI = 100
R_COIL = 0.14
FFMPEG = shutil.which("ffmpeg") or r"C:\Users\ben\AppData\Local\Microsoft\WinGet\Links\ffmpeg.exe"

plt.rcParams.update({"font.size": 15, "axes.titlesize": 16, "axes.labelsize": 15, "font.family": "DejaVu Sans",
                     "axes.spines.top": False, "axes.spines.right": False})


def load():
    sw = json.loads((ROOT / "results" / "loading_sweep.json").read_text(encoding="utf-8"))
    maps = json.loads((ROOT / "runs" / "loading" / "maps.json").read_text())
    by = {c["name"]: c for c in sw["cases"]}
    return sw, by, maps


def grid_img(xy, v, h=0.004):
    """Probe values on a regular grid -> 2D image (NaN outside), extent in cm."""
    xy = np.asarray(xy)
    v = np.asarray(v, float)
    x0, y0 = xy[:, 0].min(), xy[:, 1].min()
    ix = np.rint((xy[:, 0] - x0) / h).astype(int)
    iy = np.rint((xy[:, 1] - y0) / h).astype(int)
    img = np.full((iy.max() + 1, ix.max() + 1) + v.shape[1:], np.nan)
    img[iy, ix] = v
    ext = [100 * (x0 - h / 2), 100 * (x0 + h * ix.max() + h / 2), 100 * (y0 - h / 2), 100 * (y0 + h * iy.max() + h / 2)]
    return img, ext


def name(r=0.08, eps=78, sig=0.6):
    return f"r{round(r * 1000):03d}_e{eps:g}_s{sig:g}".replace(".", "p")


def _bare(ax):
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def axial(ax, case_maps, key="b1p_axial", norm=None, cmap=B_CMAP, r_ph=0.08, scale=1e6):
    img, ext = grid_img(case_maps["axial_xy"], np.asarray(case_maps[key]) * scale)
    im = ax.imshow(img, origin="lower", extent=ext, cmap=cmap, norm=norm, interpolation="bilinear")
    if r_ph:
        ax.add_patch(Circle((0, 0), 100 * r_ph, fill=False, ec="white", lw=1.2, ls="--", alpha=0.9))
    phi = np.pi / 16 + 2 * np.pi * np.arange(16) / 16
    ax.plot(100 * R_COIL * np.cos(phi), 100 * R_COIL * np.sin(phi), "s", ms=4, color=INK2)
    ax.set_xlim(-15, 15)
    ax.set_ylim(-15, 15)
    _bare(ax)
    return im


def coronal(ax, case_maps, key="b1p_coronal", norm=None, cmap=B_CMAP, r_ph=0.08, scale=1e6):
    img, ext = grid_img(case_maps["coronal_xz"], np.asarray(case_maps[key]) * scale)
    im = ax.imshow(img.T, origin="lower", extent=[ext[2], ext[3], ext[0], ext[1]], cmap=cmap, norm=norm,
                   interpolation="bilinear")
    if r_ph:
        ax.add_patch(Rectangle((-10, -100 * r_ph), 20, 200 * r_ph, fill=False, ec="white", lw=1.2, ls="--"))
    _bare(ax)
    return im


def foot(fig, extra=""):
    fig.text(0.5, 0.018, HONEST + (" · " + extra if extra else ""), ha="center", va="bottom", fontsize=12, color=INK2)


def save(fig, out, name_):
    fig.savefig(out / name_, dpi=DPI, facecolor="white")
    plt.close(fig)
    print("wrote", out / name_)


def coil_txt():
    return "16-rung low-pass birdcage, 28 cm, 3 T · C fixed at the empty-coil tune · phantom r 8 cm, length 20 cm"


# ------------------------------------------------------------------ cards (1080 x 1350)
def card_b1(by, maps, out):
    e, l, eo, so = name(0.08, 1, 0), name(0.08, 78, 0.6), name(0.08, 78, 0), name(0.08, 1, 0.6)
    vmax = max(np.nanmax(np.asarray(maps[k]["b1p_axial"])) for k in (e, l, eo, so) if k in maps) * 1e6
    norm = Normalize(0, vmax)
    fig = plt.figure(figsize=(10.8, 13.5), dpi=DPI)
    fig.text(0.5, 0.955, "Quadrature |B1+| at the centre slice", ha="center", fontsize=30, weight="bold", color=INK)
    fig.text(0.5, 0.925, coil_txt(), ha="center", fontsize=12.5, color=INK2)
    for i, (k, t) in enumerate(((e, "Empty coil"), (l, "Loaded: εr 78, σ 0.6 S/m"))):
        ax = fig.add_axes([0.04 + 0.48 * i, 0.47, 0.44, 0.40])
        im = axial(ax, maps[k], norm=norm, r_ph=0.08 if k == l else 0)
        ax.set_title(t, fontsize=20, color=INK, pad=10)
        c = by[k]
        fig.text(0.26 + 0.48 * i, 0.44, f"{c['b1plus_roi_cv'] * 100:.1f} %", ha="center", fontsize=42, weight="bold", color=INK)
        fig.text(0.26 + 0.48 * i, 0.415, "ROI std/mean of |B1+|", ha="center", fontsize=14, color=INK2)
        fig.text(0.26 + 0.48 * i, 0.39, f"centre / ROI mean {c['b1plus_centre_over_mean']:.2f}", ha="center", fontsize=14, color=INK2)
    cax = fig.add_axes([0.25, 0.355, 0.5, 0.014])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_label("|B1+| (µT) at 1 µT mean over the ROI · one scale for all four maps", fontsize=13)
    if eo in maps and so in maps:
        fig.text(0.5, 0.285, "Each property alone", ha="center", fontsize=20, weight="bold", color=INK)
        for i, (k, t) in enumerate(((eo, "εr 78 only (σ 0)"), (so, "σ 0.6 S/m only (εr 1)"))):
            ax = fig.add_axes([0.1 + 0.44 * i, 0.075, 0.36, 0.18])
            axial(ax, maps[k], norm=norm)
            c = by[k]
            ax.set_title(f"{t}: {c['b1plus_roi_cv'] * 100:.1f} %", fontsize=15, color=INK, pad=6)
    fig.text(0.5, 0.05, "ROI: r < 8.4 cm, |z| < 4.2 cm · rungs as grey squares · dashed: phantom", ha="center", fontsize=12, color=INK2)
    foot(fig)
    save(fig, out, "card_b1_unloaded_vs_loaded.png")


def lorentz(f, f0, q):
    return 1 / np.sqrt(1 + (2 * q * (f / f0 - 1)) ** 2)


def card_resonance(by, out):
    u, l = by[name(0.08, 1, 0)], by[name(0.08, 78, 0.6)]
    fu, fl = u["f_full_mhz"], l["f_full_mhz"]
    qu, ql = u["Q_L"], l["Q_L"]
    f = np.linspace(124.5, 130.5, 3000)
    fig = plt.figure(figsize=(10.8, 13.5), dpi=DPI)
    fig.text(0.5, 0.955, "Resonance of the k = 1 mode", ha="center", fontsize=30, weight="bold", color=INK)
    fig.text(0.5, 0.925, coil_txt(), ha="center", fontsize=12.5, color=INK2)
    fig.text(0.5, 0.902, "copper strips and shield + 0.68 Ω ESR per capacitor in both", ha="center", fontsize=12.5, color=INK2)
    for i, (fc, q, t, lab, col) in enumerate(((fu, qu, "Empty coil", "Q_U", "C0"), (fl, ql, "Loaded: εr 78, σ 0.6 S/m", "Q_L", "C1"))):
        ax = fig.add_axes([0.08 + 0.48 * i, 0.42, 0.40, 0.42])
        ax.plot(f, lorentz(f, fc, q), color=col, lw=2.5)
        ax.set_ylim(0, 1.08)
        ax.set_xlim(f[0], f[-1])
        ax.set_xticks([125, 127, 129])
        ax.set_xlabel("frequency (MHz)")
        if i == 0:
            ax.set_ylabel("|response| (normalised)")
        else:
            ax.set_yticklabels([])
        ax.set_title(t, fontsize=20, color=INK, pad=10)
        bw = fc / q
        ax.annotate("", xy=(fc - bw / 2, 0.707), xytext=(fc + bw / 2, 0.707),
                    arrowprops=dict(arrowstyle="<->", color=INK2, lw=1.2, shrinkA=0, shrinkB=0))
        ax.text(fc, 0.64, f"−3 dB width {bw * 1000:.0f} kHz", ha="center", fontsize=13, color=INK2,
                bbox=dict(fc="white", ec="none", pad=1))
        fig.text(0.28 + 0.48 * i, 0.33, f"{lab} = {q:.0f}", ha="center", fontsize=46, weight="bold", color=INK)
        fig.text(0.28 + 0.48 * i, 0.30, f"f = {fc:.2f} MHz", ha="center", fontsize=16, color=INK2)
    fig.text(0.5, 0.20, f"Q_U / Q_L = {qu / ql:.1f}", ha="center", fontsize=34, weight="bold", color=INK)
    fig.text(0.5, 0.165, f"the phantom takes {100 * (1 - ql / qu):.0f} % of the power fed to the coil", ha="center", fontsize=16, color=INK2)
    fig.text(0.5, 0.075, "Curves: single-pole response from Palace's complex eigenfrequency, Q = Re f / (2 Im f).\n"
             "Copper loss on the zero-thickness strips is an upper bound (one face carries the current), so Q_U is on the low side.",
             ha="center", fontsize=11.5, color=INK2)
    foot(fig)
    save(fig, out, "card_resonance_q.png")


# ------------------------------------------------------------------ sweep plots
def sweeps(by):
    sig = sorted([c for c in by.values() if c["r"] == 0.08 and c["eps"] == 78], key=lambda c: c["sigma"])
    eps = sorted([c for c in by.values() if c["r"] == 0.08 and c["sigma"] == 0], key=lambda c: c["eps"])
    rad = sorted([c for c in by.values() if c["eps"] == 78 and c["sigma"] == 0.6], key=lambda c: c["r"])
    return sig, eps, rad


def plot_q_fshift(by, out, fname="loading_q_vs_sigma_fshift_vs_eps.png", size=(10.8, 10.8)):
    sig, eps, _ = sweeps(by)
    u = by[name(0.08, 1, 0)]
    fig, (a, b) = plt.subplots(2, 1, figsize=size, dpi=DPI)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.86, bottom=0.12, hspace=0.42)
    fig.text(0.5, 0.95, "Phantom σ and εr: what moves Q and what moves f", ha="center", fontsize=24, weight="bold")
    fig.text(0.5, 0.915, coil_txt(), ha="center", fontsize=11.5, color=INK2)
    sq = [c for c in sig if "Q_L" in c]
    a.plot([c["sigma"] for c in sq], [c["Q_L"] for c in sq], "o-", color="C0", lw=2, label="Q_L (phantom + copper + ESR)")
    a.plot([c["sigma"] for c in sig if c["sigma"] > 0], [c["Q_phantom"][0] for c in sig if c["sigma"] > 0], "s--", color="C1",
           lw=1.5, ms=5, label="phantom loss only")
    a.axhline(u["Q_L"], color="C2", ls=":", lw=2, label=f"Q_U, empty coil = {u['Q_L']:.0f}")
    e_q = [c["Q_L"] for c in eps if "Q_L" in c]
    a.set_yscale("log")
    a.set_xlabel("phantom conductivity σ (S/m), εr 78")
    a.set_ylabel("Q of k = 1")
    a.set_ylim(10, 3000)
    a.legend(fontsize=12, frameon=False, loc="upper right")
    a.set_title(f"σ 0 → {sq[-1]['sigma']:g} S/m: Q_L {sq[0]['Q_L']:.0f} → {sq[-1]['Q_L']:.0f};   "
                f"εr 1 → 78 at σ 0: Q_L {min(e_q):.0f}–{max(e_q):.0f}", fontsize=14, loc="left")
    b.plot([c["eps"] for c in eps], [100 * c["f_shift_rel"] for c in eps], "o-", color="C0", lw=2, label="εr sweep (σ 0)")
    b.plot([78 for _ in sig], [100 * c["f_shift_rel"] for c in sig], "_", color="C1", ms=26, mew=2,
           label=f"σ 0 to {sig[-1]['sigma']:g} S/m at εr 78")
    b.axhline(0, color=INK2, lw=0.8)
    b.set_xlabel("phantom permittivity εr")
    b.set_ylabel("k = 1 shift from empty coil (%)")
    b.legend(fontsize=12, frameon=False, loc="lower left")
    fs = [100 * c["f_shift_rel"] for c in sig]
    b.set_title(f"εr 1 → {eps[-1]['eps']:g} (σ 0): {100 * eps[-1]['f_shift_rel']:+.2f} %;   σ 0 → {sig[-1]['sigma']:g} S/m at εr 78: "
                f"{fs[0]:+.2f} → {fs[-1]:+.2f} %", fontsize=14, loc="left")
    b.annotate("conductivity pulls k = 1\nback towards the empty coil", xy=(78, fs[-1]), xytext=(42, fs[-1] - 0.15),
               fontsize=12, color="C1", arrowprops=dict(arrowstyle="->", color="C1"))
    foot(fig)
    save(fig, out, fname)


def plot_supporting(by, out):
    """2 x 2: Q_U/Q_L, shift and homogeneity against sigma and eps_r; and the radius sweep if present."""
    sig, eps, rad = sweeps(by)
    fig, ax = plt.subplots(2, 2, figsize=(10.8, 10.8), dpi=DPI)
    fig.subplots_adjust(left=0.09, right=0.91, top=0.9, bottom=0.1, hspace=0.38, wspace=0.62)
    fig.text(0.5, 0.95, "Loading sweeps: Q, frequency and B1+ homogeneity", ha="center", fontsize=20, weight="bold")
    for col, (xs, xl) in enumerate(((sig, "σ (S/m), εr 78"), (eps, "εr, σ 0"))):
        x = [c["sigma"] if col == 0 else c["eps"] for c in xs]
        a = ax[0, col]
        a.plot(x, [c.get("QU_over_QL", np.nan) for c in xs], "o-", color="C0", lw=2)
        a.set_ylabel("Q_U / Q_L", color="C0")
        a.set_xlabel(xl)
        a2 = a.twinx()
        a2.plot(x, [100 * c["f_shift_rel"] for c in xs], "s--", color="C1", lw=1.5, ms=5)
        a2.set_ylabel("k = 1 shift (%)", color="C1")
        a2.spines["right"].set_visible(True)
        b = ax[1, col]
        b.plot(x, [100 * c["b1plus_roi_cv"] for c in xs], "o-", color="C0", lw=2, label="ROI std/mean")
        b.plot(x, [100 * c["b1plus_phantom_slice_cv"] for c in xs], "s--", color="C1", lw=1.5, ms=5, label="z = 0 slice in the phantom")
        b.set_ylabel("|B1+| std/mean (%)")
        b.set_xlabel(xl)
        b.legend(fontsize=11, frameon=False)
    foot(fig)
    save(fig, out, "loading_sweeps_2x2.png")
    if len(rad) >= 3:
        fig, a = plt.subplots(figsize=(10.8, 6.5), dpi=DPI)
        fig.subplots_adjust(left=0.1, right=0.88, top=0.86, bottom=0.17)
        r = [100 * c["r"] for c in rad]
        a.plot(r, [c.get("QU_over_QL", np.nan) for c in rad], "o-", color="C0", lw=2)
        a.set_xlabel("phantom radius (cm), εr 78, σ 0.6 S/m, length 20 cm")
        a.set_ylabel("Q_U / Q_L", color="C0")
        a2 = a.twinx()
        a2.plot(r, [100 * c["b1plus_phantom_slice_cv"] for c in rad], "s--", color="C1", lw=1.5)
        a2.set_ylabel("|B1+| std/mean, z = 0 slice in the phantom (%)", color="C1")
        a2.spines["right"].set_visible(True)
        fig.text(0.5, 0.93, "Phantom radius", ha="center", fontsize=20, weight="bold")
        foot(fig)
        save(fig, out, "loading_radius.png")


def strip_b1(by, maps, out, fname="loading_b1_strip.png"):
    sig_pick = [0, 0.3, 0.6, 1.0, 2.0]
    eps_pick = [1, 10, 30, 50, 78]
    rows = [[name(0.08, 78, s) for s in sig_pick], [name(0.08, e, 0) for e in eps_pick]]
    allk = [k for r in rows for k in r if k in maps]
    vmax = max(np.nanmax(np.asarray(maps[k]["b1p_axial"])) for k in allk) * 1e6
    norm = Normalize(0, vmax)
    fig = plt.figure(figsize=(10.8, 5.6), dpi=DPI)
    fig.text(0.5, 0.94, "Quadrature |B1+| at the centre slice, one colour scale", ha="center", fontsize=18, weight="bold")
    im = None
    for r, (keys, lab, vals) in enumerate(((rows[0], "σ sweep, εr 78", sig_pick), (rows[1], "εr sweep, σ 0", eps_pick))):
        fig.text(0.015, 0.66 - 0.4 * r, lab, rotation=90, va="center", fontsize=13, color=INK)
        for i, (k, v) in enumerate(zip(keys, vals)):
            if k not in maps:
                continue
            ax = fig.add_axes([0.04 + 0.172 * i, 0.47 - 0.4 * r, 0.165, 0.36])
            im = axial(ax, maps[k], norm=norm, r_ph=0.08 if (r == 0 or v != 1) else 0)
            t = f"σ {v:g} S/m" if r == 0 else f"εr {v:g}"
            ax.set_title(f"{t} · {by[k]['b1plus_roi_cv'] * 100:.0f} %", fontsize=12, pad=3)
    cax = fig.add_axes([0.905, 0.12, 0.015, 0.7])
    fig.colorbar(im, cax=cax).set_label("|B1+| (µT), 1 µT ROI mean", fontsize=11)
    fig.text(0.5, 0.005, HONEST + " · % = ROI std/mean", ha="center", va="bottom", fontsize=10.5, color=INK2)
    save(fig, out, fname)


def e_maps(by, maps, out, key=None):
    key = key or name(0.08, 78, 0.6)
    c, m = by[key], maps[key]
    rp, lp = c["r"], 0.10
    xy, xz = np.asarray(m["axial_xy"]), np.asarray(m["coronal_xz"])
    ina = np.hypot(xy[:, 0], xy[:, 1]) < rp - 0.003   # same 3 mm margin as the peak values
    inc = (np.abs(xz[:, 0]) < rp - 0.003) & (np.abs(xz[:, 1]) < lp - 0.003)
    Ea, Ec = np.asarray(m["E_axial"]), np.asarray(m["E_coronal"])
    mm = {"axial_xy": m["axial_xy"], "coronal_xz": m["coronal_xz"],
          "Ea": np.where(ina, Ea, np.nan), "Ec": np.where(inc, Ec, np.nan),
          "Sa": np.where(ina, 0.5 * c["sigma"] * Ea ** 2, np.nan), "Sc": np.where(inc, 0.5 * c["sigma"] * Ec ** 2, np.nan)}
    en = Normalize(0, np.nanmax([np.nanmax(mm["Ea"]), np.nanmax(mm["Ec"])]))
    sn = Normalize(0, np.nanmax([np.nanmax(mm["Sa"]), np.nanmax(mm["Sc"])]))
    fig = plt.figure(figsize=(10.8, 10.8), dpi=DPI)
    fig.text(0.5, 0.955, "Electric field and σ|E|²/2 inside the saline phantom", ha="center", fontsize=24, weight="bold")
    fig.text(0.5, 0.925, f"εr {c['eps']:g}, σ {c['sigma']:g} S/m, r {100 * rp:g} cm, length 20 cm · quadrature drive, 1 µT mean |B1+| over the ROI",
             ha="center", fontsize=12.5, color=INK2)
    rows = (("Ea", "Ec", en, E_CMAP, "|E| (V/m)", f"peak |E| {c['E_peak_phantom_vpm']:.0f} V/m"),
            ("Sa", "Sc", sn, "inferno", "σ|E|²/2 (W/m³), no tissue density", f"peak σ|E|²/2 {c['sar_like_peak_wpm3']:.0f} W/m³"))
    for r, (ka, kc, nrm, cm, lab, pk) in enumerate(rows):
        y0 = 0.5 - 0.43 * r
        a1 = fig.add_axes([0.05, y0, 0.33, 0.36])
        im = axial(a1, mm, ka, nrm, cm, r_ph=0, scale=1)
        a1.set_xlim(-9, 9)
        a1.set_ylim(-9, 9)
        a1.set_title("z = 0 (axial)", fontsize=15)
        a2 = fig.add_axes([0.42, y0, 0.42, 0.36])
        coronal(a2, mm, kc, nrm, cm, r_ph=0, scale=1)
        a2.set_xlim(-11, 11)
        a2.set_ylim(-9, 9)
        a2.set_title("y = 0 (coronal, z across)", fontsize=15)
        cax = fig.add_axes([0.87, y0 + 0.03, 0.018, 0.30])
        fig.colorbar(im, cax=cax).set_label(lab, fontsize=12)
        fig.text(0.5, y0 - 0.025, pk + " (probes > 3 mm inside the phantom)", ha="center", fontsize=14, weight="bold", color=INK)
    foot(fig, "z = 0 is a symmetry plane, so E there points along z")
    save(fig, out, "loading_e_sar_saline.png")


# ------------------------------------------------------------------ animations
def frames_to_media(frames, out, base, fps, gif_colors=128):
    tmp = out / f"_{base}_frames"
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir()
    for i, f in enumerate(frames):
        f.save(tmp / f"f{i:04d}.png")
    gif = out / f"{base}.gif"
    pal = tmp / "pal.png"
    subprocess.run([FFMPEG, "-y", "-loglevel", "error", "-framerate", str(fps), "-i", str(tmp / "f%04d.png"),
                    "-vf", f"palettegen=max_colors={gif_colors}:stats_mode=full", str(pal)], check=True)
    subprocess.run([FFMPEG, "-y", "-loglevel", "error", "-framerate", str(fps), "-i", str(tmp / "f%04d.png"), "-i", str(pal),
                    "-lavfi", "paletteuse=dither=bayer:bayer_scale=4", "-loop", "0", str(gif)], check=True)
    mp4 = out / f"{base}.mp4"
    n_loop = max(1, int(round(8 * fps / len(frames))))  # MP4: about 8 s of loops
    subprocess.run([FFMPEG, "-y", "-loglevel", "error", "-stream_loop", str(n_loop - 1), "-framerate", str(fps), "-i",
                    str(tmp / "f%04d.png"), "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18", "-movflags", "+faststart",
                    str(mp4)], check=True)
    frames[0].save(out / f"{base}_still.png")
    shutil.rmtree(tmp)
    print("wrote", gif, f"{gif.stat().st_size / 1e6:.2f} MB", mp4)


def fig_to_img(fig):
    fig.canvas.draw()
    im = Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[..., :3].copy())
    plt.close(fig)
    return im


def anim_rotating(by, maps, out, nfr=36, fps=12):
    """Card: instantaneous transverse B over one RF cycle, empty vs loaded, same colour scale."""
    keys = (name(0.08, 1, 0), name(0.08, 78, 0.6))
    data = []
    for k in keys:
        m = maps[k]
        Bq = np.asarray(m["Bq_axial_re"]) + 1j * np.asarray(m["Bq_axial_im"])
        data.append((m, Bq * 1e6))
    vmax = max(np.nanmax(np.hypot(np.abs(Bq[:, 0]), np.abs(Bq[:, 1]))) for _, Bq in data)
    norm = Normalize(0, vmax)
    frames = []
    for t in range(nfr):
        ph = np.exp(1j * 2 * np.pi * t / nfr)
        fig = plt.figure(figsize=(10.8, 13.5), dpi=DPI)
        fig.text(0.5, 0.955, "Rotating B1 over one RF cycle", ha="center", fontsize=30, weight="bold", color=INK)
        fig.text(0.5, 0.925, coil_txt(), ha="center", fontsize=12.5, color=INK2)
        for i, ((m, Bq), k, t_) in enumerate(zip(data, keys, ("Empty coil", "Loaded: εr 78, σ 0.6 S/m"))):
            inst = np.real(Bq * ph)
            mag = np.hypot(inst[:, 0], inst[:, 1])
            ax = fig.add_axes([0.04 + 0.48 * i, 0.42, 0.44, 0.44])
            mm = {"axial_xy": m["axial_xy"], "inst": mag}
            im = axial(ax, mm, "inst", norm, B_CMAP, r_ph=0.08 if i else 0, scale=1)
            xy = np.asarray(m["axial_xy"])
            sel = (np.abs(np.rint(xy[:, 0] / 0.004)) % 5 == 0) & (np.abs(np.rint(xy[:, 1] / 0.004)) % 5 == 0)
            ax.quiver(100 * xy[sel, 0], 100 * xy[sel, 1], inst[sel, 0], inst[sel, 1], color="white", scale=vmax * 14,
                      width=0.005, alpha=0.9)
            ax.set_title(t_, fontsize=20, color=INK, pad=10)
            c = by[k]
            fig.text(0.26 + 0.48 * i, 0.37, f"{c['b1plus_roi_cv'] * 100:.1f} %", ha="center", fontsize=40, weight="bold", color=INK)
            fig.text(0.26 + 0.48 * i, 0.345, "ROI std/mean of |B1+|", ha="center", fontsize=14, color=INK2)
            fig.text(0.26 + 0.48 * i, 0.315, f"|B1−| / |B1+| = {c['b1minus_over_b1plus_roi'] * 100:.0f} %", ha="center", fontsize=14, color=INK2)
            fig.text(0.26 + 0.48 * i, 0.285, f"f = {c['f_mean_mhz']:.2f} MHz", ha="center", fontsize=14, color=INK2)
        cax = fig.add_axes([0.25, 0.22, 0.5, 0.014])
        fig.colorbar(im, cax=cax, orientation="horizontal").set_label(
            "instantaneous |B_t| (µT), 1 µT mean |B1+| · same scale", fontsize=13)
        fig.text(0.5, 0.10, "Ideal quadrature drive of the two k = 1 modes (each from its own symmetry run).\n"
                 "Empty coil: the field turns as one block. Loaded: it turns unevenly and the\n"
                 "counter-rotating part |B1−| grows, so |B1+| is no longer flat.", ha="center", fontsize=13, color=INK2)
        foot(fig)
        frames.append(fig_to_img(fig))
    frames_to_media(frames, out, "card_fields_unloaded_vs_loaded", fps)


def anim_sweep(by, maps, out, which, fps=10):
    sig, eps, _ = sweeps(by)
    xs = [c for c in (sig if which == "sigma" else eps) if c["name"] in maps]
    vmax = max(np.nanmax(np.asarray(maps[c["name"]]["b1p_axial"])) for c in (sig + eps) if c["name"] in maps) * 1e6
    norm = Normalize(0, vmax)
    xv = [c["sigma"] if which == "sigma" else c["eps"] for c in xs]
    if which == "sigma":
        yv, ylab, title = [c.get("Q_L", np.nan) for c in xs], "Q_L", "σ sweep at εr 78"
    else:
        yv, ylab, title = [100 * c["f_shift_rel"] for c in xs], "k = 1 shift (%)", "εr sweep at σ 0"
    imgs_a = [grid_img(maps[c["name"]]["axial_xy"], np.asarray(maps[c["name"]]["b1p_axial"]) * 1e6)[0] for c in xs]
    imgs_c = [grid_img(maps[c["name"]]["coronal_xz"], np.asarray(maps[c["name"]]["b1p_coronal"]) * 1e6)[0] for c in xs]
    m0 = maps[xs[0]["name"]]
    ext_a = grid_img(m0["axial_xy"], np.zeros(len(m0["axial_xy"])))[1]
    ext_c = grid_img(m0["coronal_xz"], np.zeros(len(m0["coronal_xz"])))[1]
    seq = []  # (index a, index b, blend weight)
    for i in range(len(xs)):
        seq += [(i, i, 0.0), (i, i, 0.0), (i, (i + 1) % len(xs), 0.5)]
    frames = []
    for ia, ib, w in seq:
        cur = xs[ia]
        fig = plt.figure(figsize=(10.8, 10.8), dpi=DPI)
        fig.text(0.5, 0.955, f"|B1+| through the {title}", ha="center", fontsize=26, weight="bold")
        fig.text(0.5, 0.925, coil_txt(), ha="center", fontsize=11.5, color=INK2)
        a = fig.add_axes([0.03, 0.43, 0.42, 0.45])
        im = a.imshow((1 - w) * imgs_a[ia] + w * imgs_a[ib], origin="lower", extent=ext_a, cmap=B_CMAP, norm=norm,
                      interpolation="bilinear")
        if not (which == "eps" and cur["eps"] == 1):
            a.add_patch(Circle((0, 0), 8, fill=False, ec="white", lw=1.2, ls="--"))
        a.set_xlim(-15, 15)
        a.set_ylim(-15, 15)
        a.axis("off")
        a.set_title("z = 0 (axial)", fontsize=15)
        c = fig.add_axes([0.5, 0.43, 0.47, 0.45])
        ic = (1 - w) * imgs_c[ia] + w * imgs_c[ib]
        c.imshow(ic.T, origin="lower", extent=[ext_c[2], ext_c[3], ext_c[0], ext_c[1]], cmap=B_CMAP, norm=norm,
                 interpolation="bilinear")
        if not (which == "eps" and cur["eps"] == 1):
            c.add_patch(Rectangle((-10, -8), 20, 16, fill=False, ec="white", lw=1.2, ls="--"))
        c.axis("off")
        c.set_title("y = 0 (coronal, z across)", fontsize=15)
        cax = fig.add_axes([0.25, 0.405, 0.5, 0.012])
        fig.colorbar(im, cax=cax, orientation="horizontal").set_label(
            "|B1+| (µT), 1 µT ROI mean · one scale for both sweeps", fontsize=11)
        p = fig.add_axes([0.13, 0.09, 0.8, 0.19])
        p.plot(xv, yv, "o-", color="C0", lw=2)
        p.plot([xv[ia]], [yv[ia]], "o", ms=16, mfc="none", mec="C3", mew=3)
        p.set_xlabel("σ (S/m)" if which == "sigma" else "εr")
        p.set_ylabel(ylab)
        if which == "sigma":
            p.set_yscale("log")
        val = f"σ = {cur['sigma']:g} S/m" if which == "sigma" else f"εr = {cur['eps']:g}"
        fig.text(0.5, 0.325, f"{val}   ·   ROI std/mean {cur['b1plus_roi_cv'] * 100:.1f} %   ·   "
                 + (f"Q_L {cur.get('Q_L', float('nan')):.0f}" if which == "sigma" else f"shift {100 * cur['f_shift_rel']:+.2f} %"),
                 ha="center", fontsize=17, weight="bold", color=INK)
        foot(fig, "frames between computed points are blends")
        frames.append(fig_to_img(fig))
    frames_to_media(frames, out, f"sweep_{which}_b1", fps)


def index_html(out, by, sw):
    u, l = by[name(0.08, 1, 0)], by[name(0.08, 78, 0.6)]
    sig, eps, rad = sweeps(by)
    nan = float("nan")
    rows = "".join(
        f"<tr><td>{c['sigma']:g}</td><td>{c['f_mean_mhz']:.3f}</td><td>{100 * c['f_shift_rel']:+.2f}</td>"
        f"<td>{(c['Q_phantom'][0] if c['sigma'] > 0 else nan):.1f}</td>"
        f"<td>{c.get('Q_L', nan):.1f}</td><td>{c.get('QU_over_QL', nan):.2f}</td><td>{100 * c['b1plus_roi_cv']:.1f}</td>"
        f"<td>{c['E_peak_phantom_vpm']:.0f}</td><td>{c['sar_like_peak_wpm3']:.0f}</td></tr>" for c in sig)
    rows2 = "".join(
        f"<tr><td>{c['eps']:g}</td><td>{c['f_mean_mhz']:.3f}</td><td>{100 * c['f_shift_rel']:+.2f}</td><td>{c.get('Q_L', nan):.0f}</td>"
        f"<td>{c.get('QU_over_QL', nan):.3f}</td><td>{100 * c['b1plus_roi_cv']:.1f}</td><td>{100 * c['b1plus_phantom_slice_cv']:.1f}</td>"
        f"<td>{c['b1plus_centre_over_mean']:.2f}</td><td>{c['E_peak_phantom_vpm']:.0f}</td></tr>" for c in eps)
    rows3 = "".join(
        f"<tr><td>{100 * c['r']:g}</td><td>{c['f_mean_mhz']:.3f}</td><td>{100 * c.get('f_shift_rel', nan):+.2f}</td><td>{c.get('Q_L', nan):.1f}</td>"
        f"<td>{c.get('QU_over_QL', nan):.2f}</td><td>{100 * c['b1plus_phantom_slice_cv']:.1f}</td></tr>" for c in rad)
    items = [("card_b1_unloaded_vs_loaded.png", "Card 1 (1080 x 1350): quadrature |B1+|, empty vs loaded (εr 78, σ 0.6), one colour scale; the bottom row splits the two properties."),
             ("card_resonance_q.png", "Card 2 (1080 x 1350): k = 1 resonance, Q_U vs Q_L on the same frequency axis."),
             ("card_fields_unloaded_vs_loaded.gif", "Card 3 (1080 x 1350): rotating B1 over one RF cycle, empty vs loaded, same scale. MP4: card_fields_unloaded_vs_loaded.mp4."),
             ("loading_q_vs_sigma_fshift_vs_eps.png", "Supporting (1080 x 1080): Q vs σ and frequency shift vs εr."),
             ("loading_b1_strip.png", "Supporting: |B1+| across the σ sweep (top) and the εr sweep (bottom), one colour scale."),
             ("sweep_sigma_b1.gif", "Supporting GIF (1080 x 1080): |B1+| through the σ sweep. MP4: sweep_sigma_b1.mp4."),
             ("sweep_eps_b1.gif", "Supporting GIF (1080 x 1080): |B1+| through the εr sweep. MP4: sweep_eps_b1.mp4."),
             ("loading_e_sar_saline.png", "Supporting (1080 x 1080): |E| and σ|E|²/2 inside the saline phantom."),
             ("loading_sweeps_2x2.png", "Supporting: Q_U/Q_L, shift and homogeneity against σ and εr."),
             ("loading_radius.png", "Supporting: phantom radius 5 to 11 cm at saline values.")]
    figs = "".join(f'<figure><a href="{f}"><img src="{f}" loading="lazy"></a><figcaption><b>{f}</b> &middot; {t}</figcaption></figure>'
                   for f, t in items if (out / f).exists())
    html = f"""<!doctype html><html><head><meta charset="utf-8"><meta name=viewport content="width=device-width,initial-scale=1">
<title>Birdcage loading</title><style>
:root{{--bg:#fff;--ink:#0b0b0b;--ink2:#52514e;--line:#e4e3df}}
@media (prefers-color-scheme:dark){{:root:not([data-theme=light]){{--bg:#151515;--ink:#eee;--ink2:#aaa;--line:#333}}}}
:root[data-theme=dark]{{--bg:#151515;--ink:#eee;--ink2:#aaa;--line:#333}}
body{{background:var(--bg);color:var(--ink);font:15px/1.45 system-ui,sans-serif;margin:0 auto;max-width:1300px;padding:20px 16px}}
h1{{font-size:24px;margin:0 0 4px}} h2{{font-size:18px;margin:22px 0 4px}} .sub{{color:var(--ink2)}}
.tw{{overflow-x:auto}} table{{border-collapse:collapse;font-size:13.5px;margin:6px 0 12px}} td,th{{border-bottom:1px solid var(--line);padding:5px 8px;text-align:right}}
th{{font-size:12px;color:var(--ink2);font-weight:600}}
.grid{{display:grid;grid-template-columns:repeat(auto-fill,minmax(min(100%,360px),1fr));gap:18px}}
figure{{margin:0;border:1px solid var(--line);border-radius:6px;padding:8px;background:#fff}} img{{width:100%;display:block}}
figcaption{{font-size:13px;color:#333;margin-top:6px}}
</style></head><body>
<h1>Birdcage loading: conductivity vs permittivity</h1>
<p class=sub>{HONEST}. {sw['what']}<br>Empty coil (same mesh, phantom region εr 1): k = 1 {u['f_mean_mhz']:.3f} MHz, Q_U {u['Q_L']:.0f}.
Saline (εr 78, σ 0.6 S/m): {l['f_mean_mhz']:.3f} MHz, Q_L {l['Q_L']:.1f}, Q_U/Q_L {l['QU_over_QL']:.1f}.
Data: birdcage-palace/results/loading_sweep.json; method and cost: docs/loading-study.md. Not published; for review.</p>
<div class=grid>{figs}</div>
<h2>σ sweep (εr 78, r 8 cm)</h2><div class=tw><table><tr><th>σ (S/m)</th><th>f (MHz)</th><th>shift (%)</th><th>Q, phantom loss only</th><th>Q_L</th><th>Q_U/Q_L</th><th>B1+ ROI std/mean (%)</th><th>peak |E| (V/m)</th><th>peak σ|E|²/2 (W/m³)</th></tr>{rows}</table></div>
<h2>εr sweep (σ 0, r 8 cm)</h2><div class=tw><table><tr><th>εr</th><th>f (MHz)</th><th>shift (%)</th><th>Q_L</th><th>Q_U/Q_L</th><th>B1+ ROI std/mean (%)</th><th>B1+ z = 0 slice in phantom (%)</th><th>centre / ROI mean</th><th>peak |E| in phantom (V/m)</th></tr>{rows2}</table></div>
<h2>Phantom radius (εr 78, σ 0.6)</h2><div class=tw><table><tr><th>r (cm)</th><th>f (MHz)</th><th>shift (%)</th><th>Q_L</th><th>Q_U/Q_L</th><th>B1+ z = 0 slice in phantom (%)</th></tr>{rows3}</table></div>
<p class=sub>|E| and σ|E|²/2 at 1 µT mean |B1+| over the ROI (r &lt; 8.4 cm, |z| &lt; 4.2 cm), on the z = 0 and y = 0 probe planes (4 mm grid).</p>
</body></html>"""
    (out / "index.html").write_text(html, encoding="utf-8")
    print("wrote", out / "index.html")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--docs", action="store_true", help="only the page/README subset")
    ap.add_argument("--no-anim", action="store_true")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    sw, by, maps = load()
    if a.docs:
        plot_q_fshift(by, out, "loading_q_fshift.png", size=(10.8, 9.0))
        strip_b1(by, maps, out, "loading_b1_strip.png")
        return
    card_b1(by, maps, out)
    card_resonance(by, out)
    plot_q_fshift(by, out)
    plot_supporting(by, out)
    strip_b1(by, maps, out)
    e_maps(by, maps, out)
    if not a.no_anim:
        anim_rotating(by, maps, out)
        anim_sweep(by, maps, out, "sigma")
        anim_sweep(by, maps, out, "eps")
    index_html(out, by, sw)


if __name__ == "__main__":
    main()
