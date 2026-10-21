"""Stills comparing low-pass / high-pass / band-pass birdcages (one solved case each).

    python tools/types_figs.py --lp runs/lp16_final --hp runs/hp16_3T --bp runs/bp16_3T --out <dir>

writes coils_types.png (a), spectra_types.png (b), b1_types.png (c), netlist_text.png (backup of the .cnl as text)
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from PIL import Image  # noqa: E402
from scipy.interpolate import griddata  # noqa: E402

from birdcage_palace import plots  # noqa: E402
from birdcage_palace.modes import load, quadrature  # noqa: E402

NAMES = {"lowpass": ("Low-pass", "C in the rungs"), "highpass": ("High-pass", "C in the end rings"),
         "bandpass": ("Band-pass", "C in rungs and end rings")}
INK, INK2, MUTED = plots.INK, plots.INK2, plots.MUTED
CAPC = "#f2d13b"


def ctx(case: Path) -> dict:
    case = Path(case)
    res = json.loads((case / "modes.json").read_text())
    return {"case": case, "res": res, "des": res["design"], "an": json.loads((case / "analytic.json").read_text()),
            "info": json.loads((case / "build.json").read_text())}


def cvals(des) -> str:
    c = []
    if des.get("c_leg"):
        c.append(f"rung {des['c_leg'] * 1e12:.2f} pF")
    if des.get("c_ring"):
        c.append(f"ring {des['c_ring'] * 1e12:.2f} pF")
    return "C " + " · ".join(c)


def mesh_label(c) -> str:
    des, st = c["des"], c["res"].get("solve_stats") or {}
    kind = "production mesh" if des.get("h_coil", 1) <= 0.0051 and des.get("order") == 2 else "coarse mesh"
    return f"{kind}, order {des.get('order')}, {st.get('unknowns', 0) / 1e6:.2f} M unknowns"


# ------------------------------------------------------------------------------------------------------- (a) coils
def render_coil(c, out: Path, window=(1400, 1400)) -> Path:
    import pyvista as pv
    case = c["case"]
    s = json.loads((case / "surfaces.json").read_text())
    des = c["des"]
    pl = pv.Plotter(off_screen=True, window_size=window)
    pl.set_background("white")

    def mesh(lbl):
        d = s[lbl]
        if not d["nodes"]:
            return None
        tri = np.array(d["tris"]).reshape(-1, 3)
        return pv.PolyData(np.array(d["nodes"]) * 100, np.hstack([np.full((len(tri), 1), 3), tri]).ravel())
    pl.add_mesh(mesh("coil"), color="#d98a4a", specular=0.5, specular_power=20, ambient=0.25)
    caps = mesh("caps")
    if caps is not None:
        pl.add_mesh(caps, color=plots.BLUE)
    R, l = des["diameter"] / 2 * 100, des["length"] * 100
    pts = []
    for p in c["info"]["ports"]:
        z = 0.0 if p["kind"] == "leg" else (l / 2 if p["kind"] == "ring_top" else -l / 2)
        pts.append([R * np.cos(p["phi"]), R * np.sin(p["phi"]), z])
    glyph = pv.PolyData(np.array(pts)).glyph(geom=pv.Sphere(radius=1.15, theta_resolution=24, phi_resolution=24),
                                             scale=False, orient=False)
    pl.add_mesh(glyph, color=CAPC, specular=0.3, ambient=0.3)
    Rs = des["shield_diameter"] / 2 * 100
    Ls = (des.get("shield_length") or des["length"] + 0.20) * 100
    az = np.arctan2(-60.0, 95.0)
    th = np.linspace(az + np.pi / 2, az + 3 * np.pi / 2, 90)
    T, Z = np.meshgrid(th, np.array([-Ls / 2, Ls / 2]))
    pl.add_mesh(pv.StructuredGrid(Rs * np.cos(T), Rs * np.sin(T), Z), color="#b8bcc4", opacity=0.22, smooth_shading=True)
    pl.camera_position = [(95, -60, 55), (0, 0, 0), (0, 0, 1)]
    pl.camera.zoom(1.3)
    pl.enable_anti_aliasing("ssaa")
    pl.screenshot(str(out))
    pl.close()
    return out


def fig_coils(cs, out: Path, work: Path):
    imgs = [render_coil(c, work / f"coilcap_{c['des']['coil_type']}.png") for c in cs]
    fig = plt.figure(figsize=(12, 6.2), dpi=150, facecolor="white")
    fig.text(0.03, 0.95, "16-rung birdcage at 3 T: three places to put the capacitors", fontsize=17, weight="bold",
             color=INK, va="top")
    fig.text(0.03, 0.885, "Same cage (28 cm diameter, 25 cm long, 34 cm shield). Yellow: capacitor gaps. Values tune "
             "the k = 1 mode to 127.7 MHz (Palace eigenmode).", fontsize=10.5, color=INK2, va="top")
    for i, (c, p) in enumerate(zip(cs, imgs)):
        ax = fig.add_axes([0.01 + i * 0.33, 0.14, 0.32, 0.64])
        im = Image.open(p)
        w, h = im.size
        ax.imshow(im.crop((int(0.06 * w), int(0.04 * h), int(0.94 * w), int(0.92 * h))))
        ax.axis("off")
        name, where = NAMES[c["des"]["coil_type"]]
        n = len(c["info"]["ports"])
        fig.text(0.17 + i * 0.33, 0.80, name, ha="center", fontsize=15, weight="bold", color=INK)
        fig.text(0.17 + i * 0.33, 0.115, f"{where} ({n} capacitors)", ha="center", fontsize=11, color=INK)
        fig.text(0.17 + i * 0.33, 0.075, cvals(c["des"]) + f" · k = 1 at {c['res']['homogeneous_pair']['f_mean_mhz']:.2f} MHz",
                 ha="center", fontsize=10, color=INK2)
    fig.text(0.97, 0.02, "Palace FEM · unloaded · simulation only", ha="right", fontsize=9.5, color=INK2, style="italic")
    fig.savefig(out, dpi=150, facecolor="white")
    return out


# ----------------------------------------------------------------------------------------------------- (b) spectra
def fig_spectra(cs, out: Path):
    fig, axs = plt.subplots(len(cs), 1, figsize=(12, 2.35 * len(cs) + 1.3), dpi=150, sharex=True, facecolor="white")
    f0 = cs[0]["an"]["f_target_hz"] / 1e6
    fmax = 0
    for ax, c in zip(axs, cs):
        res = c["res"]
        name, where = NAMES[c["des"]["coil_type"]]
        rows = res["comparison"]
        pts_ = []
        for r in rows:
            fp = r["f_palace_mhz"]
            if not fp:
                continue
            k = r["k"]
            fpm = float(np.mean(fp))
            fl = r["f_ladder_mhz"]
            fmax = max(fmax, fpm, fl)
            ax.plot([fpm], [1.0], "o", ms=8, color=plots.BLUE, zorder=3)
            ax.plot([fl], [0.0], "D", ms=6.5, mfc="white", mec=plots.ORANGE, mew=1.6, zorder=3)
            ax.plot([fl, fpm], [0.0, 1.0], color=MUTED, lw=0.7, zorder=1)
            pts_.append((fpm, k, f"{k}" if not (k == 0 and c["des"]["coil_type"] != "lowpass") else "0*"))
        # labels: modes closer than 6 MHz form a cluster; a cluster of more than two gets one "a-b" label
        pts_.sort()
        clusters, cur = [], []
        for q in pts_:
            if cur and q[0] - cur[-1][0] > 6.0:
                clusters.append(cur); cur = []
            cur.append(q)
        if cur:
            clusters.append(cur)
        for cl in clusters:
            if len(cl) > 2:
                ks = sorted(q[1] for q in cl)
                ax.text(np.mean([q[0] for q in cl]), 1.22, f"{ks[0]}-{ks[-1]}", ha="center", va="bottom", fontsize=9,
                        color=INK)
            else:
                for i, (fpm, k, lab) in enumerate(cl):
                    ax.text(fpm + (i - (len(cl) - 1) / 2) * 7, 1.22, lab, ha="center", va="bottom", fontsize=9,
                            color=INK, weight="bold" if k == 1 else "normal")
        # Palace modes that are not in the ladder comparison (cap-free, co-rotating ring, other)
        used = {round(f, 6) for r in rows for f in r["f_palace_mhz"]}
        other = [m for m in res["modes"] if round(m["f_mhz"], 6) not in used]
        last_o, lev_o, last_t = -1e9, 0, None
        for m in sorted(other, key=lambda m: m["f_mhz"]):
            fmax = max(fmax, m["f_mhz"])
            ax.plot([m["f_mhz"]], [1.0], "x", ms=7, mew=1.5, color=MUTED, zorder=2)
            if m.get("cap_free"):
                t = "s"
            elif m.get("k") == 0 and m.get("ring_current_frac", 1) < 0.5:
                t = "0"
            elif m.get("ring_bot_over_top", 1) < 0 and m.get("k") is not None:
                t = f"{m['k']}c"
            elif m.get("k") == 0:
                t = "0"
            else:
                t = "?"
            if m["f_mhz"] - last_o < 6.0 and t == last_t:
                continue  # one label per degenerate pair
            lev_o = (lev_o + 1) % 2 if m["f_mhz"] - last_o < 6.0 else 0
            last_o, last_t = m["f_mhz"], t
            ax.text(m["f_mhz"] + 5 * lev_o, 0.78, t, ha="center", va="top", fontsize=8.5, color=MUTED)
        if c["an"]["spectrum"].get("f_corot_hz"):
            fc = c["an"]["spectrum"]["f_corot_hz"] / 1e6
            ax.plot([fc], [0.0], "D", ms=6.5, mfc="white", mec=MUTED, mew=1.4)
        ax.axvline(f0, color=plots.ORANGE, lw=1.2, ls="--", zorder=0)
        ax.set_ylim(-0.6, 2.05)
        ax.set_yticks([0, 1])
        ax.set_yticklabels(["ladder", "Palace"], fontsize=9.5)
        ax.text(0.0, 1.02, f"{name}  ", transform=ax.transAxes, fontsize=12, weight="bold", color=INK, va="bottom")
        ax.text({"Low-pass": 0.088, "High-pass": 0.093, "Band-pass": 0.097}[name], 1.02, f"{where} · {cvals(c['des'])} · {mesh_label(c)}",
                transform=ax.transAxes, fontsize=9.2, color=INK2, va="bottom")
        ax.grid(axis="x", color=plots.GRID, lw=0.6)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
    axs[-1].set_xlabel("frequency (MHz)")
    axs[-1].set_xlim(0, fmax * 1.04)
    axs[0].text(f0, 1.9, " 127.7 MHz (1H at 3 T)", color=plots.ORANGE, fontsize=9, va="top")
    fig.suptitle("Mode spectra: where k = 1 sits for each capacitor placement", x=0.03, ha="left", y=0.985,
                 fontsize=15, weight="bold", color=INK)
    fig.text(0.03, 0.935, "Numbers above the markers: mode index k (pairs averaged). Filled circles: Palace eigenmodes. "
             "Open diamonds: ladder model at the same C (grey: co-rotating end-ring mode).\n"
             "Grey crosses (label below): Palace modes outside the ladder model. kc: co-rotating end-ring family with k periods; "
             "s: no capacitor current (stray-C resonance);\n0: rungs in phase. 0*: end-ring k = 0, counter-rotating.", fontsize=9, color=INK2, va="top")
    fig.text(0.97, 0.01, "Palace FEM · unloaded · simulation only", ha="right", fontsize=9, color=INK2, style="italic")
    fig.subplots_adjust(left=0.07, right=0.98, top=0.82, bottom=0.08, hspace=0.55)
    fig.savefig(out, dpi=150, facecolor="white")
    return out


# --------------------------------------------------------------------------------------------------- (c) B1 maps
def fig_b1(cs, out: Path, lo=0.8, hi=1.2):
    fig = plt.figure(figsize=(12, 5.4), dpi=150, facecolor="white")
    fig.text(0.03, 0.955, "Quadrature |B1+| in the axial plane z = 0", fontsize=16, weight="bold", color=INK, va="top")
    fig.text(0.03, 0.895, "k = 1 pair combined 90° apart. Each map divided by its own mean over r < 0.6 R in this plane; "
             "same colour scale for all three.", fontsize=10.2, color=INK2, va="top")
    for i, c in enumerate(cs):
        res, case = c["res"], c["case"]
        pr = json.loads((case / "probes.json").read_text())
        raw = load(case)
        na, nc = pr["axial"], pr["coronal"]
        roi = slice(na + nc, None)
        hp = res["homogeneous_pair"]
        q = quadrature(raw["B"][hp["m"][0] - 1], raw["B"][hp["m"][1] - 1], roi)
        pts = np.array(pr["points"])[:na, :2]
        v = np.abs(q["b1p"][:na])
        R = c["des"]["diameter"] / 2
        g = np.linspace(-R, R, 300)
        X, Y = np.meshgrid(g, g)
        img = griddata(pts, v, (X, Y), method="cubic")
        img[X ** 2 + Y ** 2 > (0.88 * R) ** 2] = np.nan
        img = img / np.nanmean(img[X ** 2 + Y ** 2 < (0.6 * R) ** 2])
        ax = fig.add_axes([0.03 + i * 0.31, 0.17, 0.27, 0.62])
        im = ax.imshow(img, origin="lower", extent=(-R * 100, R * 100, -R * 100, R * 100), cmap="magma", vmin=lo,
                       vmax=hi, interpolation="bilinear")
        ci = np.where(X ** 2 + Y ** 2 < (0.8 * R) ** 2, img, np.nan)
        ax.contour(g * 100, g * 100, ci, levels=[0.95, 1.05], colors="white", linewidths=0.8)
        phr = 2 * np.pi * np.arange(c["des"]["legs"]) / c["des"]["legs"] + c["info"].get("phi0", 0.0)
        ax.plot(R * 100 * np.cos(phr), R * 100 * np.sin(phr), "o", ms=4, color=plots.ORANGE, mec="0.2", mew=0.5)
        t = np.linspace(0, 2 * np.pi, 200)
        ax.plot(60 * R * np.cos(t), 60 * R * np.sin(t), ":", color="white", lw=0.9)
        ax.set_aspect("equal")
        ax.axis("off")
        name, where = NAMES[c["des"]["coil_type"]]
        ax.set_title(f"{name}\n{where}", fontsize=12, color=INK)
        ax.text(0.5, -0.06, f"3D ROI std/mean of |B1+|: {hp['b1plus_roi_cv'] * 100:.1f} %  ·  "
                f"|B1-|/|B1+|: {hp['b1minus_over_b1plus_roi'] * 100:.1f} %", transform=ax.transAxes, ha="center",
                va="top", fontsize=9, color=INK)
        ax.text(0.5, -0.13, mesh_label(c), transform=ax.transAxes, ha="center", va="top", fontsize=8, color=MUTED)
    cax = fig.add_axes([0.945, 0.22, 0.012, 0.52])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("normalized |B1+|", fontsize=9)
    cb.ax.tick_params(labelsize=8)
    fig.text(0.03, 0.025, "Orange dots: rungs. Dotted circle: r = 0.6 R. White lines: 0.95 and 1.05 (drawn for r < 0.8 R). 3D ROI: r < 0.6 R, |z| < L/6.", fontsize=9,
             color=INK2)
    fig.text(0.97, 0.025, "Palace FEM · unloaded · simulation only", ha="right", fontsize=9, color=INK2, style="italic")
    fig.savefig(out, dpi=150, facecolor="white")
    return out


def fig_netlist(cnl: Path, out: Path):
    lines = cnl.read_text().splitlines()
    fig = plt.figure(figsize=(8.5, 0.205 * len(lines) + 0.7), dpi=170, facecolor="#fbfbfa")
    fig.text(0.02, 1 - 0.18 / (0.205 * len(lines) + 0.7), f"{cnl.parent.parent.name}/circuitrf/{cnl.name}",
             fontsize=8.5, color=INK2, family="monospace", weight="bold")
    for i, ln in enumerate(lines):
        y = 1 - (0.5 + 0.205 * i) / (0.205 * len(lines) + 0.7)
        if ln.startswith(";"):
            col = MUTED
        elif " = " in ln:
            col = plots.BLUE
        else:
            col = INK
        fig.text(0.02, y, ln, fontsize=7.6, color=col, family="monospace", va="top")
    fig.savefig(out, dpi=170, facecolor="#fbfbfa")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lp", required=True)
    ap.add_argument("--hp", required=True)
    ap.add_argument("--bp", default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cnl", default=None)
    ap.add_argument("--only", default="a,b,c")
    a = ap.parse_args()
    out = Path(a.out)
    work = out / "_work"
    work.mkdir(parents=True, exist_ok=True)
    cs = [ctx(Path(p)) for p in (a.lp, a.hp, a.bp) if p]
    only = a.only.split(",")
    if "a" in only:
        print(fig_coils(cs, out / "coils_types.png", work))
    if "b" in only:
        print(fig_spectra(cs, out / "spectra_types.png"))
    if "c" in only:
        print(fig_b1(cs, out / "b1_types.png"))
    if a.cnl:
        print(fig_netlist(Path(a.cnl), out / "porttune_netlist_text.png"))


if __name__ == "__main__":
    main()
