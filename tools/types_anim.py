"""Animated k = 1 mode of low-pass / high-pass / band-pass birdcages, driven in quadrature.

For each solved case the homogeneous k = 1 pair (a, b) is combined as a + i s b (s from modes.quadrature, the sense
that maximises B1+), and the instantaneous field Re{(a + i s b) e^{i w t}} is drawn over one RF cycle:

  * top: the cage in 3D, every rung and end-ring segment coloured by its instantaneous current (Palace capacitor
    currents; the conductors without a capacitor follow from Kirchhoff's current law), arrows along the current
  * bottom: the axial plane z = 0 (Palace field probes, 6 mm grid, r < 0.92 R), colour = instantaneous |B_t|
    normalised to its RMS mean in the central disc r < 0.6 R, white arrows = direction of B_t

    python tools/types_anim.py --case LP=runs/lp16_final --case HP=runs/hp16_tuned --out <dir>/birdcage_types_fields
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Line3DCollection  # noqa: E402
from PIL import Image  # noqa: E402
from scipy.interpolate import griddata  # noqa: E402

from birdcage_palace.modes import load  # noqa: E402

NAMES = {"lowpass": ("Low-pass", "C in the rungs"), "highpass": ("High-pass", "C in the end rings"),
         "bandpass": ("Band-pass", "C in rungs and end rings")}
FIELD_CMAP = "magma"
CUR_CMAP = "RdBu_r"
INK, INK2 = "#0b0b0b", "#52514e"


def pair_fields(case: Path) -> dict:
    """Real, orthonormalised k = 1 pair (fields at the probes and capacitor currents), phase-locked to each other."""
    case = Path(case)
    res = json.loads((case / "modes.json").read_text())
    pr = json.loads((case / "probes.json").read_text())
    info = json.loads((case / "build.json").read_text())
    des = res["design"]
    raw = load(case)
    na, nc = pr["axial"], pr["coronal"]
    roi = slice(na + nc, None)
    pts = np.array(pr["points"])
    ports = info["ports"]
    N, R = des["legs"], des["diameter"] / 2
    phi_r = 2 * np.pi * np.arange(N) / N + info.get("phi0", 0.0)
    hp = res["homogeneous_pair"]
    members = []
    for m in hp["m"][:2]:
        B, I = raw["B"][m - 1].copy(), raw["I"][m - 1].copy()
        v = B[roi].sum(0)
        j = int(np.argmax(np.abs(v)))
        ph = np.exp(-1j * np.angle(v[j]))
        B, I = B * ph, I * ph
        rung, top, bot = conductor_currents(I, ports, N)
        # tie the current phase/sign to the field: thin-wire Biot-Savart at the centre must match the Palace field there
        c0 = B[:na][np.argmin(np.sum(pts[:na, :2] ** 2, 1))][:2]
        pred = np.array([np.sum(rung * np.sin(phi_r)), np.sum(-rung * np.cos(phi_r))])
        c = np.vdot(pred, c0)
        f = np.conj(c) / abs(c)
        members.append({"B": np.real(B), "rung": np.real(rung * f), "top": np.real(top * f), "bot": np.real(bot * f),
                        "imag_frac": float(np.abs(np.imag(B)).max() / np.abs(B).max())})
    a, b = members
    na_ = np.linalg.norm(a["B"][roi])
    for k in ("B", "rung", "top", "bot"):
        a[k] = a[k] / na_
    cab = float(np.sum(a["B"][roi] * b["B"][roi]))
    for k in ("B", "rung", "top", "bot"):
        b[k] = b[k] - cab * a[k]
    nb_ = np.linalg.norm(b["B"][roi])
    for k in ("B", "rung", "top", "bot"):
        b[k] = b[k] / nb_
    best = None
    for s in (1, -1):
        Bq = a["B"] + s * 1j * b["B"]
        bp = np.abs((Bq[roi, 0] + 1j * Bq[roi, 1]) / 2).sum()
        if best is None or bp > best[0]:
            best = (bp, s)
    s = best[1]
    ax_pts = pts[:na, :2]
    j0 = int(np.argmin(np.sum(ax_pts ** 2, 1)))
    cy = a["B"][j0, 1] + s * 1j * b["B"][j0, 1]
    theta = -float(np.angle(cy))
    return {"theta": theta, "case": case, "res": res, "des": des, "info": info, "phi_r": phi_r, "R": R, "N": N,
            "ax_pts": ax_pts, "a": a, "b": b, "s": s, "na": na,
            "imag_frac": max(a.get("imag_frac", 0), b.get("imag_frac", 0))}


def conductor_currents(I: np.ndarray, ports: list, N: int):
    """Rung currents (+z) and top/bottom ring-segment currents (+phi), segment n between rungs n and n+1."""
    kinds = {k: [i for i, p in enumerate(ports) if p["kind"] == k] for k in ("leg", "ring_top", "ring_bot")}
    def by_n(k):
        idx = sorted(kinds[k], key=lambda i: ports[i]["n"])
        sgn = []
        for i in idx:
            p = ports[i]
            ph = np.array([-np.sin(p["phi"]), np.cos(p["phi"]), 0.0])
            sgn.append(np.sign(np.dot(p["dir"], ph)) if k != "leg" else np.sign(p["dir"][2]))
        return I[idx] * np.array(sgn)
    rung = by_n("leg") if kinds["leg"] else None
    top = by_n("ring_top") if kinds["ring_top"] else None
    bot = by_n("ring_bot") if kinds["ring_bot"] else None
    if rung is None:  # KCL at the top node of rung n: I_n = J_n - J_{n-1}
        rung = top - np.roll(top, 1)
    if top is None:
        top = np.cumsum(rung)
        top = top - top.mean()
    if bot is None:
        bot = -top
    return rung, top, bot


def field_grid(d: dict, n: int = 260):
    R = d["R"]
    rmax = 0.88 * R
    g = np.linspace(-R, R, n)
    X, Y = np.meshgrid(g, g)
    out = {}
    for key in ("a", "b"):
        B = d[key]["B"][: d["na"]]
        out[key] = np.stack([griddata(d["ax_pts"], B[:, i], (X, Y), method="cubic") for i in (0, 1)], -1)
    mask = X ** 2 + Y ** 2 > rmax ** 2
    for key in out:
        out[key][mask] = np.nan
    disc = X ** 2 + Y ** 2 < (0.6 * R) ** 2
    rms = np.sqrt((np.sum(out["a"] ** 2, -1) + np.sum(out["b"] ** 2, -1)) / 2)
    d["grid"] = {"X": X, "Y": Y, "g": g, "a": out["a"], "b": out["b"], "norm": float(np.nanmean(rms[disc]))}
    # coarse arrow grid
    q = np.arange(-0.75 * R, 0.75 * R + 1e-9, 0.15 * R)
    QX, QY = np.meshgrid(q, q)
    keep = QX ** 2 + QY ** 2 <= (0.80 * R) ** 2
    qa = np.stack([griddata(d["ax_pts"], d["a"]["B"][: d["na"], i], (QX, QY), method="cubic") for i in (0, 1)], -1)
    qb = np.stack([griddata(d["ax_pts"], d["b"]["B"][: d["na"], i], (QX, QY), method="cubic") for i in (0, 1)], -1)
    d["quiv"] = {"X": QX[keep], "Y": QY[keep], "a": qa[keep], "b": qb[keep]}
    return d


def inst(d, key, wt):
    """Re{(a + i s b) e^{i w t}} = a cos(wt) - s b sin(wt)."""
    return d["a"][key] * np.cos(wt) - d["s"] * d["b"][key] * np.sin(wt) if key in d["a"] else None


def draw_cage(ax, d, wt, cmax, zoom=1.2):
    des = d["des"]
    R, l, N, phr = d["R"] * 100, des["length"] * 100, d["N"], d["phi_r"]
    rung = d["a"]["rung"] * np.cos(wt) - d["s"] * d["b"]["rung"] * np.sin(wt)
    top = d["a"]["top"] * np.cos(wt) - d["s"] * d["b"]["top"] * np.sin(wt)
    bot = d["a"]["bot"] * np.cos(wt) - d["s"] * d["b"]["bot"] * np.sin(wt)
    cmap = plt.get_cmap(CUR_CMAP)
    nrm = Normalize(-cmax, cmax)
    segs, cols = [], []
    for n in range(N):
        x, y = R * np.cos(phr[n]), R * np.sin(phr[n])
        segs.append([(x, y, -l / 2), (x, y, l / 2)])
        cols.append(rung[n])
    for z, J in ((l / 2, top), (-l / 2, bot)):
        for n in range(N):
            t = np.linspace(phr[n], phr[n] + 2 * np.pi / N, 8)
            p = np.c_[R * np.cos(t), R * np.sin(t), np.full_like(t, z)]
            for i in range(len(t) - 1):
                segs.append([tuple(p[i]), tuple(p[i + 1])])
                cols.append(J[n])
    # draw back conductors first (lighter), so the front ones sit on top
    ax.computed_zorder = False
    lc = Line3DCollection(segs, colors=cmap(nrm(np.array(cols))), linewidths=3.4, zorder=1)
    ax.add_collection3d(lc)
    # arrows: rungs at z = +-l/4, ring segments at their midpoints; length ~ current
    L = 0.32 * l / cmax
    xs, ys, zs, us, vs, ws, cs = [], [], [], [], [], [], []
    for n in range(N):
        x, y = R * np.cos(phr[n]), R * np.sin(phr[n])
        for z in (l / 4, -l / 4):
            xs += [x]; ys += [y]; zs += [z - 0.5 * L * rung[n]]; us += [0]; vs += [0]; ws += [L * rung[n]]; cs += [rung[n]]
    cs = np.array(cs)
    big = np.abs(cs) > 0.12 * cmax
    ax.quiver(np.array(xs)[big], np.array(ys)[big], np.array(zs)[big], np.array(us)[big], np.array(vs)[big],
              np.array(ws)[big], color="0.15", linewidth=0.9, arrow_length_ratio=0.4, normalize=False, zorder=2)
    # capacitors
    caps = [(p, d) for p, d in ((p, None) for p in d["info"]["ports"])]
    cx, cy, cz = [], [], []
    for p, _ in caps:
        ph = p["phi"] + (d["info"].get("phi0", 0.0) if p["kind"] == "leg" and False else 0.0)
        z = 0.0 if p["kind"] == "leg" else (l / 2 if p["kind"] == "ring_top" else -l / 2)
        cx.append(R * np.cos(ph)); cy.append(R * np.sin(ph)); cz.append(z)
    ax.scatter(cx, cy, cz, marker="s", s=14, c="#f2d13b", edgecolors="0.1", linewidths=0.5, depthshade=False, zorder=3)
    # slice plane outline
    t = np.linspace(0, 2 * np.pi, 120)
    ax.plot(0.88 * R * np.cos(t), 0.88 * R * np.sin(t), 0 * t, color="0.45", lw=0.7, ls="--", zorder=0)
    lim = R * 1.02
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-l / 2 * 1.15, l / 2 * 1.15)
    ax.set_box_aspect((1, 1, (l * 1.15) / (2 * lim)), zoom=zoom)
    ax.view_init(elev=22, azim=-62)
    ax.set_axis_off()


def draw_slice(ax, d, wt, vmax=2.0, width=0.008):
    G = d["grid"]
    Bx = G["a"][..., 0] * np.cos(wt) - d["s"] * G["b"][..., 0] * np.sin(wt)
    By = G["a"][..., 1] * np.cos(wt) - d["s"] * G["b"][..., 1] * np.sin(wt)
    mag = np.sqrt(Bx ** 2 + By ** 2) / G["norm"]
    R = d["R"] * 100
    im = ax.imshow(mag, origin="lower", extent=(-R, R, -R, R), cmap=FIELD_CMAP, vmin=0, vmax=vmax,
                   interpolation="bilinear")
    Q = d["quiv"]
    qx = Q["a"][:, 0] * np.cos(wt) - d["s"] * Q["b"][:, 0] * np.sin(wt)
    qy = Q["a"][:, 1] * np.cos(wt) - d["s"] * Q["b"][:, 1] * np.sin(wt)
    qm = np.hypot(qx, qy) / G["norm"]
    sc = np.minimum(qm, 1.6) / np.maximum(qm, 1e-12)
    ax.quiver(Q["X"] * 100, Q["Y"] * 100, qx * sc / G["norm"], qy * sc / G["norm"], color="white", scale=1 / 1.75,
              scale_units="xy", angles="xy", width=width, headwidth=3.4, headlength=3.6, pivot="middle", alpha=0.95)
    # rungs at the coil radius, coloured by their instantaneous current
    rung = d["a"]["rung"] * np.cos(wt) - d["s"] * d["b"]["rung"] * np.sin(wt)
    cm = d["cmax"]
    ax.scatter(R * np.cos(d["phi_r"]), R * np.sin(d["phi_r"]), s=34, c=rung, cmap=CUR_CMAP, vmin=-cm, vmax=cm,
               edgecolors="0.15", linewidths=0.6, zorder=5)
    t = np.linspace(0, 2 * np.pi, 200)
    ax.plot(0.6 * R * np.cos(t), 0.6 * R * np.sin(t), ":", color="white", lw=0.8, alpha=0.7)
    ax.set_xlim(-R * 1.08, R * 1.08); ax.set_ylim(-R * 1.08, R * 1.08)
    ax.set_aspect("equal")
    ax.axis("off")
    return im


def label(d):
    des, res = d["des"], d["res"]
    name, where = NAMES[des["coil_type"]]
    hp = res["homogeneous_pair"]
    c = []
    if des.get("c_leg"):
        c.append(f"C rung {des['c_leg'] * 1e12:.2f} pF")
    if des.get("c_ring"):
        c.append(f"C ring {des['c_ring'] * 1e12:.2f} pF")
    st = res.get("solve_stats") or {}
    o = des.get("order")
    kind = "production mesh" if des.get("h_coil", 1) <= 0.0051 and o == 2 else "coarse mesh"
    ml = f"{kind} · {st.get('unknowns', 0) / 1e6:.2f} M unknowns" if st.get("unknowns") else kind
    return name, where, f"k = 1 at {hp['f_mean_mhz']:.2f} MHz", " · ".join(c), ml


def render(datas, out: Path, nfr: int = 30, width_in: float = 5.4, dpi: int = 200):
    out.parent.mkdir(parents=True, exist_ok=True)
    fdir = out.parent / "_work" / (out.name + "_frames")
    fdir.mkdir(parents=True, exist_ok=True)
    for f in fdir.glob("*.png"):
        f.unlink()
    ncol = len(datas)
    for d in datas:
        d["cmax"] = max(np.max(np.abs(d["a"][k] + 1j * d["b"][k])) for k in ("rung", "top", "bot"))
    if ncol == 1:  # 4:5, cage above slice
        H = width_in * 1.25
        L = dict(title=0.966, sub=0.936, name=0.895, where=0.872, cage=[0.10, 0.585, 0.80, 0.27],
                 sl=[0.20, 0.215, 0.60, 0.36], f=0.192, c=0.170, mesh=0.150, cb=0.085, fs=(12, 9, 9, 8))
    else:  # square, one column per coil
        H = width_in
        L = dict(title=0.962, sub=0.924, name=0.872, where=0.845, cage=[0, 0.535, 1, 0.30],
                 sl=[0.005, 0.235, -0.01, 0.30], f=0.205, c=0.180, mesh=0.163, cb=0.105, fs=(10.5, 7.6, 7.6, 6.6))
    paths = []
    for fi in range(nfr):
        wt = 2 * np.pi * fi / nfr
        fig = plt.figure(figsize=(width_in, H), dpi=dpi, facecolor="white")
        title = ("16-rung birdcage at 3 T: the k = 1 mode in quadrature" if ncol == 1
                 else "16-rung birdcage at 3 T: where the capacitors go")
        fig.text(0.04, L["title"], title, fontsize=11.5, weight="bold", va="top", color=INK)
        fig.text(0.04, L["sub"], "k = 1 pair driven 90° apart, one RF cycle. Top: currents. Bottom: B in the axial "
                 "plane z = 0.", fontsize=6.6, color=INK2, va="top")
        w = 0.94 / ncol
        fs = L["fs"]
        for i, d in enumerate(datas):
            x0 = 0.03 + i * w
            name, where, fl, cl, ml = label(d)
            fig.text(x0 + w / 2, L["name"], name, ha="center", fontsize=fs[0], weight="bold", color=INK)
            fig.text(x0 + w / 2, L["where"], where, ha="center", fontsize=fs[1], color=INK2)
            if ncol == 1:
                a3 = fig.add_axes(L["cage"], projection="3d")
                a2 = fig.add_axes(L["sl"])
            else:
                c3, s2 = L["cage"], L["sl"]
                a3 = fig.add_axes([x0 + c3[0] * w, c3[1], c3[2] * w, c3[3]], projection="3d")
                a2 = fig.add_axes([x0 + s2[0], s2[1], w + s2[2], s2[3]])
            a3.set_facecolor("white")
            draw_cage(a3, d, wt + d["theta"], d["cmax"], zoom=1.4 if ncol == 1 else 1.2)
            draw_slice(a2, d, wt + d["theta"], width=0.006 if ncol == 1 else 0.009)
            fig.text(x0 + w / 2, L["f"], fl, ha="center", fontsize=fs[2], color=INK, weight="bold")
            fig.text(x0 + w / 2, L["c"], cl, ha="center", fontsize=fs[3], color=INK2)
            fig.text(x0 + w / 2, L["mesh"], ml, ha="center", fontsize=fs[3] - 0.8, color="#8a8984")
        cb_y = L["cb"]
        cax = fig.add_axes([0.08, cb_y, 0.36, 0.012 * 5.4 * 1.25 / H])
        cb = fig.colorbar(plt.cm.ScalarMappable(Normalize(0, 2.0), FIELD_CMAP), cax=cax, orientation="horizontal")
        cb.set_ticks([0, 0.5, 1, 1.5, 2]); cb.ax.tick_params(labelsize=5.8, length=2, pad=1)
        cb.outline.set_linewidth(0.4)
        fig.text(0.08, cb_y + 0.022 * 6.75 / H, "|B| (instantaneous) / RMS mean over r < 0.6 R", fontsize=5.8, color=INK2)
        cax2 = fig.add_axes([0.56, cb_y, 0.36, 0.012 * 5.4 * 1.25 / H])
        cb2 = fig.colorbar(plt.cm.ScalarMappable(Normalize(-1, 1), CUR_CMAP), cax=cax2, orientation="horizontal")
        cb2.set_ticks([-1, 0, 1]); cb2.ax.tick_params(labelsize=5.8, length=2, pad=1)
        cb2.outline.set_linewidth(0.4)
        fig.text(0.56, cb_y + 0.022 * 6.75 / H, "current / coil max (rungs +z, rings anticlockwise)", fontsize=5.8,
                 color=INK2)
        fig.text(0.08, 0.045 * 6.75 / H, "Yellow squares: capacitors. Dotted circle: r = 0.6 R. White arrows: B "
                 "direction.", fontsize=5.8, color=INK2)
        fig.text(0.08, 0.020 * 6.75 / H, "Palace FEM · unloaded · simulation only", fontsize=6.8, color=INK2,
                 style="italic")
        p = fdir / f"f{fi:03d}.png"
        fig.savefig(p, dpi=dpi, facecolor="white")
        plt.close(fig)
        paths.append(p)
    return paths


def export(paths, out: Path, width=1080, fps=10, colors=96, cycles_mp4=4):
    imgs = [Image.open(p).convert("RGB") for p in paths]
    h = int(round(imgs[0].height * width / imgs[0].width))
    imgs = [im.resize((width, h), Image.LANCZOS) for im in imgs]
    imgs[0].save(out.with_name(out.name + "_still.png"))
    pal = imgs[0].quantize(colors=colors, method=Image.Quantize.MEDIANCUT)
    small = [im.quantize(palette=pal, dither=Image.Dither.NONE) for im in imgs]
    gif = out.with_suffix(".gif")
    small[0].save(gif, save_all=True, append_images=small[1:], duration=int(1000 / fps), loop=0, optimize=True)
    fdir = paths[0].parent
    mp4 = out.with_suffix(".mp4")
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-stream_loop", str(cycles_mp4 - 1), "-framerate", str(fps),
                    "-i", str(fdir / "f%03d.png"), "-vf", f"scale={width}:{h + h % 2}:flags=lanczos,format=yuv420p",
                    "-c:v", "libx264", "-crf", "18", "-movflags", "+faststart", str(mp4)], check=True)
    return gif, mp4, gif.stat().st_size / 1e6


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", action="append", required=True, help="LABEL=path (order = panel order)")
    ap.add_argument("--out", required=True, help="output path without extension")
    ap.add_argument("--frames", type=int, default=30)
    ap.add_argument("--fps", type=int, default=10)
    ap.add_argument("--colors", type=int, default=96)
    a = ap.parse_args()
    datas = [field_grid(pair_fields(Path(c.split("=", 1)[1]))) for c in a.case]
    for d in datas:
        print(d["case"], "s =", d["s"], "imag frac", f"{d['imag_frac']:.1e}", "norm", f"{d['grid']['norm']:.3g}")
    paths = render(datas, Path(a.out), nfr=a.frames)
    gif, mp4, mb = export(paths, Path(a.out), fps=a.fps, colors=a.colors)
    print(f"wrote {gif} ({mb:.2f} MB), {mp4}")


if __name__ == "__main__":
    main()
