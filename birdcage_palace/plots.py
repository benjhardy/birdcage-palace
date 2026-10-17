"""Figures: mode spectrum, |B1| maps, rung-current patterns, 3D coil/mesh render."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

from .modes import load, quadrature  # noqa: E402

INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
SEQ = LinearSegmentedColormap.from_list("seq_blue", ["#f7f9fc", "#9cc3ee", "#2a78d6", "#123f80", "#061a3a"])
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.titlesize": 11, "axes.titleweight": "bold",
                     "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 200,
                     "savefig.facecolor": "white", "legend.fontsize": 8})


def _ctx(case: Path):
    case = Path(case)
    res = json.loads((case / "modes.json").read_text())
    pr = json.loads((case / "probes.json").read_text())
    raw = load(case)
    pts = np.array(pr["points"])
    na, nc = pr["axial"], pr["coronal"]
    return res, raw, pts, {"axial": slice(0, na), "coronal": slice(na, na + nc), "roi": slice(na + nc, None)}


def _grid(pts2: np.ndarray, vals: np.ndarray):
    xs, ys = np.unique(np.round(pts2[:, 0], 5)), np.unique(np.round(pts2[:, 1], 5))
    img = np.full((len(ys), len(xs)), np.nan)
    ix = np.searchsorted(xs, np.round(pts2[:, 0], 5))
    iy = np.searchsorted(ys, np.round(pts2[:, 1], 5))
    img[iy, ix] = vals
    h = xs[1] - xs[0]
    return img, [xs[0] - h / 2, xs[-1] + h / 2, ys[0] - h / 2, ys[-1] + h / 2]


def pick_modes(res: dict) -> list[dict]:
    """Representative modes: homogeneous k=1, k=2, k=N/2 and the k=0 / end-ring mode (where found)."""
    modes = res["modes"]
    N = res["N"]
    out = []
    hom = [m for m in modes if m["homogeneous"]]
    if hom:
        out.append(hom[0])
    for k in (2, N // 2, 0):
        c = [m for m in modes if m.get("k") == k and m.get("k_purity", 0) > 0.8 and not m["homogeneous"]]
        if c:
            out.append(sorted(c, key=lambda m: -m["k_purity"])[0])
    return out


def _bt(B):
    return np.sqrt(np.abs(B[:, 0]) ** 2 + np.abs(B[:, 1]) ** 2)


def _coil_marks(ax, des, plane):
    R, l, N = des["diameter"] / 2, des["length"], des["legs"]
    if plane == "axial":
        ph = 2 * np.pi * np.arange(N) / N
        ax.plot(R * np.cos(ph) * 100, R * np.sin(ph) * 100, "o", ms=3.5, color=ORANGE, mec="white", mew=0.6)
        t = np.linspace(0, 2 * np.pi, 200)
        ax.plot(0.6 * R * np.cos(t) * 100, 0.6 * R * np.sin(t) * 100, "--", lw=0.8, color="white", alpha=0.8)
    else:
        for x in (-R, R):
            ax.plot([x * 100] * 2, [-l / 2 * 100, l / 2 * 100], color=ORANGE, lw=1.6)
        for z in (-l / 2, l / 2):
            ax.plot([-R * 100, R * 100], [z * 100] * 2, color=ORANGE, lw=0.8, ls=":")
        ax.add_patch(plt.Rectangle((-0.6 * R * 100, -l / 6 * 100), 1.2 * R * 100, l / 3 * 100, fill=False,
                                   ec="white", ls="--", lw=0.8, alpha=0.8))


def fig_b1_maps(case: Path, out: Path, quad: bool = True) -> Path:
    res, raw, pts, sl = _ctx(case)
    des = res["design"]
    sel = pick_modes(res)
    cols = len(sel) + (1 if quad and res.get("homogeneous_pair", {}).get("m", []) and len(res["homogeneous_pair"]["m"]) == 2 else 0)
    fig, axs = plt.subplots(2, cols, figsize=(3.0 * cols, 6.6), gridspec_kw={"height_ratios": [1, 1.25]})
    axs = np.atleast_2d(axs).reshape(2, cols)
    panels = []
    for m in sel:
        B = raw["B"][m["m"] - 1]
        title = m["label"] + (" (homogeneous)" if m["homogeneous"] else "")
        panels.append((title, f"{m['f_mhz']:.1f} MHz", _bt(B), "|B$_t$|"))
    if cols > len(sel):
        a, b = res["homogeneous_pair"]["m"]
        q = quadrature(raw["B"][a - 1], raw["B"][b - 1], sl["roi"])
        panels.insert(1, ("k=1 pair, quadrature", "B$_1^+$ = (B$_x$+iB$_y$)/2", np.abs(q["b1p"]), "|B$_1^+$|"))
    for j, (title, sub, v, lab) in enumerate(panels):
        vmax = np.nanmax(v[sl["axial"]].tolist() + v[sl["coronal"]].tolist())
        for i, (plane, cut) in enumerate((("axial", (0, 1)), ("coronal", (0, 2)))):
            ax = axs[i, j]
            p = pts[sl[plane]][:, cut]
            img, ext = _grid(p, v[sl[plane]] / vmax)
            ext = [e * 100 for e in ext]
            im = ax.imshow(img, origin="lower", extent=ext, cmap=SEQ, vmin=0, vmax=1, interpolation="bilinear")
            _coil_marks(ax, des, plane)
            ax.set_aspect("equal")
            ax.tick_params(labelsize=7)
            ax.set_xlabel("x (cm)", fontsize=8)
            if j == 0:
                ax.set_ylabel("y (cm)" if plane == "axial" else "z (cm)", fontsize=8)
            for s in ax.spines.values():
                s.set_visible(False)
            if i == 0:
                ax.set_title(f"{title}\n{sub}", fontsize=9.5, color=ORANGE if "homogeneous" in title or "quadrature" in title else INK)
    cb = fig.colorbar(im, ax=axs, shrink=0.55, pad=0.015)
    cb.set_label("field magnitude / plane max", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    fig.suptitle(f"Transverse |B| per mode: axial plane z=0 (top) and coronal plane y=0 (bottom); dashed = ROI",
                 fontsize=10, color=INK2, x=0.45)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def _spectrum_axes(ax, res, title=None, legend=True, label_target=True, legend_loc=None):
    an = res["analytic"]["spectrum"]
    N = res["N"]
    lowpass = res["design"]["coil_type"] == "lowpass"
    f0 = res["design"]["b0"] * 42.577478
    sc = plt.rcParams["font.size"] / 10.0
    ks = np.array(an["k"])
    fa = np.array(an["f_hz"]) / 1e6
    keep = fa > 0
    ax.plot(ks[keep], fa[keep], "-", color=MUTED, lw=1.2, zorder=1)
    ax.plot(ks[keep], fa[keep], "o", ms=8 * sc, mfc="white", mec=INK2, mew=1.4, zorder=2, label="ladder model (analytic)")
    other = []
    for m in res["modes"]:
        if m.get("k_purity", 0) >= 0.8 and m.get("k") is not None:
            ax.plot(m["k"], m["f_mhz"], "o", ms=6 * sc, color=ORANGE if m["homogeneous"] else BLUE, mec="white", mew=0.8, zorder=3)
        else:
            other.append(m)
    ax.plot([], [], "o", ms=6, color=BLUE, label="Palace eigenmodes, k from rung currents")
    ax.plot([], [], "o", ms=6, color=ORANGE, label="homogeneous k=1 pair")
    if lowpass:
        k0 = [m for m in res["modes"] if m.get("k") == 0 and m.get("k_purity", 0) >= 0.8]
        if k0:
            ax.annotate("rungs in phase,\nshield return", (0, k0[0]["f_mhz"]), xytext=(0.25, k0[0]["f_mhz"] - 8),
                        fontsize=7 * sc, color=INK2, va="top")
    ax.axhline(f0, color=ORANGE, lw=0.8, ls="--", zorder=0)
    if label_target:
        ax.text(-0.35, f0 + 4, f"{f0:.1f} MHz (1H at {res['design']['b0']:g} T)", va="bottom", fontsize=8 * sc, color=ORANGE)
    if other:
        ax.plot([N / 2 + 1.1] * len(other), [m["f_mhz"] for m in other], "x", color=MUTED, ms=6 * sc, mew=1.4 * sc,
                label="Palace modes with no capacitor current")
    ax.set_xticks(list(range(N // 2 + 1)) + ([N / 2 + 1.1] if other else []))
    ax.set_xticklabels([str(k) for k in range(N // 2 + 1)] + (["other"] if other else []))
    ax.set_xlabel("mode index k (rung current ∝ cos 2πkn/N)")
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_xlim(-0.5, N / 2 + 1.9)
    if legend:
        if legend_loc:
            ax.legend(frameon=False, loc=legend_loc, bbox_to_anchor=(0.93, 1.02))
        elif lowpass:
            ax.legend(frameon=False, loc="lower right", bbox_to_anchor=(1.0, 0.03))
        else:
            ax.legend(frameon=False, loc="upper right", bbox_to_anchor=(0.92, 1.0))
    if title:
        ax.set_title(title, loc="left")


def fig_spectrum(case: Path, out: Path, figsize=(7.2, 4.2)) -> Path:
    res = json.loads((Path(case) / "modes.json").read_text())
    fig, ax = plt.subplots(figsize=figsize)
    t = res["design"]["coil_type"].replace("pass", "-pass")
    _spectrum_axes(ax, res, f"{res['N']}-rung {t} birdcage: mode spectrum")
    ax.set_ylabel("frequency (MHz)")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def fig_spectrum_compare(case_a: Path, case_b: Path, out: Path, figsize=(12.5, 4.8), vertical: bool = False) -> Path:
    """(a) first-guess capacitor from the ladder model, (b) after retuning in Palace; same mesh."""
    ra = json.loads((Path(case_a) / "modes.json").read_text())
    rb = json.loads((Path(case_b) / "modes.json").read_text())
    fig, axs = plt.subplots(2, 1, figsize=figsize, sharex=True) if vertical else plt.subplots(1, 2, figsize=figsize, sharey=True)
    ca, cb = ra["design"]["c_leg"] * 1e12, rb["design"]["c_leg"] * 1e12
    _spectrum_axes(axs[0], ra, f"(a) ladder first guess: C = {ca:.2f} pF per rung", legend=not vertical)
    _spectrum_axes(axs[1], rb, f"(b) retuned in Palace: C = {cb:.2f} pF per rung", legend=vertical,
                   label_target=False, legend_loc="upper right" if vertical else None)
    axs[0].set_ylabel("frequency (MHz)")
    if vertical:
        axs[1].set_ylabel("frequency (MHz)")
        axs[0].set_xlabel("")
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def fig_currents(case: Path, out: Path) -> Path:
    res = json.loads((Path(case) / "modes.json").read_text())
    sel = pick_modes(res)
    N = res["N"]
    fig, axs = plt.subplots(1, len(sel), figsize=(2.9 * len(sel), 2.7), sharey=True)
    axs = np.atleast_1d(axs)
    n = np.arange(N)
    for ax, m in zip(axs, sel):
        c = np.array(m["currents_re"]) + 1j * np.array(m["currents_im"])
        c = np.real(c) / np.max(np.abs(c))
        col = ORANGE if m["homogeneous"] else BLUE
        ax.bar(n, c, width=0.7, color=col, edgecolor="white", lw=0.5)
        k = m["k"]
        nn = np.linspace(0, N - 1, 300) if k <= N // 4 else n.astype(float)
        # least-squares fit of cos(2 pi k n / N + phi)
        A = np.stack([np.cos(2 * np.pi * k * n / N), np.sin(2 * np.pi * k * n / N)], 1)
        cf, *_ = np.linalg.lstsq(A, c, rcond=None)
        ax.plot(nn, cf[0] * np.cos(2 * np.pi * k * nn / N) + cf[1] * np.sin(2 * np.pi * k * nn / N), color=INK,
                lw=1 if k <= N // 4 else 0, marker=None if k <= N // 4 else "_", ms=9, mew=1.4)
        ax.axhline(0, color=MUTED, lw=0.6)
        ax.set_title(f"{m['label']}  {m['f_mhz']:.1f} MHz", fontsize=9, color=col if m["homogeneous"] else INK)
        ax.set_xlabel("rung n", fontsize=8)
        ax.set_xticks([0, N // 4, N // 2, 3 * N // 4, N - 1])
        ax.tick_params(labelsize=7)
        ax.text(0.02, 0.03, f"k-purity {m['k_purity']:.3f}", transform=ax.transAxes, fontsize=7, color=INK2)
    axs[0].set_ylabel("rung current (norm.)", fontsize=8)
    src = sel[0].get("pattern_source", "")
    fig.suptitle(f"Current pattern per mode (Palace, {src}); line = cos(2πkn/N + φ) fit", fontsize=9, color=INK2, y=1.02)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def render_coil(case: Path, out: Path, show_mesh: bool = False, window=(1600, 1200)) -> Path:
    import pyvista as pv
    case = Path(case)
    s = json.loads((case / "surfaces.json").read_text())
    des = json.loads((case / "design.json").read_text())
    pl = pv.Plotter(off_screen=True, window_size=window)
    pl.set_background("white")

    def mesh(lbl):
        d = s[lbl]
        if not d["nodes"]:
            return None
        tri = np.array(d["tris"]).reshape(-1, 3)
        return pv.PolyData(np.array(d["nodes"]) * 100, np.hstack([np.full((len(tri), 1), 3), tri]).ravel())
    coil, caps, outer = mesh("coil"), mesh("caps"), mesh("outer")
    pl.add_mesh(coil, color="#d98a4a", specular=0.5, specular_power=20, ambient=0.25, smooth_shading=False, show_edges=show_mesh,
                edge_color="#5a3412", line_width=0.4)
    if caps is not None:
        pl.add_mesh(caps, color=BLUE, show_edges=show_mesh, edge_color="#0d2b5c", line_width=0.4)
    if des.get("shield_diameter"):  # back half of the shield, drawn analytically (clean edges)
        Rs = des["shield_diameter"] / 2 * 100
        Ls = (des.get("shield_length") or des["length"] + 0.20) * 100
        az = np.arctan2(-60.0, 95.0)
        th = np.linspace(az + np.pi / 2, az + 3 * np.pi / 2, 90)
        T, Z = np.meshgrid(th, np.array([-Ls / 2, Ls / 2]))
        sg = pv.StructuredGrid(Rs * np.cos(T), Rs * np.sin(T), Z)
        pl.add_mesh(sg, color="#b8bcc4", opacity=0.25, smooth_shading=True)
        for z in (-Ls / 2, Ls / 2):
            pl.add_mesh(pv.lines_from_points(np.c_[Rs * np.cos(th), Rs * np.sin(th), np.full_like(th, z)]),
                        color="#8a8f99", line_width=2)
    pl.camera_position = [(95, -60, 55), (0, 0, 0), (0, 0, 1)]
    pl.camera.zoom(1.25)
    pl.enable_anti_aliasing("ssaa")
    pl.screenshot(str(out))
    pl.close()
    return out


def render_mesh_cut(case: Path, out: Path, window=(1600, 1200)) -> Path:
    """Two views of the mesh: the y=0 slice through the volume mesh (grading towards the strips and the capacitor
    gaps) and a close-up of the surface mesh on the coil around one capacitor."""
    import gmsh
    import pyvista as pv
    from PIL import Image
    case = Path(case)
    info = json.loads((case / "build.json").read_text())
    des = json.loads((case / "design.json").read_text())
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 0)
    gmsh.open(str(case / info["mesh"]))
    tags, xyz, _ = gmsh.model.mesh.getNodes()
    xyz = xyz.reshape(-1, 3) * 100
    idx = np.zeros(int(tags.max()) + 1, int)
    idx[tags.astype(int)] = np.arange(len(tags))
    _, et = gmsh.model.mesh.getElementsByType(4)
    gmsh.finalize()
    tet = idx[et.astype(int)].reshape(-1, 4)
    ug = pv.UnstructuredGrid(np.hstack([np.full((len(tet), 1), 4), tet]).ravel(), np.full(len(tet), pv.CellType.TETRA), xyz)
    sl = ug.slice(normal=(0, 1, 0), origin=(0, 1e-3, 0))
    imgs = []
    pl = pv.Plotter(off_screen=True, window_size=window)
    pl.set_background("white")
    pl.add_mesh(sl, color="#eef2f8", show_edges=True, edge_color="#55606f", line_width=0.6)
    R, l = des["diameter"] / 2 * 100, des["length"] * 100
    for x in (-R, R):
        pl.add_mesh(pv.Line((x, -0.5, -l / 2), (x, -0.5, l / 2)), color="#d98a4a", line_width=5)
    pl.view_xz(negative=False)
    pl.camera.zoom(1.05)
    pl.screenshot(str(out.with_name(out.stem + "_slice.png")))
    pl.close()
    s = json.loads((case / "surfaces.json").read_text())
    pl = pv.Plotter(off_screen=True, window_size=window)
    pl.set_background("white")
    for lbl, col in (("coil", "#d98a4a"), ("caps", BLUE)):
        d = s[lbl]
        if d["nodes"]:
            tri = np.array(d["tris"]).reshape(-1, 3)
            pm = pv.PolyData(np.array(d["nodes"]) * 100, np.hstack([np.full((len(tri), 1), 3), tri]).ravel())
            pl.add_mesh(pm, color=col, show_edges=True, edge_color="#3d2410" if lbl == "coil" else "#0b2547",
                        line_width=1.0, ambient=0.3)
    pl.camera_position = [(R + 22, -9, 3), (R, 0, 0), (0, 0, 1)]
    pl.camera.view_angle = 30
    pl.screenshot(str(out.with_name(out.stem + "_closeup.png")))
    pl.close()
    a_ = Image.open(out.with_name(out.stem + "_slice.png")).convert("RGB")
    b_ = Image.open(out.with_name(out.stem + "_closeup.png")).convert("RGB")
    W = Image.new("RGB", (a_.width + b_.width, max(a_.height, b_.height)), "white")
    W.paste(a_, (0, 0))
    W.paste(b_, (a_.width, 0))
    W.save(out)
    return out


def all_figures(case: Path, figdir: Path) -> list[Path]:
    figdir = Path(figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    outs = [fig_spectrum(case, figdir / "spectrum.png"), fig_b1_maps(case, figdir / "b1_maps.png"),
            fig_currents(case, figdir / "currents.png")]
    try:
        outs.append(render_coil(case, figdir / "coil.png"))
        outs.append(render_mesh_cut(case, figdir / "mesh.png"))
    except Exception as e:  # pyvista / OpenGL not available: figures above are still written
        print("3D render skipped:", e)
    print("figures:", ", ".join(p.name for p in outs))
    return outs


def fig_loaded_compare(case_u: Path, case_l: Path, out: Path) -> Path:
    """Quadrature |B1+| of the k=1 pair, unloaded vs phantom-loaded, normalised to the value at the coil centre."""
    rows = []
    for case in (case_u, case_l):
        res, raw, pts, sl = _ctx(case)
        hp = res["homogeneous_pair"]
        q = quadrature(raw["B"][hp["m"][0] - 1], raw["B"][hp["m"][1] - 1], sl["roi"])
        v = np.abs(q["b1p"])
        ax_pts = pts[sl["axial"]]
        c = v[sl["axial"]][np.argmin((ax_pts[:, :2] ** 2).sum(1))]
        q_val = np.mean([res["modes"][m - 1]["Q"] for m in hp["m"]])
        rows.append((res, pts, sl, v / c, hp, q_val))
    fig, axs = plt.subplots(2, 2, figsize=(7.6, 8.6), gridspec_kw={"width_ratios": [1, 0.8]})
    for i, (res, pts, sl, v, hp, qv) in enumerate(rows):
        des = res["design"]
        for j, (plane, cut) in enumerate((("axial", (0, 1)), ("coronal", (0, 2)))):
            ax = axs[i, j]
            img, ext = _grid(pts[sl[plane]][:, cut], v[sl[plane]])
            im = ax.imshow(img, origin="lower", extent=[e * 100 for e in ext], cmap=SEQ, vmin=0, vmax=1.6,
                           interpolation="bilinear")
            _coil_marks(ax, des, plane)
            ph = des.get("phantom")
            if ph:
                if plane == "axial":
                    t = np.linspace(0, 2 * np.pi, 200)
                    ax.plot(ph["radius"] * 100 * np.cos(t), ph["radius"] * 100 * np.sin(t), color=AQUA, lw=1.4)
                else:
                    ax.add_patch(plt.Rectangle((-ph["radius"] * 100, -ph["length"] * 50), ph["radius"] * 200,
                                               ph["length"] * 100, fill=False, ec=AQUA, lw=1.4))
            ax.set_aspect("equal")
            ax.tick_params(labelsize=7)
            for sp in ax.spines.values():
                sp.set_visible(False)
            if j == 0:
                lab = (f"saline phantom (εr {ph['eps_r']:g}, σ {ph['sigma']:g} S/m)" if ph else "unloaded")
                qs = f", Q = {qv:.0f}" if np.isfinite(qv) and qv < 1e6 else ""
                ax.set_title(f"{lab}\nk = 1 at {hp['f_mean_mhz']:.2f} MHz{qs}; B1+ ROI std/mean {hp['b1plus_roi_cv'] * 100:.0f} %",
                             fontsize=9, loc="left", color=INK)
    cb = fig.colorbar(im, ax=axs, shrink=0.6, pad=0.02)
    cb.set_label("quadrature |B1+| / value at the centre", fontsize=8)
    fig.suptitle("k = 1 pair driven in quadrature: axial plane z = 0 (left), coronal plane y = 0 (right), axes in cm",
                 fontsize=9.5, color=INK2, x=0.45)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out
