"""Symmetry-plane study: figures in docs/img and results/lp16_symmetry.json.

    python tools/symmetry_report.py --full runs/lp16_final --quarter runs/sym_q --eighth runs/sym_e_em \
        [--through runs/sym_qt_em --wrong-c runs/sym_qt_em_wrongC] [--conv runs/a runs/b ...] [--also <dir>]

`--quarter runs/sym_q` expects the four PEC/PMC combinations as runs/sym_q_em, _me, _mm, _ee
(first letter: plane x = 0, second: plane y = 0; e = PEC, m = PMC), all on one mesh.
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from birdcage_palace import plots  # noqa: E402
from birdcage_palace.modes import load, quadrature  # noqa: E402

PEC_C, PMC_C = "#c2410c", "#0e9f6e"   # plane colours in the renders and figures (PEC burnt orange, PMC green)
COMBOS = (("em", "PEC", "PMC"), ("me", "PMC", "PEC"), ("mm", "PMC", "PMC"), ("ee", "PEC", "PEC"))


# ------------------------------------------------------------------------------------------------ data helpers
def ctx(case: Path):
    case = Path(case)
    res = json.loads((case / "modes.json").read_text())
    pr = json.loads((case / "probes.json").read_text())
    raw = load(case)
    pts = np.array(pr["points"])
    na, nc = pr["axial"], pr["coronal"]
    return res, raw, pts, {"axial": slice(0, na), "coronal": slice(na, na + nc), "roi": slice(na + nc, None)}


def cost(case: Path) -> dict:
    case = Path(case)
    log = (case / "palace.log").read_text(errors="replace")
    st = json.loads((case / "solve_stats.json").read_text())
    mem = re.findall(r"peak per-node memory usage is:.*Total ([\d.]+)([GM])", log)
    gb = (float(mem[-1][0]) / (1024 if mem[-1][1] == "M" else 1)) if mem else None
    tot = re.findall(r"^Total\s+([\d.]+)", log, re.M)
    cfg = json.loads((case / "config.json").read_text())
    info = json.loads((case / "build.json").read_text())
    return {"unknowns": st["unknowns"], "tets": info["tets"], "palace_gb": gb, "palace_s": float(tot[-1]) if tot else None,
            "shared_machine": bool(st.get("shared_machine", False)),
            "wall_s": st["wall_s"], "ranks": st["ranks"], "order": cfg["Solver"]["Order"],
            "modes": cfg["Solver"]["Eigenmode"]["N"]}


def real_pair(raw, res):
    """The two k=1 eigenvectors of the full model, real, and rotated within the pair so mode a has B1 along y."""
    hp = res["homogeneous_pair"]["m"]
    sl_roi = slice(res["_roi0"], None)
    vs = []
    for m in hp:
        B = raw["B"][m - 1]
        v = B[sl_roi].sum(0)
        B = B * np.exp(-1j * np.angle(v[np.argmax(np.abs(v))]))
        vs.append(np.real(B))
    a, b = vs
    ma, mb = a[sl_roi].mean(0), b[sl_roi].mean(0)
    # combination c = cos t a + sin t b with zero mean Bx in the ROI
    t = np.arctan2(-ma[0], mb[0])
    y = np.cos(t) * a + np.sin(t) * b
    x = -np.sin(t) * a + np.cos(t) * b
    if y[sl_roi].mean(0)[1] < 0:
        y = -y
    if x[sl_roi].mean(0)[0] < 0:
        x = -x
    return y, x, t


def parity_check(full: Path) -> dict:
    """Mirror test of the full-model field: is the plane x = 0 / y = 0 PEC-like or PMC-like for each k=1 member?"""
    res, raw, pts, sl = ctx(full)
    res["_roi0"] = sl["roi"].start
    By, Bx, t = real_pair(raw, res)
    roi = pts[sl["roi"]] - 1.1e-5            # undo the probe nudge: the ROI grid is symmetric about x = 0 and y = 0
    key = {tuple(np.round(p, 5)): i for i, p in enumerate(roi)}
    out = {}
    for name, B in (("B1 along y", By[sl["roi"]]), ("B1 along x", Bx[sl["roi"]])):
        rec = {}
        for ax_i, ax in ((0, "x"), (1, "y")):
            M = np.ones(3)
            M[ax_i] = -1
            src, dst = [], []
            for i, p in enumerate(roi):
                q = p.copy()
                q[ax_i] = -q[ax_i]
                j = key.get(tuple(np.round(q, 5)))
                if j is not None and abs(p[ax_i]) > 1e-4:
                    src.append(i)
                    dst.append(j)
            Bs, Bd = B[src], B[dst]
            nrm = np.linalg.norm(Bs)
            pec = np.linalg.norm(Bd - Bs * M) / nrm        # PEC plane: B(Mr) = +M B(r)
            pmc = np.linalg.norm(Bd + Bs * M) / nrm        # PMC plane: B(Mr) = -M B(r)
            rec[f"{ax}=0"] = {"residual_if_pec": float(pec), "residual_if_pmc": float(pmc),
                              "verdict": "PEC" if pec < pmc else "PMC", "pairs": len(src)}
        out[name] = rec
    # rung currents of the y-polarised member, with the same phase per eigenvector as its field
    hp = res["homogeneous_pair"]["m"]
    Is = []
    for m in hp:
        v = raw["B"][m - 1][sl["roi"]].sum(0)
        Is.append(np.real(raw["I"][m - 1] * np.exp(-1j * np.angle(v[np.argmax(np.abs(v))]))))
    N = res["N"]
    phi = 2 * np.pi * np.arange(N) / N
    Iy = np.cos(t) * Is[0] + np.sin(t) * Is[1]
    A = np.stack([np.cos(phi), np.sin(phi)], 1)
    cf, *_ = np.linalg.lstsq(A, Iy, rcond=None)
    fit = A @ cf
    # 180 deg = I_n ~ -cos(phi_n) for B1 along +y: +z current at phi = 0 gives B along -y (right-hand rule)
    out["B1y_member_rung_currents"] = {"phase_deg_of_cos_fit": float(np.degrees(np.arctan2(cf[1], cf[0]))),
                                       "fit_residual": float(np.linalg.norm(Iy - fit) / np.linalg.norm(Iy))}
    return out


# ------------------------------------------------------------------------------------------------ figures
def render_model(case: Path, out: Path, mode: str = "quarter", show_mesh: bool = False, window=(1500, 1500),
                 camera=None, ghost: bool = True) -> Path:
    """3D render of a symmetry case. mode = 'quarter' (the modelled part, planes coloured), 'full' (the cage rebuilt
    by mirroring, no planes), 'cut' (full cage with the planes drawn through it)."""
    import pyvista as pv
    case = Path(case)
    s = json.loads((case / "surfaces.json").read_text())
    des = json.loads((case / "design.json").read_text())
    planes = [a for a in "xyz" if f"sym_{a}" in s]

    def poly(lbl, mirror=(1, 1, 1)):
        d = s[lbl]
        tri = np.array(d["tris"]).reshape(-1, 3)
        P = np.array(d["nodes"]) * 100 * np.array(mirror)
        return pv.PolyData(P, np.hstack([np.full((len(tri), 1), 3), tri]).ravel())

    mirrors = [(1, 1, 1)]
    if mode in ("full", "cut"):
        mirrors = [(sx, sy, sz) for sx in ((1, -1) if "x" in planes else (1,)) for sy in ((1, -1) if "y" in planes else (1,))
                   for sz in ((1, -1) if "z" in planes else (1,))]
    pl = pv.Plotter(off_screen=True, window_size=window)
    pl.set_background("white")
    for mi in mirrors:
        pl.add_mesh(poly("coil", mi), color="#d98a4a", specular=0.5, specular_power=20, ambient=0.25,
                    show_edges=show_mesh, edge_color="#5a3412", line_width=0.5)
        pl.add_mesh(poly("caps", mi), color=plots.BLUE, show_edges=show_mesh, edge_color="#0d2b5c", line_width=0.5)
    Rs = des["shield_diameter"] / 2 * 100
    Ls = (des.get("shield_length") or des["length"] + 0.20) * 100
    if mode == "quarter":
        if ghost:  # the three mirrored quarters as a faint ghost
            for mi in [(sx, sy, 1) for sx in (1, -1) for sy in (1, -1)][1:]:
                pl.add_mesh(poly("coil", mi), color="#d9d4c8", opacity=0.35)
        bc = des.get("sym_bc", {})
        for a in planes:
            col = PEC_C if bc.get(a) == "pec" else PMC_C
            pl.add_mesh(poly(f"sym_{a}"), color=col, opacity=0.45, show_edges=show_mesh, edge_color="#333333",
                        line_width=0.4, ambient=0.75, diffuse=0.3, specular=0.0)
        # quarter shield outline
        th = np.linspace(0, np.pi / 2, 40)
        for z in (-Ls / 2, Ls / 2):
            pl.add_mesh(pv.lines_from_points(np.c_[Rs * np.cos(th), Rs * np.sin(th), np.full_like(th, z)]),
                        color="#8a8f99", line_width=2)
        for x, y in ((Rs, 0), (0, Rs)):
            pl.add_mesh(pv.Line((x, y, -Ls / 2), (x, y, Ls / 2)), color="#8a8f99", line_width=2)
    else:
        az = np.arctan2(-60.0, 95.0)
        th = np.linspace(az + np.pi / 2, az + 3 * np.pi / 2, 90)
        T, Z = np.meshgrid(th, np.array([-Ls / 2, Ls / 2]))
        pl.add_mesh(pv.StructuredGrid(Rs * np.cos(T), Rs * np.sin(T), Z), color="#b8bcc4", opacity=0.25)
        for z in (-Ls / 2, Ls / 2):
            pl.add_mesh(pv.lines_from_points(np.c_[Rs * np.cos(th), Rs * np.sin(th), np.full_like(th, z)]),
                        color="#8a8f99", line_width=2)
        if mode == "cut":
            bc = des.get("sym_bc", {})
            for a in [p for p in planes if p in "xy"]:
                if a == "x":
                    P = pv.Plane(center=(0, 0, 0), direction=(1, 0, 0), i_size=Ls, j_size=2 * Rs)
                else:
                    P = pv.Plane(center=(0, 0, 0), direction=(0, 1, 0), i_size=Ls, j_size=2 * Rs)
                pl.add_mesh(P, color=PEC_C if bc.get(a) == "pec" else PMC_C, opacity=0.40, ambient=0.75, diffuse=0.3)
    pl.camera_position = camera or ([(95, -60, 55), (0, 0, 0), (0, 0, 1)] if mode != "quarter" else
                                    [(125, 100, 75), (4, 4, -2), (0, 0, 1)])
    pl.camera.zoom(1.25 if mode != "quarter" else 1.15)
    pl.enable_anti_aliasing("ssaa")
    pl.screenshot(str(out))
    pl.close()
    return out


def fig_quarter_render(case: Path, out: Path, work: Path) -> Path:
    from PIL import Image
    a = render_model(case, work / "render_quarter.png", "quarter")
    b = render_model(case, work / "render_quarter_mesh.png", "quarter", show_mesh=True, ghost=False)
    ims = [Image.open(p).convert("RGB") for p in (a, b)]
    fig, axs = plt.subplots(1, 2, figsize=(12, 6.4))
    des = json.loads((Path(case) / "design.json").read_text())
    bc = des["sym_bc"]
    for ax, im, t in zip(axs, ims, ("quarter model (ghost: the three mirrored quarters)", "its surface mesh")):
        arr = np.asarray(im)
        ys, xs = np.where((arr < 245).any(2))
        ax.imshow(im.crop((xs.min() - 10, ys.min() - 10, xs.max() + 10, ys.max() + 10)))
        ax.axis("off")
        ax.set_title(t, fontsize=11, color=plots.INK2)
    from matplotlib.patches import Patch
    hs = [Patch(color=PEC_C if bc["x"] == "pec" else PMC_C, alpha=0.6, label=f"x = 0 plane: {bc['x'].upper()}"),
          Patch(color=PEC_C if bc["y"] == "pec" else PMC_C, alpha=0.6, label=f"y = 0 plane: {bc['y'].upper()}"),
          Patch(color="#d98a4a", label="rungs and end rings (PEC strips)"), Patch(color=plots.BLUE, label="lumped capacitors")]
    fig.legend(handles=hs, loc="lower center", ncol=4, frameon=False, fontsize=10)
    fig.suptitle("16-rung low-pass birdcage, quarter model: symmetry planes between rungs (cage rotated by π/16)",
                 fontsize=11.5, color=plots.INK, y=0.97)
    fig.savefig(out, bbox_inches="tight", dpi=170)
    plt.close(fig)
    return out


def fig_costs(rows: list[tuple[str, dict]], out: Path, figsize=(11, 3.6), big: bool = False, vertical: bool = False) -> Path:
    """Bars: unknowns, Palace peak memory, Palace time for each model."""
    fs = 1.45 if big else 1.0
    fig, axs = plt.subplots(3, 1, figsize=figsize) if vertical else plt.subplots(1, 3, figsize=figsize)
    labs = [r[0] for r in rows]
    cols = [plots.INK2] + [plots.BLUE] * (len(rows) - 1)
    for ax, (key, unit, sc, fmt) in zip(axs, (("unknowns", "unknowns (M)", 1e-6, "{:.2f} M"),
                                              ("palace_gb", "peak memory (GB, Palace estimate)", 1, "{:.1f} GB"),
                                              ("palace_s", "solve time (min)", 1 / 60, "{:.1f} min"))):
        v = [r[1][key] * sc for r in rows]
        y = np.arange(len(rows))[::-1]
        shared = [key == "palace_s" and r[1].get("shared_machine") for r in rows]
        for yi, vi, ci, sh in zip(y, v, cols, shared):
            ax.barh(yi, vi, color="white" if sh else ci, edgecolor=ci, hatch="///" if sh else None, height=0.62, lw=1)
            ax.text(vi + max(v) * 0.02, yi, fmt.format(vi) + (" (CPU shared with\nanother job)" if sh else ""),
                    va="center", fontsize=(7.5 if sh else 9) * fs, color=plots.INK)
        ax.set_yticks(y)
        ax.set_yticklabels(labs if (ax is axs[0] or vertical) else [""] * len(rows), fontsize=9 * fs)
        ax.set_xlim(0, max(v) * 1.35)
        ax.set_title(unit, fontsize=10 * fs, loc="left")
        ax.tick_params(axis="x", labelsize=8 * fs)
        ax.grid(axis="x", color=plots.GRID, lw=0.6)
        ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight", dpi=170)
    plt.close(fig)
    return out


def _panel(ax, pts, v, plane, des, title, quarter_lines=False, vmax=None):
    cut = (0, 1) if plane == "axial" else (0, 2)
    img, ext = plots._grid(pts[:, cut], v)
    ext = [e * 100 for e in ext]
    im = ax.imshow(img, origin="lower", extent=ext, cmap=plots.SEQ, vmin=0, vmax=vmax or 1, interpolation="bilinear")
    plots._coil_marks(ax, des, plane)
    if quarter_lines:
        bc = des["sym_bc"]
        if plane == "axial":
            ax.plot([0, 0], [ext[2], ext[3]], color=PEC_C if bc["x"] == "pec" else PMC_C, lw=1.6)
            ax.plot([ext[0], ext[1]], [0, 0], color=PEC_C if bc["y"] == "pec" else PMC_C, lw=1.6)
            ax.text(ext[1] * 0.42, ext[3] * 0.42, "modelled", ha="center", va="center", fontsize=8, color="white")
        else:
            ax.plot([0, 0], [ext[2], ext[3]], color=PEC_C if bc["x"] == "pec" else PMC_C, lw=1.6)
    ax.set_aspect("equal")
    ax.tick_params(labelsize=7)
    for s in ax.spines.values():
        s.set_visible(False)
    if title:
        ax.set_title(title, fontsize=9.5)
    return im


def fig_b1_compare(full: Path, q_em: Path, q_me: Path, out: Path, figsize=(13.0, 7.6), big=False,
                   linear_only: bool = False) -> tuple[Path, dict]:
    """Linear k=1 member with B1 along y, full model vs mirrored quarter, and the quadrature |B1+| of each."""
    res, raw, pts, sl = ctx(full)
    res["_roi0"] = sl["roi"].start
    By, Bx, _ = real_pair(raw, res)
    rq, rawq, ptsq, slq = ctx(q_em)
    rq2, rawq2, _, _ = ctx(q_me)
    mq = [m for m in rq["modes"] if m["homogeneous"]][0]
    mq2 = [m for m in rq2["modes"] if m["homogeneous"]][0]
    Bq = rawq["B"][mq["m"] - 1]
    v = Bq[slq["roi"]].sum(0)
    Bq = np.real(Bq * np.exp(-1j * np.angle(v[np.argmax(np.abs(v))])))
    hp = res["homogeneous_pair"]["m"]
    qf = quadrature(raw["B"][hp[0] - 1], raw["B"][hp[1] - 1], sl["roi"])
    qq = quadrature(rawq["B"][mq["m"] - 1], rawq2["B"][mq2["m"] - 1], slq["roi"])
    stats = {}

    def norm_c(v, pts_, sl_):  # normalise to the value at the coil centre
        ax_pts = pts_[sl_["axial"]]
        return v / v[sl_["axial"]][np.argmin((ax_pts[:, :2] ** 2).sum(1))]
    cols = [("full model\nlinear k=1, |B$_t$|", plots._bt(By), res["design"], pts, sl, False),
            ("quarter, mirrored\nlinear k=1, |B$_t$|", plots._bt(Bq), rq["design"], ptsq, slq, True),
            ("full model\nquadrature |B$_1^+$|", np.abs(qf["b1p"]), res["design"], pts, sl, False),
            ("two quarter runs, mirrored\nquadrature |B$_1^+$|", np.abs(qq["b1p"]), rq["design"], ptsq, slq, True)]
    if linear_only:
        cols = cols[:2]
    fig, axs = plt.subplots(2, len(cols), figsize=figsize, gridspec_kw={"height_ratios": [1, 1.25]})
    for j, (title, vv, des, pts_, sl_, ql) in enumerate(cols):
        vv = norm_c(vv, pts_, sl_)
        roi = vv[sl_["roi"]]
        stats[["full_linear", "quarter_linear", "full_quad", "quarter_quad"][j]] = float(roi.std() / roi.mean())
        for i, plane in enumerate(("axial", "coronal")):
            im = _panel(axs[i, j], pts_[sl_[plane]], vv[sl_[plane]], plane, des, title if i == 0 else None, ql, vmax=1.6)
            axs[i, j].set_xlabel("x (cm)", fontsize=8)
            if j == 0:
                axs[i, j].set_ylabel("y (cm)" if plane == "axial" else "z (cm)", fontsize=8)
        axs[1, j].text(0.5, -0.24, f"ROI std/mean {stats[list(stats)[-1]] * 100:.2f} %", transform=axs[1, j].transAxes,
                       ha="center", fontsize=9 if not big else 11, color=plots.INK)
    cb = fig.colorbar(im, ax=axs, shrink=0.55, pad=0.015)
    cb.set_label("field / value at the coil centre", fontsize=8)
    if linear_only:
        fig.suptitle("one k = 1 member (B1 along y); axial plane z = 0 (top), coronal plane y = 0 (bottom)",
                     fontsize=12, color=plots.INK2, x=0.45)
    else:
        fig.suptitle("k = 1 mode, axial plane z = 0 (top) and coronal plane y = 0 (bottom). Linear: the member with B1 along y. "
                 "Lines: x = 0 PEC (orange), y = 0 PMC (green); dashed = ROI.\nThe symmetric models rotate the cage by "
                 "π/16, so their coronal plane passes between rungs instead of through two (hence the weaker field at x = ±13 cm).",
                 fontsize=10, color=plots.INK2, x=0.45)
    fig.savefig(out, bbox_inches="tight", dpi=170)
    plt.close(fig)
    return out, stats


def mode_family(m: dict, combo: str) -> str:
    if m.get("cap_free"):
        return "no capacitor current"
    k = m.get("k")
    return f"k={k}"


def fig_spectrum_combos(full: Path, qbase: str, out: Path, fmax: float, figsize=(11.5, 5.0), big=False) -> tuple[Path, dict]:
    rf = json.loads((Path(full) / "modes.json").read_text())
    rows = [("full model\n(20 modes)", [m for m in rf["modes"]], plots.INK2)]
    union = []
    table = {}
    for c, bx, by in COMBOS:
        r = json.loads((Path(f"{qbase}_{c}") / "modes.json").read_text())
        ms = [m for m in r["modes"] if m["f_mhz"] <= fmax]
        rows.append((f"x=0 {bx}, y=0 {by}", ms, plots.BLUE))
        union += [(m["f_mhz"], c, m) for m in ms]
        table[c] = [{"f_mhz": round(m["f_mhz"], 3), "label": m["label"], "k": m.get("k")} for m in r["modes"]]
    fs = 1.35 if big else 1.0
    fig, ax = plt.subplots(figsize=figsize)
    for i, (lab, ms, col) in enumerate(rows):
        y = len(rows) - 1 - i
        for m in ms:
            cf = m.get("cap_free")
            mk = "x" if cf else "o"
            hom = m.get("homogeneous")
            ax.plot(m["f_mhz"], y, mk, ms=(9 if hom else 7) * fs, color=plots.ORANGE if hom else (plots.MUTED if cf else col),
                    mec="white" if not cf else None, mew=0.8 if not cf else 1.6 * fs, zorder=3)
            if not cf:
                dx = {7: -2.2, 8: 2.2}.get(m["k"], 0) if i == 0 else 0
                ax.text(m["f_mhz"] + dx, y + 0.22, str(m["k"]), ha="center", va="bottom", fontsize=7.5 * fs, color=plots.INK2)
        ax.text(fmax * 1.01 + 8, y, "", va="center")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows][::-1], fontsize=9 * fs)
    ax.set_xlim(100, fmax + 10)
    ax.set_ylim(-0.6, len(rows) - 0.3)
    ax.axhline(len(rows) - 1.5, color=plots.GRID, lw=1)
    ax.set_xlabel("eigenfrequency (MHz)", fontsize=10 * fs)
    ax.tick_params(axis="x", labelsize=9 * fs)
    ax.grid(axis="x", color=plots.GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.plot([], [], "o", color=plots.ORANGE, label="homogeneous k=1")
    ax.plot([], [], "o", color=plots.BLUE, label="ladder mode (number = k)")
    ax.plot([], [], "x", color=plots.MUTED, mew=1.6, label="no capacitor current (ring / stray-C)")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.13), ncol=3, fontsize=8.5 * fs)
    fig.savefig(out, bbox_inches="tight", dpi=170)
    plt.close(fig)
    return out, table


def match_union(full: Path, qbase: str, fmax: float) -> dict:
    rf = json.loads((Path(full) / "modes.json").read_text())
    pool = []
    for c, _, _ in COMBOS:
        r = json.loads((Path(f"{qbase}_{c}") / "modes.json").read_text())
        pool += [dict(m, combo=c) for m in r["modes"]]
    used, rows = set(), []
    for m in rf["modes"]:
        best, bd = None, 1e9
        for j, q in enumerate(pool):
            if j in used or q.get("cap_free") != m.get("cap_free") or (not m.get("cap_free") and q.get("k") != m.get("k")):
                continue
            d = abs(q["f_mhz"] - m["f_mhz"])
            if d < bd:
                best, bd = j, d
        if best is not None:
            used.add(best)
            q = pool[best]
            rows.append({"full_f_mhz": round(m["f_mhz"], 3), "label": m["label"], "combo": q["combo"],
                         "quarter_f_mhz": round(q["f_mhz"], 3), "rel_diff": (q["f_mhz"] - m["f_mhz"]) / m["f_mhz"]})
        else:
            rows.append({"full_f_mhz": round(m["f_mhz"], 3), "label": m["label"], "combo": None})
    extra = [{"combo": q["combo"], "f_mhz": round(q["f_mhz"], 3), "label": q["label"]} for j, q in enumerate(pool)
             if j not in used and q["f_mhz"] <= fmax]
    return {"matched": rows, "unmatched_quarter_modes_below_fmax": extra}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", required=True)
    ap.add_argument("--quarter", required=True, help="prefix of the four combination cases, e.g. runs/sym_q")
    ap.add_argument("--eighth", default=None)
    ap.add_argument("--through", default=None)
    ap.add_argument("--wrong-c", default=None)
    ap.add_argument("--conv", nargs="*", default=[])
    ap.add_argument("--also", action="append", default=[])
    ap.add_argument("--work", default="runs/sym_figs")
    a = ap.parse_args()
    full, qb = Path(a.full), a.quarter
    work = Path(a.work)
    work.mkdir(parents=True, exist_ok=True)
    img = ROOT / "docs" / "img"
    rf = json.loads((full / "modes.json").read_text())
    fmax = max(m["f_mhz"] for m in rf["modes"]) + 5
    out = {"full": {"case": str(full), "cost": cost(full), "k1_mhz": rf["homogeneous_pair"]["f_mhz"],
                    "k1_mean_mhz": rf["homogeneous_pair"]["f_mean_mhz"], "roi_bt_cv": rf["homogeneous_pair"]["roi_bt_cv"],
                    "b1plus_roi_cv": rf["homogeneous_pair"]["b1plus_roi_cv"]}}
    out["parity_check_full_model"] = parity_check(full)
    out["combos"] = {}
    for c, bx, by in COMBOS:
        case = Path(f"{qb}_{c}")
        r = json.loads((case / "modes.json").read_text())
        hp = r.get("homogeneous_pair")
        out["combos"][c] = {"x0": bx, "y0": by, "cost": cost(case),
                            "k1_mhz": hp["f_mean_mhz"] if hp else None, "roi_bt_cv": hp["roi_bt_cv"][0] if hp else None,
                            "modes": [{"f_mhz": round(m["f_mhz"], 3), "label": m["label"]} for m in r["modes"]]}
    out["union_vs_full"] = match_union(full, qb, fmax)
    for key, case in (("eighth", a.eighth), ("through_rung", a.through), ("through_rung_wrong_C", a.wrong_c)):
        if case:
            r = json.loads((Path(case) / "modes.json").read_text())
            hp = r["homogeneous_pair"]
            out[key] = {"case": case, "cost": cost(Path(case)), "k1_mhz": hp["f_mean_mhz"], "roi_bt_cv": hp["roi_bt_cv"][0],
                        "modes": [{"f_mhz": round(m["f_mhz"], 3), "label": m["label"]} for m in r["modes"]]}
    conv = []
    for case in a.conv:
        case = Path(case)
        r = json.loads((case / "modes.json").read_text())
        d = json.loads((case / "design.json").read_text())
        hp = r["homogeneous_pair"]
        conv.append({"case": str(case), "h_mm": [d["h_coil"] * 1e3, d["h_gap"] * 1e3, d["h_max"] * 1e3], "order": d["order"],
                     "symmetry": d["symmetry"], "cost": cost(case), "k1_mhz": hp["f_mean_mhz"], "roi_bt_cv": hp["roi_bt_cv"][0]})
    out["convergence"] = conv
    out["quarter_prefix"] = qb
    o2 = sorted([c for c in conv if c["order"] == 2], key=lambda c: -c["h_mm"][0])
    if len(o2) >= 3:
        # Richardson on three h-levels at order 2 with the observed rate (h ratios need not be equal: solve for p)
        from scipy.optimize import brentq
        (h1, f1), (h2, f2), (h3, f3) = [(c["h_mm"][0], c["k1_mhz"]) for c in o2[-3:]]
        g = lambda p: (f1 - f2) / (f2 - f3) - (h1 ** p - h2 ** p) / (h2 ** p - h3 ** p)
        try:
            p_obs = brentq(g, 0.3, 8)
            finf = f3 - (f2 - f3) / ((h2 / h3) ** p_obs - 1)
        except ValueError:
            p_obs, finf = None, None
        out["richardson"] = {"h_mm": [h1, h2, h3], "f_mhz": [f1, f2, f3], "observed_rate": p_obs, "f_extrapolated_mhz": finf}

    figs = []
    figs.append(fig_quarter_render(Path(f"{qb}_em"), img / "symmetry_quarter_model.png", work))
    rows = [("full model", out["full"]["cost"]), ("quarter (x=0 PEC, y=0 PMC)", out["combos"]["em"]["cost"])]
    if a.eighth:
        rows.append(("eighth (+ z=0 PEC)", out["eighth"]["cost"]))
    figs.append(fig_costs(rows, img / "symmetry_costs.png"))
    p, st = fig_b1_compare(full, Path(f"{qb}_em"), Path(f"{qb}_me"), img / "symmetry_b1_compare.png")
    out["b1_maps_roi_cv_centre_normalised"] = st
    figs.append(p)
    p, _ = fig_spectrum_combos(full, qb, img / "symmetry_spectrum.png", fmax)
    figs.append(p)
    (ROOT / "results" / "lp16_symmetry.json").write_text(json.dumps(out, indent=1))
    for d in a.also:
        d = Path(d)
        d.mkdir(parents=True, exist_ok=True)
        for f in figs:
            shutil.copy(f, d / f.name)
    print(json.dumps({k: out[k] for k in ("parity_check_full_model",)}, indent=1))
    print("wrote", ", ".join(f.name for f in figs), "+ results/lp16_symmetry.json")


if __name__ == "__main__":
    main()
