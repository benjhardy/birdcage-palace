"""Silent 4:5 (1080x1350) explainer clip from a solved case: coil -> spectrum -> B1 maps per mode -> the k=1 result.

    python tools/make_video.py --case runs/lp16_final --out docs/img/birdcage_modes.mp4 --ffmpeg <path to ffmpeg>

Look: cream graph paper, cobalt ink, plots kept in their original colours on white plates with a thin cobalt border.
Every number on screen is read from the case's modes.json / design.json.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from birdcage_palace import plots  # noqa: E402
from birdcage_palace.modes import load, quadrature  # noqa: E402

W, H, FPS = 1080, 1350, 30
CREAM, GRID1, GRID2 = (246, 242, 231), (232, 226, 210), (218, 210, 190)
COBALT, COBALT2 = (31, 63, 174), (88, 110, 190)
FONTS = Path("C:/Windows/Fonts")


def font(name, size):
    for n in name.split("|"):
        p = FONTS / n
        if p.exists():
            return ImageFont.truetype(str(p), size)
    return ImageFont.load_default()


SERIF = lambda s: font("cambria.ttc|georgia.ttf", s)  # noqa: E731
SERIF_I = lambda s: font("cambriai.ttf|georgiai.ttf", s)  # noqa: E731
SANS = lambda s: font("segoeui.ttf", s)  # noqa: E731
SANS_B = lambda s: font("segoeuib.ttf", s)  # noqa: E731
MONO = lambda s: font("consola.ttf", s)  # noqa: E731


def paper() -> Image.Image:
    im = Image.new("RGB", (W, H), CREAM)
    d = ImageDraw.Draw(im)
    for x in range(0, W, 27):
        d.line([(x, 0), (x, H)], fill=GRID2 if x % 135 == 0 else GRID1, width=1)
    for y in range(0, H, 27):
        d.line([(0, y), (W, y)], fill=GRID2 if y % 135 == 0 else GRID1, width=1)
    return im


def ease(t):  # power3.out
    t = min(max(t, 0.0), 1.0)
    return 1 - (1 - t) ** 3


def plate(path: Path, box_w: int, box_h: int, pad: int = 14) -> Image.Image:
    src = Image.open(path).convert("RGB")
    # crop white margins
    a = np.asarray(src)
    mask = (a < 245).any(2)
    ys, xs = np.where(mask)
    if len(xs):
        src = src.crop((max(xs.min() - 8, 0), max(ys.min() - 8, 0), min(xs.max() + 8, src.width), min(ys.max() + 8, src.height)))
    s = min((box_w - 2 * pad) / src.width, (box_h - 2 * pad) / src.height)
    src = src.resize((int(src.width * s), int(src.height * s)), Image.LANCZOS)
    pl = Image.new("RGB", (src.width + 2 * pad, src.height + 2 * pad), "white")
    pl.paste(src, (pad, pad))
    d = ImageDraw.Draw(pl)
    d.rectangle([0, 0, pl.width - 1, pl.height - 1], outline=COBALT, width=2)
    return pl


class Scene:
    def __init__(self, dur, items):
        self.dur, self.items = dur, items  # items: (t_in, kind, payload, xy)


def draw_text(d, xy, text, f, fill=COBALT, anchor="la"):
    d.text(xy, text, font=f, fill=fill, anchor=anchor)


def mode_panel(case, res, raw, pts, sl, m, out: Path, quad=None, partner=None):
    """Axial + coronal |B_t| maps and the rung-current bars for one mode, as one figure."""
    des = res["design"]
    N = res["N"]
    fig = plt.figure(figsize=(9.6, 8.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.45], width_ratios=[1, 0.82], hspace=0.30, wspace=0.16)
    if quad is not None:
        v = np.abs(quad["b1p"])
        lab = "|B$_1^+$| (quadrature)"
    else:
        B = raw["B"][m["m"] - 1]
        v = np.sqrt(np.abs(B[:, 0]) ** 2 + np.abs(B[:, 1]) ** 2)
        lab = "|B$_t$|"
    vmax = np.nanmax(np.r_[v[sl["axial"]], v[sl["coronal"]]])
    for j, (plane, cut) in enumerate((("axial", (0, 1)), ("coronal", (0, 2)))):
        ax = fig.add_subplot(gs[0, j])
        p = pts[sl[plane]][:, cut]
        img, ext = plots._grid(p, v[sl[plane]] / vmax)
        ext = [e * 100 for e in ext]
        ax.imshow(img, origin="lower", extent=ext, cmap=plots.SEQ, vmin=0, vmax=1, interpolation="bilinear")
        plots._coil_marks(ax, des, plane)
        ax.set_aspect("equal")
        ax.set_title("axial plane z = 0" if plane == "axial" else "coronal plane y = 0", fontsize=15, color=plots.INK2)
        ax.set_xlabel("x (cm)", fontsize=13)
        ax.set_ylabel("y (cm)" if plane == "axial" else "z (cm)", fontsize=13)
        ax.tick_params(labelsize=11)
        for s in ax.spines.values():
            s.set_visible(False)
    ax = fig.add_subplot(gs[1, :])
    c = np.array(m["currents_re"]) + 1j * np.array(m["currents_im"])
    c = np.real(c) / np.max(np.abs(c))
    if m.get("cap_free"):  # show the true (tiny) capacitor current relative to the strongest mode
        ref = max(mm.get("cap_current_sq", 0) for mm in res["modes"])
        c = c * np.sqrt(m["cap_current_sq"] / ref)
    n = np.arange(N)
    col = plots.ORANGE if m["homogeneous"] else plots.BLUE
    if partner is not None:  # quadrature: show both members of the pair (cos and sin patterns)
        c2 = np.array(partner["currents_re"]) + 1j * np.array(partner["currents_im"])
        c2 = np.real(c2) / np.max(np.abs(c2))
        ax.bar(n - 0.2, c, width=0.4, color=col, edgecolor="white", label="mode a")
        ax.bar(n + 0.2, c2, width=0.4, color="#f4b08f", edgecolor="white", label="mode b")
        ax.legend(frameon=False, fontsize=11, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.25))
    else:
        ax.bar(n, c, width=0.7, color=col, edgecolor="white")
    if m.get("cap_free"):
        ax.text(N / 2 - 0.5, 0.45, f"capacitor currents ≈ 0  ({np.sqrt(m['cap_current_sq'] / ref) * 100:.1f} % of the strongest mode)",
                ha="center", fontsize=13, color=plots.INK2)
    ax.axhline(0, color=plots.MUTED, lw=0.6)
    ax.set_xlim(-0.6, N - 0.4)
    ax.set_ylim(-1.15, 1.15)
    ax.set_xticks(range(0, N, 2))
    ax.set_xlabel("rung n", fontsize=13)
    ax.set_ylabel("capacitor\ncurrent", fontsize=13)
    ax.tick_params(labelsize=11)
    fig.text(0.5, 0.955, f"{lab} / plane max  ·  dashed = ROI", ha="center", fontsize=14, color=plots.INK2)
    fig.savefig(out, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out


def build(case: Path, out: Path, ffmpeg: str, workdir: Path, first_case: Path | None = None):
    workdir.mkdir(parents=True, exist_ok=True)
    res = json.loads((case / "modes.json").read_text())
    pr = json.loads((case / "probes.json").read_text())
    raw = load(case)
    pts = np.array(pr["points"])
    na, nc = pr["axial"], pr["coronal"]
    sl = {"axial": slice(0, na), "coronal": slice(na, na + nc), "roi": slice(na + nc, None)}
    des = res["design"]
    N = res["N"]
    hp = res["homogeneous_pair"]
    f0 = des["b0"] * 42.577478
    figs = case / "figs"

    by = {}
    for m in res["modes"]:
        if m.get("cap_free"):
            by.setdefault("capfree", m)
        elif m.get("k") is not None and m.get("k_purity", 0) > 0.8:
            by.setdefault(m["k"], m)
    hom = [m for m in res["modes"] if m["homogeneous"]][0]
    seq = []
    for key, title, note in ((N // 2, f"k = {N // 2}", f"neighbouring rungs opposite: field stays at the conductors"),
                             (2, "k = 2", "two current periods: field is zero at the centre"),
                             (0, "k = 0", "all rungs in phase, returning through the shield"),
                             ("capfree", "ring mode", "capacitors carry no current: an end-ring / stray-C resonance")):
        if key in by:
            seq.append((by[key], title, note))
    panels = []
    for i, (m, title, note) in enumerate(seq):
        panels.append((mode_panel(case, res, raw, pts, sl, m, workdir / f"mode_{i}.png"), m, title, note))
    hom_png = mode_panel(case, res, raw, pts, sl, hom, workdir / "mode_k1.png")
    q = quadrature(raw["B"][hp["m"][0] - 1], raw["B"][hp["m"][1] - 1], sl["roi"])
    partner = [m for m in res["modes"] if m["homogeneous"] and m["m"] != hom["m"]]
    quad_png = mode_panel(case, res, raw, pts, sl, hom, workdir / "mode_quad.png", quad=q, partner=partner[0] if partner else None)

    base = paper()
    cache = {}

    def get_plate(p, bw, bh):
        k = (str(p), bw, bh)
        if k not in cache:
            cache[k] = plate(Path(p), bw, bh)
        return cache[k]

    def below(p, bw, bh, y, gap=34):
        return y + get_plate(p, bw, bh).height + gap

    # video-sized spectrum (bigger type than the README figure)
    with plt.rc_context({"font.size": 15, "axes.titlesize": 16, "legend.fontsize": 13}):
        if first_case:
            spec_png = plots.fig_spectrum_compare(first_case, case, workdir / "spectrum_video.png", figsize=(8.4, 10.0), vertical=True)
        else:
            spec_png = plots.fig_spectrum(case, workdir / "spectrum_video.png", figsize=(8.2, 7.4))

    scenes = []
    coil_png = figs / "coil.png"
    y1 = below(coil_png, 940, 830, 250)
    scenes.append(Scene(5.0, [
        (0.0, "text", ("A birdcage coil has many resonant modes.", SERIF(48)), (70, 92)),
        (0.35, "text", ("Only one gives a uniform transverse B1.", SERIF_I(48)), (70, 156)),
        (0.8, "plate", (coil_png, 940, 830), (70, 250)),
        (1.3, "text", (f"{N}-rung low-pass head coil · Ø {des['diameter'] * 100:.0f} cm · RF shield Ø {des['shield_diameter'] * 100:.0f} cm", SANS(28)), (70, y1)),
        (1.5, "text", (f"target: ¹H at {des['b0']:g} T = {f0:.1f} MHz", SANS(28)), (70, y1 + 44)),
    ]))
    nmodes = len(res["modes"])
    y2 = below(spec_png, 940, 940, 230, gap=26)
    scenes.append(Scene(7.0, [
        (0.0, "text", (f"Palace eigenmode solve: {nmodes} modes", SERIF(48)), (70, 92)),
        (0.3, "text", ("ladder-network model (open) vs Palace FEM (filled)", SANS(28)), (70, 160)),
        (0.6, "plate", (spec_png, 940, 940), (70, 230)),
        (1.2, "text", ("k = harmonic of the rung currents, Iₙ ∝ cos(2πkn/N)", SERIF(30)), (70, y2)),
        (1.4, "text", ("every k between 1 and N/2 − 1 is a degenerate pair", SERIF(30)), (70, y2 + 46)),
    ]))
    for png, m, title, note in panels:
        y3 = below(png, 940, 920, 225)
        scenes.append(Scene(3.0, [
            (0.0, "text", (f"{title}", SERIF(58)), (70, 80)),
            (0.0, "text", (f"{m['f_mhz']:.1f} MHz", MONO(42)), (W - 70, 96, "ra")),
            (0.15, "text", (note, SANS(28)), (70, 162)),
            (0.25, "plate", (png, 940, 920), (70, 225)),
            (0.5, "text", (f"ROI std/mean of |B_t|: {m['roi_bt_cv'] * 100:.0f} %", MONO(32)), (70, y3)),
        ]))
    y4 = below(hom_png, 940, 920, 225)
    scenes.append(Scene(4.5, [
        (0.0, "text", ("k = 1", SERIF(58)), (70, 80)),
        (0.0, "text", (f"{hom['f_mhz']:.1f} MHz", MONO(42)), (W - 70, 96, "ra")),
        (0.15, "text", ("sinusoidal rung currents: uniform transverse field", SANS(28)), (70, 162)),
        (0.25, "plate", (hom_png, 940, 920), (70, 225)),
        (0.5, "text", (f"ROI std/mean of |B_t|: {hom['roi_bt_cv'] * 100:.1f} %", MONO(32)), (70, y4)),
    ]))
    cl = des.get("c_leg")
    fg = des.get("notes", {}).get("c_leg_first_guess")
    y5 = below(quad_png, 940, 700, 380, gap=28)
    scenes.append(Scene(7.0, [
        (0.0, "text", ("The homogeneous pair, tuned", SERIF(48)), (70, 92)),
        (0.3, "big", (f"{hp['f_mean_mhz']:.1f} MHz", f"k = 1 pair · target {f0:.1f} MHz"), (70, 175)),
        (0.5, "big", (f"{cl * 1e12:.2f} pF" if cl else "", f"per rung · ladder first guess {fg * 1e12:.2f} pF" if fg else "per rung"), (560, 175)),
        (0.7, "plate", (quad_png, 940, 700), (70, 380)),
        (1.0, "text", (f"quadrature B1+: ROI std/mean {hp['b1plus_roi_cv'] * 100:.1f} %, |B1−|/|B1+| {hp['b1minus_over_b1plus_roi'] * 100:.1f} %", MONO(26)), (70, y5)),
        (1.1, "text", (f"pair split {hp['split_mhz'] * 1e3:.0f} kHz · ROI r < 0.6 R, |z| < L/6", MONO(26)), (70, y5 + 40)),
        (1.3, "text", ("Simulation only · Palace FEM eigenmode · PEC strips, lumped C, closed shield, unloaded", SANS(22)), (70, y5 + 92)),
    ]))
    scenes.append(Scene(3.0, [
        (0.0, "text", ("birdcage-palace", SERIF(68)), (W // 2, 560, "ma")),
        (0.4, "text", ("an open-source birdcage coil builder for Palace", SANS(32)), (W // 2, 664, "ma")),
    ]))

    cmd = [ffmpeg, "-y", "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", f"{W}x{H}", "-r", str(FPS), "-i", "-",
           "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18", "-preset", "slow", "-movflags", "+faststart", str(out)]
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stderr=subprocess.DEVNULL)
    fade = 0.35
    stills = {}
    for si, sc in enumerate(scenes):
        nf = int(sc.dur * FPS)
        for fi in range(nf):
            t = fi / FPS
            im = base.copy()
            layer = Image.new("RGBA", (W, H), (0, 0, 0, 0))
            for t_in, kind, payload, xy in sc.items:
                a = ease((t - t_in) / 0.6)
                if a <= 0:
                    continue
                dy = int((1 - a) * 18)
                tmp = Image.new("RGBA", (W, H), (0, 0, 0, 0))
                d = ImageDraw.Draw(tmp)
                if kind == "text":
                    txt, f = payload
                    anchor = xy[2] if len(xy) > 2 else "la"
                    d.text((xy[0], xy[1] + dy), txt, font=f, fill=COBALT + (255,), anchor=anchor)
                elif kind == "big":
                    big, small = payload
                    d.text((xy[0], xy[1] + dy), big, font=SERIF(72), fill=COBALT + (255,))
                    d.text((xy[0], xy[1] + 92 + dy), small, font=MONO(22), fill=COBALT2 + (255,))
                elif kind == "plate":
                    pl = get_plate(*payload)
                    tmp.paste(pl, (xy[0] + (payload[1] - pl.width) // 2, xy[1] + dy))
                if a < 1:
                    alpha = tmp.getchannel("A").point(lambda v, a=a: int(v * a))
                    tmp.putalpha(alpha)
                layer = Image.alpha_composite(layer, tmp)
            im = Image.alpha_composite(im.convert("RGBA"), layer)
            d = ImageDraw.Draw(im)
            d.text((40, H - 50), "BIRDCAGE-PALACE", font=MONO(22), fill=COBALT + (255,))
            d.text((W - 40, H - 50), f"{si + 1:02d} / {len(scenes):02d}", font=MONO(22), fill=COBALT + (255,), anchor="ra")
            # fade out at the scene end (not the last)
            if si < len(scenes) - 1 and t > sc.dur - fade:
                k = (t - (sc.dur - fade)) / fade
                im = Image.blend(im.convert("RGB"), base, k * 0.85).convert("RGBA")
            rgb = im.convert("RGB")
            if (si == len(scenes) - 1 and fi == nf - 1) or (si < len(scenes) - 1 and fi == int((sc.dur - fade) * FPS) - 1):
                stills[si] = rgb
            proc.stdin.write(rgb.tobytes())
    proc.stdin.close()
    proc.wait()
    for si, img in stills.items():
        img.save(workdir / f"still_{si + 1:02d}.png")
    total = sum(s.dur for s in scenes)
    print(f"wrote {out} ({total:.1f} s, {len(scenes)} scenes)")
    return stills


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ffmpeg", default="ffmpeg")
    ap.add_argument("--work", default=None)
    ap.add_argument("--first-case", default=None, help="same mesh solved with the ladder first-guess C (spectrum panel a)")
    a = ap.parse_args()
    case = Path(a.case)
    build(case, Path(a.out), a.ffmpeg, Path(a.work) if a.work else case / "video", Path(a.first_case) if a.first_case else None)

