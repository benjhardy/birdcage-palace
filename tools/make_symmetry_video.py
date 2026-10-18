"""Silent 4:5 (1080x1350) clip for the symmetry-plane follow-up: full cage -> mirror planes -> quarter mesh ->
same homogeneous mode -> cost drop -> converged result. Every number comes from results/lp16_symmetry.json.

    python tools/make_symmetry_video.py --out docs/img/birdcage_symmetry.mp4 --ffmpeg <ffmpeg>

Same look as tools/make_video.py (cream graph paper, cobalt ink, plots in their original colours on white plates).
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import make_video as mv  # noqa: E402
import symmetry_report as sr  # noqa: E402
from birdcage_palace import plots  # noqa: E402

W, H, FPS = mv.W, mv.H, mv.FPS


def build(out: Path, ffmpeg: str, work: Path):
    work.mkdir(parents=True, exist_ok=True)
    R = json.loads((ROOT / "results" / "lp16_symmetry.json").read_text())
    full, q = R["full"], R["combos"]["em"]
    fc, qc = full["cost"], q["cost"]
    ec = R.get("eighth", {}).get("cost")
    qcase = Path(R["quarter_prefix"] + "_em")
    # renders
    full_png = sr.render_model(qcase, work / "v_full.png", "full", window=(1400, 1400))
    cut_png = sr.render_model(qcase, work / "v_cut.png", "cut", window=(1400, 1400))
    qmesh_png = sr.render_model(qcase, work / "v_quarter_mesh.png", "quarter", show_mesh=True, ghost=False, window=(1400, 1400))
    with plt.rc_context({"font.size": 14}):
        b1_png, _ = sr.fig_b1_compare(Path(R["full"]["case"]), qcase, Path(R["quarter_prefix"] + "_me"), work / "v_b1.png",
                                      figsize=(8.6, 9.4), big=True, linear_only=True)
    rows = [("full", fc), ("quarter", qc)]  # the quarter was timed on a free machine; full = the production run
    cost_png = sr.fig_costs(rows, work / "v_costs.png", figsize=(9.6, 6.6), big=True, vertical=True)
    conv_png = (sr.fig_convergence(R["convergence"], work / "v_conv.png", R, figsize=(8.6, 6.6), fs=1.35)
                if R.get("convergence") else None)

    S, T, TI, SB, MO = mv.SANS, mv.SERIF, mv.SERIF_I, mv.SANS_B, mv.MONO
    Scene = mv.Scene
    scenes = []
    scenes.append(Scene(4.5, [
        (0.0, "text", ("A full birdcage model", T(50)), (70, 92)),
        (0.3, "text", ("16-rung low-pass head coil, Palace eigenmode solve", S(30)), (70, 162)),
        (0.6, "plate", (full_png, 940, 860), (70, 230)),
        (1.1, "text", (f"{fc['unknowns'] / 1e6:.2f} M unknowns · {fc['palace_gb']:.1f} GB · {fc['palace_s'] / 60:.1f} min on {fc['ranks']} ranks", MO(30)), (70, 1130)),
    ]))
    scenes.append(Scene(6.0, [
        (0.0, "text", ("Two mirror planes", T(50)), (70, 92)),
        (0.3, "text", ("for one member of the k = 1 pair (B1 along y)", S(30)), (70, 162)),
        (0.6, "plate", (cut_png, 940, 800), (70, 225)),
        (1.2, "text", ("x = 0: rung currents antisymmetric → PEC, tangential E = 0", S(30)), (70, 1060)),
        (1.6, "text", ("y = 0: rung currents symmetric → PMC, tangential H = 0", S(30)), (70, 1108)),
        (2.0, "text", ("planes between rungs (cage rotated by π/16): no capacitor is cut", S(26)), (70, 1166)),
    ]))
    scenes.append(Scene(4.5, [
        (0.0, "text", ("The quarter model", T(50)), (70, 92)),
        (0.3, "text", ("PEC set explicitly; PMC is Palace's natural boundary", S(30)), (70, 162)),
        (0.6, "plate", (qmesh_png, 940, 860), (70, 230)),
        (1.0, "text", (f"{qc['unknowns'] / 1e6:.2f} M unknowns, same mesh density", MO(30)), (70, 1130)),
    ]))
    df = (q["k1_mhz"] / full["k1_mean_mhz"] - 1) * 100
    scenes.append(Scene(6.0, [
        (0.0, "text", ("The same homogeneous mode", T(50)), (70, 92)),
        (0.3, "big", (f"{q['k1_mhz']:.2f} MHz", f"quarter · full model {full['k1_mean_mhz']:.2f} MHz ({df:+.2f} %)"), (70, 170)),
        (0.5, "big", (f"{q['roi_bt_cv'] * 100:.1f} %", f"ROI std/mean |B_t| · full {full['roi_bt_cv'][0] * 100:.1f} %"), (640, 170)),
        (0.8, "plate", (b1_png, 940, 860), (70, 330)),
    ]))
    tq = fc["factor_s"] / qc["factor_s"]
    scenes.append(Scene(5.0, [
        (0.0, "text", ("Cost of one eigen solve", T(50)), (70, 92)),
        (0.3, "text", (f"quarter: {fc['unknowns'] / qc['unknowns']:.1f}× fewer unknowns, {fc['palace_gb'] / qc['palace_gb']:.1f}× less memory, {tq:.0f}× faster factorisation", S(28)), (70, 162)),
        (0.5, "text", ("full: 20 modes · quarter: 10 modes (k = 1 is its lowest) · free machine, 4 ranks each", S(24)), (70, 900)),
        (0.6, "plate", (cost_png, 940, 900), (70, 230)),
    ]))
    if conv_png and R.get("converged"):
        cv = R["converged"]
        scenes.append(Scene(5.5, [
            (0.0, "text", ("Spending the savings", T(50)), (70, 92)),
            (0.3, "text", ("eighth model (z = 0 is PEC too): up to 1.0 M unknowns, ~8 M in the full coil", S(28)), (70, 162)),
            (0.6, "plate", (conv_png, 940, 760), (70, 225)),
            (1.0, "big", (f"{cv['k1_mhz']:.2f} MHz", "k = 1 extrapolated from the strip-edge refinement"), (70, 1010)),
            (1.3, "text", (cv["note"], MO(25)), (70, 1150)),
        ]))
    scenes.append(Scene(3.0, [
        (0.0, "text", ("birdcage-palace", T(68)), (W // 2, 540, "ma")),
        (0.4, "text", ("symmetry planes for Palace eigenmode solves", S(32)), (W // 2, 644, "ma")),
        (0.7, "text", ("Simulation only", S(26)), (W // 2, 720, "ma")),
    ]))
    render(scenes, out, ffmpeg, work)


def render(scenes, out, ffmpeg, work):
    base = mv.paper()
    cache = {}

    def get_plate(p, bw, bh):
        k = (str(p), bw, bh)
        if k not in cache:
            cache[k] = mv.plate(Path(p), bw, bh)
        return cache[k]
    cmd = [ffmpeg, "-y", "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", f"{W}x{H}", "-r", str(FPS), "-i", "-",
           "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18", "-preset", "slow", "-movflags", "+faststart", str(out)]
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stderr=subprocess.DEVNULL)
    fade = 0.35
    for si, sc in enumerate(scenes):
        nf = int(sc.dur * FPS)
        for fi in range(nf):
            t = fi / FPS
            layer = Image.new("RGBA", (W, H), (0, 0, 0, 0))
            for t_in, kind, payload, xy in sc.items:
                a = mv.ease((t - t_in) / 0.6)
                if a <= 0:
                    continue
                dy = int((1 - a) * 18)
                tmp = Image.new("RGBA", (W, H), (0, 0, 0, 0))
                d = ImageDraw.Draw(tmp)
                if kind == "text":
                    txt, f = payload
                    d.text((xy[0], xy[1] + dy), txt, font=f, fill=mv.COBALT + (255,), anchor=xy[2] if len(xy) > 2 else "la")
                elif kind == "big":
                    big, small = payload
                    d.text((xy[0], xy[1] + dy), big, font=mv.SERIF(72), fill=mv.COBALT + (255,))
                    d.text((xy[0], xy[1] + 92 + dy), small, font=mv.MONO(22), fill=mv.COBALT2 + (255,))
                elif kind == "plate":
                    pl = get_plate(*payload)
                    tmp.paste(pl, (xy[0] + (payload[1] - pl.width) // 2, xy[1] + dy))
                if a < 1:
                    tmp.putalpha(tmp.getchannel("A").point(lambda v, a=a: int(v * a)))
                layer = Image.alpha_composite(layer, tmp)
            im = Image.alpha_composite(base.convert("RGBA"), layer)
            d = ImageDraw.Draw(im)
            d.text((40, H - 50), "BIRDCAGE-PALACE", font=mv.MONO(22), fill=mv.COBALT + (255,))
            d.text((W - 40, H - 50), f"{si + 1:02d} / {len(scenes):02d}", font=mv.MONO(22), fill=mv.COBALT + (255,), anchor="ra")
            if si < len(scenes) - 1 and t > sc.dur - fade:
                k = (t - (sc.dur - fade)) / fade
                im = Image.blend(im.convert("RGB"), base, k * 0.85).convert("RGBA")
            rgb = im.convert("RGB")
            if fi == int((sc.dur - fade) * FPS) - 1 or (si == len(scenes) - 1 and fi == nf - 1):
                rgb.save(work / f"still_{si + 1:02d}.png")
            proc.stdin.write(rgb.tobytes())
    proc.stdin.close()
    proc.wait()
    print(f"wrote {out} ({sum(s.dur for s in scenes):.1f} s, {len(scenes)} scenes)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--ffmpeg", default="ffmpeg")
    ap.add_argument("--work", default=str(ROOT / "runs" / "sym_video"))
    a = ap.parse_args()
    build(Path(a.out), a.ffmpeg, Path(a.work))
