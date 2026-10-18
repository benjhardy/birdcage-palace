"""Silent 4:5 (1080x1350) clip of the port-tuning step: the 16-port of the bare coil, then circuitRF's minimax
iterations pulling the drive-port |S11| and |S21| onto 127.7 MHz, then the tuned values.

    python tools/make_porttune_video.py --tune runs/ptl_e/circuitrf/tune_loaded.json --out docs/img/birdcage_porttune.mp4

Each frame's curves are the 2-port response at circuitRF's iterate, evaluated from the same 16-port S-matrix by
tools/port_tune.python_two_port (agrees with circuitRF's own sweep to ~1e-5).
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
import port_tune as pt  # noqa: E402
import symmetry_report as sr  # noqa: E402
from birdcage_palace import plots  # noqa: E402

W, H, FPS = mv.W, mv.H, mv.FPS


def frame_plot(f, S2, f0, out: Path, it_label: str, band=None) -> Path:
    with plt.rc_context({"font.size": 15}):
        fig, ax = plt.subplots(figsize=(8.6, 6.4))
        ax.plot(f / 1e6, 20 * np.log10(np.abs(S2[:, 0, 0])), color=plots.BLUE, lw=2.2, label="|S11|, drive port 1")
        ax.plot(f / 1e6, 20 * np.log10(np.abs(S2[:, 1, 0]) + 1e-12), color=plots.ORANGE, lw=1.8, label="|S21|, port 1 to port 2")
        ax.axvline(f0 / 1e6, color=plots.MUTED, ls=":", lw=1)
        ax.axhline(-20, color=plots.MUTED, ls=":", lw=1)
        if band:
            ax.axvspan(band[0] / 1e6, band[1] / 1e6, color=plots.GRID, alpha=0.8, lw=0)
        ax.set_ylim(-50, 2)
        ax.set_xlim(f[0] / 1e6, f[-1] / 1e6)
        ax.set_xlabel("frequency (MHz)")
        ax.set_ylabel("dB")
        ax.grid(color=plots.GRID, lw=0.6)
        ax.legend(frameon=False, loc="lower left", fontsize=13)
        ax.set_title(it_label, loc="left", fontsize=15, color=plots.INK2)
        fig.savefig(out, dpi=130, bbox_inches="tight", facecolor="white")
        plt.close(fig)
    return out


def build(tune_json: Path, out: Path, ffmpeg: str, work: Path, coil_case: Path):
    work.mkdir(parents=True, exist_ok=True)
    R = json.loads(tune_json.read_text())
    f0 = R["f0_hz"]
    fe, Se, r = pt.read_snp(Path(R["snp"]).with_name(Path(R["snp"]).stem + "_sym" + Path(R["snp"]).suffix))
    f = np.arange(f0 - 2.5e6, f0 + 2.5e6, 5e3)
    band = (f0 - R["band_khz"] * 1e3, f0 + R["band_khz"] * 1e3)
    traj = R["trajectory"]
    frames = []
    v0 = R["start"]  # tune()'s starting point
    while len(traj) > 2 and traj[-1]["vals"] == traj[-2]["vals"]:  # the optimizer has stopped moving
        traj = traj[:-1]
    seq = [("start", v0)] + [(f"iteration {t['iter']}", t["vals"]) for t in traj]
    for i, (lab, vals) in enumerate(seq):
        S2 = pt.python_two_port(f, (fe, Se), r, vals, R["esr_ohm"])
        txt = (f"{lab}:  C {vals['Ct'] * 1e12:.3f} pF · tap {vals['Cp1'] * 1e12:.1f} pF · series {vals['Cm1'] * 1e12:.1f} pF"
               if "Cp1" in vals else
               f"{lab}:  C {vals['Ct'] * 1e12:.3f} pF · drive {vals['Cd1'] * 1e12:.3f} pF · match {vals['Cm1'] * 1e12:.3f} pF")
        frames.append(frame_plot(f, S2, f0, work / f"it_{i:02d}.png", txt, band))
    coil = sr.render_model(coil_case, work / "v_coil.png", "full", window=(1300, 1300))
    S, T, MO = mv.SANS, mv.SERIF, mv.MONO
    Scene = mv.Scene
    scenes = [Scene(5.0, [
        (0.0, "text", ("Port tuning", T(56)), (70, 92)),
        (0.3, "text", ("every capacitor gap becomes a port: a 16-port of the bare coil", S(30)), (70, 168)),
        (0.6, "plate", (coil, 940, 760), (70, 240)),
        (1.0, "text", ("Palace driven sweep, 100-160 MHz, 4 eighth models × 4 ports", MO(27)), (70, 1050)),
        (1.3, "text", ("mirror symmetry rebuilds the full 16 × 16 S-matrix", MO(27)), (70, 1094)),
    ])]
    dt = max(0.55, min(1.2, 9.0 / len(frames)))
    items = [(0.0, "text", ("circuitRF closes the gaps", T(50)), (70, 92)),
             (0.2, "text", ("capacitors (ESR), capacitive-divider feed at gaps 0 and 4, 50 Ω", S(28)), (70, 160))]
    scenes.append(Scene(dt * len(frames) + 2.0, items, ))
    scenes[-1].frames = [(i * dt + 0.5, p) for i, p in enumerate(frames)]
    v = R["raw"]["vals"]
    gl = {g["name"]: g["value"] for g in R["raw"]["goals"]}
    scenes.append(Scene(5.5, [
        (0.0, "text", ("Tuned", T(56)), (70, 92)),
        (0.3, "big", (f"{v['Ct'] * 1e12:.3f} pF", f"12 rungs, loaded · the unloaded eigenmode loop gave {R['c_eigen'] * 1e12:.3f} pF"), (70, 190)),
        (0.5, "big", (f"{v['Cd1'] * 1e12:.2f} + {v['Cp1'] * 1e12:.1f}", "pF divider at drive gap 1 (top + tap)"), (70, 390)),
        (0.7, "big", (f"{v['Cm1'] * 1e12:.1f} pF", "series capacitor to the 50 Ω port"), (70, 590)),
        (0.9, "text", (f"worst over 127.7 MHz ± {R['band_khz']:g} kHz:", MO(26)), (70, 820)),
        (0.9, "text", (f"|S11| {gl['match1']:.1f} dB · |S22| {gl['match2']:.1f} dB · |S21| {gl['isol']:.1f} dB", MO(26)), (70, 860)),
        (1.0, "text", ("last step: per-port trims on the matrix as solved", S(26)), (70, 920)),
        (1.1, "text", ("Simulation only · saline phantom load · lossless conductors", S(26)), (70, 970)),
    ]))
    scenes.append(Scene(3.0, [
        (0.0, "text", ("birdcage-palace", T(68)), (W // 2, 540, "ma")),
        (0.4, "text", ("port tuning with Palace and circuitRF", S(32)), (W // 2, 644, "ma")),
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
        fr = getattr(sc, "frames", None)
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
            if fr:  # stepped plot frames (no motion within a step)
                cur = [p for ts, p in fr if ts <= t]
                if cur:
                    pl = get_plate(cur[-1], 940, 980)
                    layer.paste(pl, (70 + (940 - pl.width) // 2, 240))
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
    ap.add_argument("--tune", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ffmpeg", default="ffmpeg")
    ap.add_argument("--coil-case", default=str(ROOT / "runs" / "sym_q_em"))
    ap.add_argument("--work", default=str(ROOT / "runs" / "porttune_video"))
    a = ap.parse_args()
    build(Path(a.tune), Path(a.out), a.ffmpeg, Path(a.work), Path(a.coil_case))
