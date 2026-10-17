"""Copy the figures of a solved case into docs/img (cropped), and write results/<name>.json with the key numbers.

    python tools/export_docs.py --case runs/lp16_final --first-case runs/lp16_m2 --name lp16_3T [--also <dir>]
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def crop(src: Path, dst: Path, pad: int = 24) -> None:
    im = Image.open(src).convert("RGB")
    a = np.asarray(im)
    ys, xs = np.where((a < 245).any(2))
    if len(xs):
        im = im.crop((max(xs.min() - pad, 0), max(ys.min() - pad, 0), min(xs.max() + pad, im.width), min(ys.max() + pad, im.height)))
    im.save(dst)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True)
    ap.add_argument("--first-case", default=None)
    ap.add_argument("--name", required=True)
    ap.add_argument("--also", action="append", default=[])
    a = ap.parse_args()
    case = Path(a.case)
    img = ROOT / "docs" / "img"
    img.mkdir(parents=True, exist_ok=True)
    from birdcage_palace import plots
    figs = case / "figs"
    if a.first_case:
        plots.fig_spectrum_compare(Path(a.first_case), case, figs / "spectrum_compare.png")
    outs = []
    for name in ("coil.png", "mesh.png", "spectrum_compare.png", "spectrum.png", "b1_maps.png", "currents.png"):
        if (figs / name).exists():
            crop(figs / name, img / name)
            outs.append(img / name)
    vid = case / "birdcage_modes.mp4"
    if vid.exists():
        shutil.copy(vid, img / vid.name)
        outs.append(img / vid.name)
        th = case / "video" / "still_08.png"
        if th.exists():
            shutil.copy(th, img / "birdcage_modes_thumbnail.png")
            outs.append(img / "birdcage_modes_thumbnail.png")
    res = json.loads((case / "modes.json").read_text())
    summ = {"design": res["design"], "homogeneous_pair": res["homogeneous_pair"], "comparison": res["comparison"],
            "solve_stats": res["solve_stats"],
            "modes": [{k: m.get(k) for k in ("m", "f_mhz", "label", "k", "k_purity", "roi_bt_cv", "roi_transverse_frac",
                                              "cap_free", "homogeneous")} for m in res["modes"]]}
    for extra in ("tune_history_part1.json", "tune_history.json"):
        if (case / extra).exists():
            summ.setdefault("tune_history", []).extend(json.loads((case / extra).read_text()))
    if a.first_case:
        r0 = json.loads((Path(a.first_case) / "modes.json").read_text())
        summ["first_guess"] = {"c_leg": r0["design"]["c_leg"], "comparison": r0["comparison"],
                               "homogeneous_pair": r0["homogeneous_pair"], "solve_stats": r0["solve_stats"]}
    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results" / f"{a.name}.json").write_text(json.dumps(summ, indent=1))
    for d in a.also:
        d = Path(d)
        d.mkdir(parents=True, exist_ok=True)
        for p in outs:
            shutil.copy(p, d / p.name)
    print("exported:", ", ".join(p.name for p in outs), "+ results/" + a.name + ".json")


if __name__ == "__main__":
    main()
