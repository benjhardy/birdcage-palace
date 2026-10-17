"""Hero image (1600x900): coil render next to the k=1 quadrature |B1+| map and the k=2 map for contrast."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from PIL import Image  # noqa: E402

from birdcage_palace import plots  # noqa: E402
from birdcage_palace.modes import load, quadrature  # noqa: E402


def main(case: Path, out: Path):
    res = json.loads((case / "modes.json").read_text())
    pr = json.loads((case / "probes.json").read_text())
    raw = load(case)
    pts = np.array(pr["points"])
    na = pr["axial"]
    ax_sl = slice(0, na)
    roi = slice(na + pr["coronal"], None)
    hp = res["homogeneous_pair"]
    q = quadrature(raw["B"][hp["m"][0] - 1], raw["B"][hp["m"][1] - 1], roi)
    k2 = [m for m in res["modes"] if m.get("k") == 2][0]
    B2 = raw["B"][k2["m"] - 1]
    bt2 = np.sqrt(np.abs(B2[:, 0]) ** 2 + np.abs(B2[:, 1]) ** 2)
    fig = plt.figure(figsize=(16, 9), dpi=100)
    a0 = fig.add_axes([0.02, 0.02, 0.44, 0.80])
    a0.imshow(Image.open(case.parents[1] / "docs" / "img" / "coil.png"))
    a0.axis("off")
    for i, (v, title) in enumerate(((np.abs(q["b1p"]), f"k = 1 pair, quadrature |B1+|\n{hp['f_mean_mhz']:.1f} MHz"),
                                    (bt2, f"k = 2, |B_t|\n{k2['f_mhz']:.1f} MHz"))):
        ax = fig.add_axes([0.48 + i * 0.26, 0.12, 0.24, 0.66])
        img, ext = plots._grid(pts[ax_sl][:, :2], v[ax_sl] / np.nanmax(v[ax_sl]))
        ax.imshow(img, origin="lower", extent=[e * 100 for e in ext], cmap=plots.SEQ, vmin=0, vmax=1, interpolation="bilinear")
        plots._coil_marks(ax, res["design"], "axial")
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(title, fontsize=17, color=plots.ORANGE if i == 0 else plots.INK, pad=10)
    fig.text(0.03, 0.93, "birdcage-palace", fontsize=30, color="#1f3fae", family="serif")
    fig.text(0.03, 0.885, "16-rung low-pass head coil at 3 T, solved with the Palace FEM eigenmode solver",
             fontsize=15, color=plots.INK2)
    fig.text(0.48, 0.05, "axial plane z = 0, field / plane max · dashed circle = ROI · simulation only",
             fontsize=12, color=plots.INK2)
    fig.savefig(out, dpi=100, facecolor="white")
    print("wrote", out)


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
