"""Field uniformity per mode: the quasi-static model (birdcage_palace/uniformity.py and its port in the web builder)
against itself and against the Palace eigenmode results of the three 16-rung 3 T reference coils.

    python -m pytest tests/test_uniformity.py -s      # -s prints the model-vs-Palace table

Metric: std/mean of the transverse |B| over the ROI (r < 0.6 R, |z| < length/6) and the transverse fraction
sum |B_t|^2 / sum |B|^2, as birdcage_palace.modes computes them from Palace's probes. Palace numbers are the mean
over each degenerate pair (results/lp16_3T.json, hp16_3T.json, bp16_3T.json; production mesh, order 2, closed can).

Result (model / Palace LP, HP, BP), CV in %:
    k = 1           3.2 / 4.0, 3.4, 3.4
    k = 2          36.0 / 35.4, 36.1, 36.0
    k = 3          58.1 / 57.2, 57.8, 57.9
    k = 4          75.5 / 74.6, 75.0, 75.1
    k = 5..8       90-126 / 89-126 (within 2.5 % relative)
    ring k = 0 counter-rotating   33.9 / -, 34.2, 34.4   (transverse fraction 0.64 / 0.62)
    ring co-rotating              72.9 / -, 72.8, 72.2   (transverse fraction 0.42 / 0.065: the model gets
                                  "mostly B_z" right but overestimates B_t; the image rule is a line-current rule)
The quasi-static model gives the same currents for every coil type, so it cannot see the LP-vs-HP difference at
k = 1 (4.0 vs 3.4 %: stray capacitance along the low-pass rungs). Not in the ladder model and not coloured on the
page: the end-ring k = 1 family of HP/BP (co-rotating rings, 263-267 MHz, Palace 60 %) and the capacitor-free
modes (stray C).
"""
import json
import math
import random
import shutil
import subprocess
from pathlib import Path

import pytest

from birdcage_palace import uniformity

ROOT = Path(__file__).resolve().parents[1]
PAGE = ROOT / "web" / "birdcage-builder" / "index.html"
REFS = {"lp": "lp16_3T", "hp": "hp16_3T", "bp": "bp16_3T"}

JS = r"""
const fs = require('fs'), vm = require('vm');
const html = fs.readFileSync(process.argv[process.argv.length - 1], 'utf8');
const core = html.match(/<script id="core">([\s\S]*?)<\/script>/)[1];
const ctx = {}; vm.createContext(ctx); vm.runInContext(core + '\nthis.Uniformity = Uniformity; this.PALACE_U = PALACE_U;', ctx);
const designs = JSON.parse(fs.readFileSync(0, 'utf8'));
process.stdout.write(JSON.stringify({ out: designs.map(d => ctx.Uniformity.spectrum(d)), palace: ctx.PALACE_U }));
"""


def run_js(designs):
    js = subprocess.run([shutil.which("node"), "-e", JS, str(PAGE)], input=json.dumps(designs), capture_output=True,
                        text=True, check=True)
    return json.loads(js.stdout)


def palace_table(name: str) -> dict:
    """Mean CV and transverse fraction per mode family of one Palace result file."""
    acc = {}
    for m in json.loads((ROOT / "results" / f"{name}.json").read_text())["modes"]:
        lab, rt = m["label"], m.get("ring_bot_over_top", 1.0)
        if "rungs in phase" in lab:
            key = "rungs_in_phase"
        elif lab.startswith("k="):
            key = ("endring_k" if rt < 0 else "k") + lab[2:]
        elif "counter-rotating" in lab:
            key = "ring0_counter"
        elif "co-rotating" in lab:
            key = "ring0_corot"
        else:  # no capacitor current (stray-C modes), one key per family
            key = f"stray{int(m['f_mhz'])}"
        acc.setdefault(key, []).append((m["roi_bt_cv"], m["roi_transverse_frac"], m["f_mhz"]))
    return {k: {"cv": sum(v[0] for v in vs) / len(vs), "tf": sum(v[1] for v in vs) / len(vs),
                "f_mhz": [min(v[2] for v in vs), max(v[2] for v in vs)]} for k, vs in acc.items()}


def ref_design(name):
    return json.loads((ROOT / "results" / f"{name}.json").read_text())["design"]


def random_designs(n=10, seed=20261023):
    rng = random.Random(seed)
    out = []
    for i in range(n):
        dia = rng.uniform(0.1, 0.6)
        out.append(dict(legs=rng.choice(range(8, 33, 2)), diameter=dia, length=rng.uniform(0.6, 1.3) * dia,
                        shield_diameter=None if i % 4 == 3 else dia + rng.uniform(0.02, 0.14),
                        shield_length=None if i % 2 else dia * rng.uniform(1.2, 2.0),
                        shield_ends="absorbing" if i % 3 == 2 else "pec"))
    return out


@pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")
def test_js_matches_python():
    designs = random_designs() + [ref_design(n) for n in REFS.values()]
    js = run_js(designs)["out"]
    worst = 0.0
    for d, rj in zip(designs, js):
        rp = uniformity.spectrum_uniformity(d)
        assert rj["n_points"] == rp["n_points"]
        for a, b in zip(rp["modes"] + [rp["corot"]], rj["modes"] + [rj["corot"]]):
            for key in ("bt_cv", "transverse_frac"):
                worst = max(worst, abs(a[key] - b[key]) / max(abs(a[key]), 1e-12))
    print(f"\n{len(designs)} designs, worst relative difference JS vs Python {worst:.1e}")
    assert worst < 1e-9


def test_roi_matches_palace_probes():
    from birdcage_palace.design import CoilDesign
    from birdcage_palace.palace import probe_points
    import numpy as np
    d = ref_design("lp16_3T")
    des = CoilDesign(**{k: v for k, v in d.items() if k in CoilDesign.__dataclass_fields__})
    a, b = np.array(probe_points(des)["roi"]), uniformity.roi_points(d)
    assert a.shape == b.shape and abs(a - b).max() < 1e-6


def test_mode_currents_are_kirchhoff():
    """Rung n carries cos(k phi_n) and every node conserves current (rung = difference of the two ring segments)."""
    N = 16
    for k in range(1, N // 2 + 1):
        w = uniformity.mode_weights(N, k)
        rung, top, bot = w[:N], w[N:2 * N], w[2 * N:]
        for n in range(N):
            assert rung[n] == pytest.approx(math.cos(2 * math.pi * k * n / N), abs=1e-12)
            assert rung[n] == pytest.approx(top[n] - top[n - 1], abs=1e-12)
            assert -rung[n] == pytest.approx(bot[n] - bot[n - 1], abs=1e-12)


def test_model_vs_palace():
    rows = []
    for t, name in REFS.items():
        pal, mod = palace_table(name), uniformity.spectrum_uniformity(ref_design(name))
        cv = {m["k"]: m["bt_cv"] for m in mod["modes"]}
        tf = {m["k"]: m["transverse_frac"] for m in mod["modes"]}
        # k = 1 is the most uniform mode by far, and the CV grows with k up to N/2
        assert cv[1] < 0.2 * min(v for k, v in cv.items() if k != 1)
        assert all(cv[k] < cv[k + 1] for k in range(1, 8))
        assert abs(cv[1] - pal["k1"]["cv"]) < 0.01                    # 1 percentage point
        for k in range(2, 9):
            assert cv[k] == pytest.approx(pal[f"k{k}"]["cv"], rel=0.05)
            assert tf[k] > 0.95 and pal[f"k{k}"]["tf"] > 0.95
        if t != "lp":
            assert cv[0] == pytest.approx(pal["ring0_counter"]["cv"], rel=0.05)
            assert tf[0] == pytest.approx(pal["ring0_counter"]["tf"], abs=0.05)
            assert mod["corot"]["bt_cv"] == pytest.approx(pal["ring0_corot"]["cv"], rel=0.05)
            assert mod["corot"]["transverse_frac"] < 0.5 and pal["ring0_corot"]["tf"] < 0.5   # mostly along z in both
        rows.append((t, pal, cv, mod["corot"]))
    print("\nCV of |B_t| in the ROI, % (model / Palace)")
    for t, pal, cv, co in rows:
        cells = " ".join(f"k{k} {cv[k] * 100:5.1f}/{pal[f'k{k}']['cv'] * 100:5.1f}" for k in range(1, 9))
        extra = "" if t == "lp" else f"  ring0 {cv[0] * 100:.1f}/{pal['ring0_counter']['cv'] * 100:.1f}  corot {co['bt_cv'] * 100:.1f}/{pal['ring0_corot']['cv'] * 100:.1f}"
        print(f"  {t}: {cells}{extra}")


@pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")
def test_page_palace_table_matches_results():
    """The page's "Measured in Palace" table (PALACE_U) is the results/*.json numbers, rounded to 0.1 %."""
    page = run_js([])["palace"]
    for t, name in REFS.items():
        pal = palace_table(name)
        assert set(page[t]) == set(pal), t
        for key, v in pal.items():
            assert page[t][key]["cv"] == pytest.approx(v["cv"] * 100, abs=0.051), (t, key)
            assert page[t][key]["tf"] == pytest.approx(v["tf"], abs=0.0051), (t, key)
            assert page[t][key]["f"] == pytest.approx(v["f_mhz"][0], abs=0.051), (t, key)
