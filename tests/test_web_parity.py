"""The web builder's JavaScript ladder model (web/birdcage-builder/index.html, <script id="core">) against ladder.py.

Random designs across all three coil types, shielded and unshielded, with and without given capacitors. Every
inductance, first-guess capacitor and mode frequency must agree to 1e-6 relative.

    python -m pytest tests/test_web_parity.py -s      # -s prints the worst relative difference
"""
import json
import math
import random
import shutil
import subprocess
from pathlib import Path

import pytest

from birdcage_palace import ladder
from birdcage_palace.design import CoilDesign

PAGE = Path(__file__).resolve().parents[1] / "web" / "birdcage-builder" / "index.html"
TOL = 1e-6

JS = r"""
const fs = require('fs'), vm = require('vm');
const html = fs.readFileSync(process.argv[process.argv.length - 1], 'utf8');
const core = html.match(/<script id="core">([\s\S]*?)<\/script>/)[1];
const ctx = {}; vm.createContext(ctx); vm.runInContext(core + '\nthis.Ladder = Ladder;', ctx);
const designs = JSON.parse(fs.readFileSync(0, 'utf8'));
process.stdout.write(JSON.stringify(designs.map(d => ctx.Ladder.analytic(d))));
"""


def random_designs(n=24, seed=20261021):
    rng = random.Random(seed)
    types = ["lowpass", "highpass", "bandpass"]
    out = []
    for i in range(n):
        legs = rng.choice(range(8, 33, 2))
        dia = rng.uniform(0.12, 0.6)
        d = dict(legs=legs, diameter=dia, length=rng.uniform(0.6, 1.3) * dia, leg_width=rng.uniform(0.004, 0.02),
                 ring_width=rng.uniform(0.006, 0.025), coil_type=types[i % 3], b0=rng.choice([1.5, 3.0, 7.0, rng.uniform(0.5, 9.4)]),
                 shield_diameter=None if i % 5 == 4 else dia + rng.uniform(0.02, 0.14),
                 bandpass_ratio=rng.uniform(0.3, 3.0), c_leg=None, c_ring=None)
        if d["leg_width"] * legs >= 0.7 * math.pi * dia:
            d["leg_width"] = 0.5 * math.pi * dia / legs
        if i % 4 == 3:  # given capacitors instead of the first guess
            d["c_leg"] = rng.uniform(1, 60) * 1e-12
            d["c_ring"] = rng.uniform(5, 80) * 1e-12
        out.append(d)
    return out


def _flat(r):
    sp = r["spectrum"]
    v = {"c_leg": r["c_leg"], "c_ring": r["c_ring"], "f_corot_hz": sp.get("f_corot_hz"),
         "Lleg_self_nH": r["Lleg_self_nH"], "Lring_seg_self_nH": r["Lring_seg_self_nH"]}
    for i, f in enumerate(sp["f_hz"]):
        v[f"f_k{i}"] = f
    for i, (a, b) in enumerate(zip(r["Lleg_k_nH"], r["Lring_k_nH"])):
        v[f"Lleg_k{i}"], v[f"Lring_k{i}"] = a, b
    return v


def _rel(a, b):
    if a is None or b is None:
        return 0.0 if a is None and b is None else math.inf
    if (isinstance(a, float) and math.isnan(a)) or (isinstance(b, float) and math.isnan(b)):
        return 0.0 if (isinstance(a, float) and math.isnan(a)) and (isinstance(b, float) and math.isnan(b)) else math.inf
    scale = max(abs(a), abs(b))
    return 0.0 if scale == 0 else abs(a - b) / scale


def run_parity(designs):
    node = shutil.which("node")
    js = subprocess.run([node, "-e", JS, str(PAGE)], input=json.dumps(designs), capture_output=True, text=True, check=True)
    out_js = json.loads(js.stdout.replace("NaN", "null"))
    rows = []
    for d, rj in zip(designs, out_js):
        des = CoilDesign(**{k: v for k, v in d.items()})
        rp = json.loads(json.dumps(ladder.analytic(des)).replace("NaN", "null"))
        fp, fj = _flat(rp), _flat(rj)
        assert fp.keys() == fj.keys()
        worst = max(fp, key=lambda k: _rel(fp[k], fj[k]))
        rows.append((d, worst, _rel(fp[worst], fj[worst]), len(fp)))
    return rows


@pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")
def test_js_ladder_matches_python():
    designs = random_designs()
    assert {d["coil_type"] for d in designs} == {"lowpass", "highpass", "bandpass"}
    rows = run_parity(designs)
    worst = max(r[2] for r in rows)
    n_vals = sum(r[3] for r in rows)
    print(f"\n{len(rows)} designs, {n_vals} values compared, worst relative difference {worst:.2e}")
    for d, k, e, _ in rows:
        print(f"  {d['coil_type']:<9} N={d['legs']:>2} D={d['diameter']:.3f} shield={d['shield_diameter'] and round(d['shield_diameter'], 3)} "
              f"C given={d['c_leg'] is not None}: worst {k} {e:.1e}")
    assert worst <= TOL


if __name__ == "__main__":
    test_js_ladder_matches_python()
