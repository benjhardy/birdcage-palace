"""Shield choice (with / without / both), the unshielded build, compare, and the load-study config options. No Palace."""
import json
import math
import shutil
import subprocess
from pathlib import Path

import pytest

from birdcage_palace import cli
from birdcage_palace.design import CoilDesign, Phantom
from birdcage_palace.palace import EPS0, materials

PAGE = Path(__file__).resolve().parents[1] / "web" / "birdcage-builder" / "index.html"

JS = r"""
const fs = require('fs'), vm = require('vm');
const html = fs.readFileSync(process.argv[process.argv.length - 1], 'utf8');
const core = html.match(/<script id="core">([\s\S]*?)<\/script>/)[1];
const ctx = {}; vm.createContext(ctx);
vm.runInContext(core + '\nthis.Ladder = Ladder; this.Spec = Spec; this.Pkg = Pkg; this.Cost = Cost;', ctx);
const out = {};
for (const mode of ['with', 'without', 'both']) for (const ph of [false, true]) {
  const S = Object.assign({}, ctx.Spec.DEFAULTS, { shield_mode: mode, phantom: ph, mesh: 'coarse', h_coil: 0.008, h_gap: 0.004, h_max: 0.04 });
  S.shielded = mode !== 'without';
  const lad = ctx.Ladder.design(ctx.Spec.toDesign(S)), est = ctx.Cost.estimate(S);
  const cases = ctx.Spec.cases(S).map(c => ({ shielded: c.shielded, name: ctx.Spec.caseName(c), args: ctx.Spec.buildArgs(c) }));
  out[mode + (ph ? '_ph' : '')] = { cases, files: Object.keys(ctx.Pkg.files(S, lad, est, '2026-10-21')),
    runsh: ctx.Pkg.files(S, lad, est, '2026-10-21')['run.sh'], plan: ctx.Pkg.plan(S) };
}
process.stdout.write(JSON.stringify(out));
"""


@pytest.fixture(scope="module")
def js():
    node = shutil.which("node")
    if node is None:
        pytest.skip("node not installed")
    r = subprocess.run([node, "-e", JS, str(PAGE)], capture_output=True, text=True, check=True)
    return json.loads(r.stdout)


def _argv(args, out):
    v = ["build"]
    for k, x in args:
        v += [f"--{k}"] if x is None else [f"--{k}", str(x)]
    return v + ["--out", out]


def test_page_shield_modes_give_cli_cases(js):
    w, wo, b = js["with"], js["without"], js["both"]
    assert [c["shielded"] for c in w["cases"]] == [True]
    assert [c["shielded"] for c in wo["cases"]] == [False]
    assert [c["shielded"] for c in b["cases"]] == [True, False]
    shields = [cli.parse_args(_argv(c["args"], "o")).shield for c in b["cases"]]
    assert shields == [0.34, 0]
    assert [cli.parse_args(_argv(c["args"], "o")).modes for c in b["cases"]] == [20, 2]   # unshielded: k = 1 pair only
    assert b["cases"][1]["name"].endswith("_noshield") and not b["cases"][0]["name"].endswith("_noshield")
    assert "spec_noshield.json" in b["files"] and "spec_noshield.json" not in w["files"]
    assert b["plan"]["compare"] and "compare --cases" in b["runsh"] and w["plan"]["compare"] is None
    assert b["runsh"].count(" tune  --case ") == 2


def test_page_phantom_tunes_air_then_loads(js):
    p = js["with_ph"]
    a = cli.parse_args(_argv(p["cases"][0]["args"], "o"))
    assert a.phantom and a.phantom_eps == 1 and a.phantom_sigma == 0   # tuned with the region filled with air
    blk = p["plan"]["blocks"][0]
    assert "load  --case" in blk["load"] and "--phantom-eps 78" in blk["load"] and "--esr 0.68" in blk["load"]
    lv = blk["load"].replace('"$NP"', "4").split()[1:]
    la = cli.parse_args(lv)
    assert la.cmd == "load" and la.phantom_sigma == 0.6 and la.copper_sigma == 5.8e7
    assert len(js["both_ph"]["plan"]["blocks"]) == 2 and all(b["load"] for b in js["both_ph"]["plan"]["blocks"])


def test_unshielded_build_has_absorbing_box(tmp_path):
    pytest.importorskip("gmsh")
    out = tmp_path / "u"
    cli.main(["build", "--legs", "8", "--diameter", "0.2", "--length", "0.16", "--shield", "0", "--h-coil", "0.02",
              "--h-gap", "0.01", "--h-max", "0.07", "--order", "1", "--modes", "6", "--out", str(out)])
    cfg = json.loads((out / "config.json").read_text())
    info = json.loads((out / "build.json").read_text())
    assert info["attrs"]["shield"] is None and cfg["Boundaries"]["Absorbing"]["Attributes"] == [12]
    assert cfg["Boundaries"]["PEC"]["Attributes"] == [10]
    assert cfg["Solver"]["Eigenmode"]["Tol"] == 1e-6 and cfg["Solver"]["Eigenmode"]["N"] == 2   # k = 1 pair only
    gap = max(0.15, 299792458 / (42.577478e6 * 3.0) / 4)   # lambda0 / 4 = 0.587 m at 3 T
    assert info["domain_radius"] == pytest.approx(0.1 + gap) and info["domain_length"] == pytest.approx(0.16 + 2 * gap)
    assert cfg["Boundaries"]["Absorbing"]["Order"] == 2
    an = json.loads((out / "analytic.json").read_text())
    assert an["c_leg"] > 0
    # a larger air box on request
    out2 = tmp_path / "u2"
    cli.main(["build", "--legs", "8", "--diameter", "0.2", "--length", "0.16", "--shield", "0", "--air-radius", "0.4",
              "--air-length", "1.0", "--h-coil", "0.02", "--h-gap", "0.01", "--h-max", "0.08", "--order", "1",
              "--out", str(out2)])
    i2 = json.loads((out2 / "build.json").read_text())
    assert (i2["domain_radius"], i2["domain_length"]) == (0.4, 1.0)


def test_losstan_and_coil_loss_config(tmp_path):
    des = CoilDesign(phantom=Phantom(eps_r=78, sigma=0.6, loss_model="losstan"))
    m = materials(des, {"air": 1, "phantom": 2})[1]
    assert m["LossTan"] == pytest.approx(0.6 / (2 * math.pi * des.f_larmor * EPS0 * 78))
    assert "Conductivity" not in m
    des.phantom.loss_model = "conductivity"
    assert materials(des, {"air": 1, "phantom": 2})[1]["Conductivity"] == 0.6
    des.phantom.sigma = 0
    assert set(materials(des, {"air": 1, "phantom": 2})[1]) == {"Attributes", "Permittivity", "Permeability"}
    pytest.importorskip("gmsh")
    out = tmp_path / "c"
    cli.main(["build", "--legs", "8", "--diameter", "0.2", "--length", "0.16", "--shield", "0.28", "--phantom",
              "--phantom-radius", "0.05", "--phantom-length", "0.1", "--phantom-eps", "1", "--phantom-sigma", "0",
              "--h-coil", "0.02", "--h-gap", "0.01", "--h-max", "0.06", "--order", "1", "--modes", "2", "--out", str(out)])
    v = tmp_path / "v"
    cli.main(["variant", "--case", str(out), "--out", str(v), "--phantom-eps", "50", "--phantom-sigma", "0.3",
              "--copper-sigma", "5.8e7", "--esr", "0.68"])
    cfg = json.loads((v / "config.json").read_text())
    assert cfg["Domains"]["Materials"][1]["Permittivity"] == 50 and cfg["Domains"]["Materials"][1]["Conductivity"] == 0.3
    cond = cfg["Boundaries"]["Conductivity"]
    assert cond[0]["Attributes"] == [10] and cond[1] == {"Attributes": [11], "Conductivity": 5.8e7, "External": True}
    assert 10 not in cfg["Boundaries"]["PEC"]["Attributes"]
    des = CoilDesign.from_json(v / "design.json")
    lp = cfg["Boundaries"]["LumpedPort"][0]
    w0c = 2 * math.pi * des.f_larmor * des.c_leg
    assert lp["R"] == pytest.approx(0.68 * (1 + (1 / (w0c * 0.68)) ** 2))   # series ESR -> parallel R at f0


def test_compare_two_cases(tmp_path):
    def case(name, shield, hist, cv):
        d = tmp_path / name
        d.mkdir()
        CoilDesign(shield_diameter=shield, c_leg=hist[-1][0]).to_json(d / "design.json")
        (d / "tune_history.json").write_text(json.dumps([{"c_leg": c, "c_ring": None, "f_k1_mhz": f} for c, f in hist]))
        (d / "tune_summary.json").write_text(json.dumps({"c_leg_pf": hist[-1][0] * 1e12, "f_k1_mhz": hist[-1][1],
                                                          "c_leg_first_guess_pf": hist[0][0] * 1e12,
                                                          "f_k1_at_first_guess_mhz": hist[0][1], "b1plus_roi_cv": cv,
                                                          "solves": len(hist)}))
        (d / "modes.json").write_text(json.dumps({"homogeneous_pair": {"m": [1, 2]},
                                                   "modes": [{"m": 1, "Q": 1e15 if shield else 900.0},
                                                             {"m": 2, "Q": 1e15 if shield else 910.0}]}))
        return d
    a = case("s", 0.34, [(2.8e-12, 111.0), (1.83e-12, 127.73)], 0.023)
    # unshielded: 1/f^2 = L (C + Cs) with L, Cs chosen so the fit is exact
    L, Cs = 4.0e-4 / 1e-12, 1.0e-12
    f = lambda c: (L * (c + Cs)) ** -0.5  # noqa: E731
    b = case("u", None, [(2.5e-12, f(2.5e-12)), (2.2e-12, f(2.2e-12))], 0.025)
    out = cli.cmd_compare(cli.parse_args(["compare", "--cases", str(a), str(b), "--out", str(tmp_path / "cmp.json")]))
    assert out["k1_shift"]["f_case2_at_same_c_mhz"] == pytest.approx(f(1.83e-12), rel=1e-9)
    assert out["cases"][0]["Q_k1"] is None and out["cases"][1]["Q_k1"] == pytest.approx(905)
    assert (tmp_path / "cmp.md").exists()
