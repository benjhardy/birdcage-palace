"""CLI checks that need no Palace: spec files, the Palace launcher option, a small coarse build."""
import json

import pytest

from birdcage_palace import cli
from birdcage_palace.design import CoilDesign


def _spec(tmp_path, d):
    p = tmp_path / "spec.json"
    p.write_text(json.dumps(d))
    return p


def test_spec_sets_build_options_and_flags_override(tmp_path):
    p = _spec(tmp_path, {"schema": "birdcage-palace/spec-1",
                         "build": {"legs": 12, "diameter": 0.24, "type": "highpass", "leg-width": 0.01, "out": "runs/x"}})
    a = cli.parse_args(["build", "--spec", str(p)])
    assert (a.legs, a.diameter, a.type, a.leg_width, a.out) == (12, 0.24, "highpass", 0.01, "runs/x")
    assert a.length == 0.25  # untouched options keep the CLI defaults
    a = cli.parse_args(["build", "--spec", str(p), "--legs", "8", "--out", "runs/y"])
    assert (a.legs, a.out, a.type) == (8, "runs/y", "highpass")


def test_flat_spec_and_errors(tmp_path):
    a = cli.parse_args(["build", "--spec", str(_spec(tmp_path, {"legs": 20, "phantom": True})), "--out", "o"])
    assert a.legs == 20 and a.phantom is True
    with pytest.raises(SystemExit):
        cli.parse_args(["build", "--spec", str(_spec(tmp_path, {"rungs": 12})), "--out", "o"])
    with pytest.raises(SystemExit):
        cli.parse_args(["build", "--spec", str(_spec(tmp_path, {"type": "notch"})), "--out", "o"])
    with pytest.raises(SystemExit):
        cli.parse_args(["build", "--legs", "12"])  # no --out anywhere


def test_palace_cmd_option_and_env(monkeypatch):
    monkeypatch.delenv("BIRDCAGE_PALACE_CMD", raising=False)
    assert cli.parse_args(["solve", "--case", "c"]).palace_cmd == "run_palace {cfg} {np}"
    monkeypatch.setenv("BIRDCAGE_PALACE_CMD", "palace -np {np} {cfg}")
    for sub in ("solve", "tune"):
        assert cli.parse_args([sub, "--case", "c"]).palace_cmd == "palace -np {np} {cfg}"
    assert cli.parse_args(["tune", "--case", "c", "--palace-cmd", "x {cfg} {np}"]).palace_cmd == "x {cfg} {np}"


def test_tune_summary_fields():
    des = CoilDesign(coil_type="highpass", c_ring=37.8e-12, notes={"c_leg_first_guess": None, "c_ring_first_guess": 36.9e-12})
    res = {"homogeneous_pair": {"f_mean_mhz": 127.73, "split_mhz": 0.004, "roi_bt_cv": [0.034, 0.034],
                                "b1plus_roi_cv": 0.028, "b1minus_over_b1plus_roi": 0.023}}
    hist = [{"f_k1_mhz": 135.2, "solve": {"wall_s": 60}}, {"f_k1_mhz": 127.73, "solve": {"wall_s": 58}}]
    s = cli.tune_summary(des, res, hist)
    assert s["c_ring_pf"] == pytest.approx(37.8) and s["c_leg_pf"] is None
    assert s["f_k1_at_first_guess_mhz"] == 135.2 and s["solves"] == 2
    assert abs(s["rel_err"]) < 1e-3


def test_coarse_build_from_spec(tmp_path):
    pytest.importorskip("gmsh")
    out = tmp_path / "case"
    p = _spec(tmp_path, {"build": {"legs": 8, "diameter": 0.2, "length": 0.16, "type": "highpass", "shield": 0.28,
                                   "h_coil": 0.02, "h_gap": 0.01, "h_max": 0.06, "order": 1, "modes": 6}})
    cli.main(["build", "--spec", str(p), "--out", str(out)])
    cfg = json.loads((out / "config.json").read_text())
    an = json.loads((out / "analytic.json").read_text())
    ports = cfg["Boundaries"]["LumpedPort"]
    assert len(ports) == 16  # 8 gaps in each end ring
    assert all(pt["C"] == pytest.approx(an["c_ring"]) for pt in ports)
    assert cfg["Solver"]["Order"] == 1 and cfg["Solver"]["Eigenmode"]["N"] == 6
    assert cfg["Solver"]["Linear"]["Type"] == "MUMPS" and cfg["Solver"]["Linear"]["MGMaxLevels"] == 1
    des = CoilDesign.from_json(out / "design.json")
    assert des.legs == 8 and des.coil_type == "highpass"
