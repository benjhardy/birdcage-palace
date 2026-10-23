"""Port tuning without Palace: how circuitRF is found and called (native Linux or WSL), Touchstone I/O, the N-rung
layout, the netlist, and the resonance algebra on a synthetic lossless N-port."""
import json
import os
import platform
import stat

import numpy as np
import pytest

from birdcage_palace import cli, palace, porttune as pt, sparams


# ------------------------------------------------------------------------------------------- finding circuitRF
def _fake_crf(path, body='print("1.0.9")'):
    path.write_text("#!/usr/bin/env python3\nimport sys, json\n" + body + "\n")
    path.chmod(path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return path


def test_linux_runs_natively_without_wsl(tmp_path, monkeypatch):
    exe = _fake_crf(tmp_path / "circuitRF")
    monkeypatch.delenv("BIRDCAGE_CIRCUITRF", raising=False)
    c = pt.CircuitRF(crf=str(exe), system="Linux", environ={})
    argv = c.argv(["opt", "x.cnl", "--json"], tmp_path)
    assert argv == [str(exe), "opt", "x.cnl", "--json"]
    assert not any("wsl" in a for a in argv)


def test_crf_option_env_and_auto(tmp_path, monkeypatch):
    a, b = _fake_crf(tmp_path / "a"), _fake_crf(tmp_path / "b")
    assert pt.CircuitRF(system="Linux", environ={"BIRDCAGE_CIRCUITRF": str(a)}).native_exe() == str(a)
    assert pt.CircuitRF(crf=str(b), system="Linux", environ={"BIRDCAGE_CIRCUITRF": str(a)}).native_exe() == str(b)
    # auto: circuitRF / circuitrf on the PATH, else ~/.local/share/circuitRF/current/circuitRF
    bindir = tmp_path / "bin"
    bindir.mkdir()
    lower = _fake_crf(bindir / ("circuitrf.exe" if platform.system() == "Windows" else "circuitrf"))
    monkeypatch.setenv("PATH", str(bindir))
    assert os.path.samefile(pt.CircuitRF(system="Linux", environ={}).native_exe(), lower)
    lower.unlink()
    home = tmp_path / "home"
    inst = home / ".local" / "share" / "circuitRF" / "current"
    inst.mkdir(parents=True)
    _fake_crf(inst / "circuitRF")
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("USERPROFILE", str(home))
    assert os.path.samefile(pt.CircuitRF(system="Linux", environ={}).native_exe(), inst / "circuitRF")
    (inst / "circuitRF").unlink()
    with pytest.raises(FileNotFoundError, match="circuitrf-tuning"):
        pt.CircuitRF(system="Linux", environ={}).native_exe()
    with pytest.raises(FileNotFoundError):
        pt.CircuitRF(crf=str(tmp_path / "nope"), system="Linux", environ={}).native_exe()


def test_windows_goes_through_wsl(tmp_path):
    auto = pt.CircuitRF(system="Windows", environ={}).argv(["opt", "t.cnl", "--json"], tmp_path)
    assert auto[:4] == ["wsl.exe", "--exec", "bash", "-lc"]  # auto: the WSL default distribution, no -d
    assert "command -v circuitRF" in auto[-1] and auto[-1].endswith(" opt t.cnl --json")
    named = pt.CircuitRF(distro="Ubuntu-24.04", system="Windows", environ={}).argv(["--version"], tmp_path)
    assert named[:3] == ["wsl.exe", "-d", "Ubuntu-24.04"]
    env = pt.CircuitRF(system="Windows", environ={"BIRDCAGE_WSL_DISTRO": "Debian", "BIRDCAGE_CIRCUITRF": "~/crf/circuitRF"})
    argv = env.argv(["--version"], tmp_path)
    assert argv[:3] == ["wsl.exe", "-d", "Debian"] and '"$HOME"/crf/circuitRF' in argv[-1]


def test_palace_wsl_prefix():
    assert palace.wsl_prefix("auto") == ["wsl"] and palace.wsl_prefix(None) == ["wsl"]
    assert palace.wsl_prefix("Ubuntu-24.04") == ["wsl", "-d", "Ubuntu-24.04"]


@pytest.mark.skipif(platform.system() == "Windows", reason="native run: Linux/macOS only")
def test_native_opt_and_sparam_run(tmp_path):
    """A stand-in circuitRF on Linux: the JSON on stdout is parsed and the -o file is where tune expects it."""
    body = ('a = sys.argv[1:]\n'
            'if a[0] == "opt":\n'
            '    print(json.dumps({"result": {"optimize": {"variables": [{"key": "Ct", "best": "2.5 pF"}]}}}))\n'
            'elif a[0] == "sparam":\n'
            '    open(a[a.index("-o") + 1], "w").write("# MHz S RI R 50\\n1 0 0 1 0 0 0 0 0\\n")\n')
    exe = _fake_crf(tmp_path / "circuitRF", body)
    c = pt.CircuitRF(crf=str(exe), environ={})
    cnl = tmp_path / "t.cnl"
    cnl.write_text("; x\n")
    o = pt.crf_opt(tmp_path, cnl, c)
    assert pt.best_vals(o, {"Ct": 1e-12})["Ct"] == pytest.approx(2.5e-12)
    out = pt.crf_sparam(tmp_path, cnl, "o.s2p", crf=c)
    f, S, z = pt.read_touchstone(out)
    assert z == 50 and S[0, 1, 0] == 1 and S[0, 0, 1] == 0  # 2-port order S11 S21 S12 S22


# ------------------------------------------------------------------------------------------- CLI
def test_cli_porttune_options(monkeypatch):
    monkeypatch.setenv("BIRDCAGE_WSL_DISTRO", "Ubuntu-24.04")
    a = cli.parse_args(["porttune", "tune", "--snp", "x.s12p", "--name", "u", "--crf", "/opt/crf/circuitRF"])
    assert (a.pt_cmd, a.crf, a.wsl_distro, a.esr, a.band_khz) == ("tune", "/opt/crf/circuitRF", None, 0.68, 2.0)
    c = pt.CircuitRF(a.crf, a.wsl_distro, system="Windows")
    assert c.distro == "Ubuntu-24.04"  # env used when the flag is absent
    a = cli.parse_args(["porttune", "tune", "--snp", "x.s12p", "--name", "u", "--wsl-distro", "Other", "--c-ref", "3e-12"])
    assert a.wsl_distro == "Other" and a.c_eigen == 3e-12
    for sub in ("assemble --cases p --out o.s12p", "resonance --snp o.s12p", "verify --snp o --ref-case a --tuned-case b"):
        assert cli.parse_args(["porttune"] + sub.split()).fn
    monkeypatch.delenv("BIRDCAGE_WSL_DISTRO")
    assert cli.parse_args(["solve", "--case", "c"]).wsl_distro == "auto"


# ------------------------------------------------------------------------------------------- Touchstone, layout, netlist
def test_touchstone_round_trip_12_port(tmp_path):
    rng = np.random.default_rng(1)
    S = rng.normal(size=(5, 12, 12)) + 1j * rng.normal(size=(5, 12, 12))
    f = np.linspace(100e6, 160e6, 5)
    p = sparams.write_touchstone(tmp_path / "x.s12p", f, S, 500.0, ["test"])
    f2, S2, r = pt.read_touchstone(p)
    assert r == 500 and np.allclose(f2, f) and np.abs(S2 - S).max() < 1e-9


def test_layout_for_n_rungs():
    assert pt.layout(16) == ((0, 4), (2, 10))
    assert pt.layout(12) == ((0, 3), (1, 2, 7, 8))
    assert pt.layout(8) == ((0, 2), (1, 5))
    with pytest.raises(ValueError):
        pt.layout(10)
    assert pt.range_scale(1.825e-12) == 1.0 and pt.range_scale(3.1e-12) == 1.7


def test_netlist_12_rungs(tmp_path):
    snp = sparams.write_touchstone(tmp_path / "c.s12p", np.array([1e8, 2e8]), np.zeros((2, 12, 12)), 500.0)
    vals = {"Ct": 3e-12, "Ci": 3e-12, "Cd1": 3.5e-12, "Cd2": 3.5e-12, "Cp1": 30e-12, "Cp2": 30e-12, "Cm1": 30e-12, "Cm2": 30e-12}
    drive, iso = pt.layout(12)
    (tmp_path / "w").mkdir()
    cnl = pt.write_cnl(tmp_path / "w", snp, "u", 127.73e6, 0.68, vals, ["Ct", "Ci", "Cd1", "Cp1", "Cm1"], drive=drive,
                       iso=iso, n_ports=12, band_khz=2, tie={"Cd2": "Cd1", "Cp2": "Cp1", "Cm2": "Cm1"}, scale=1.7)
    t = cnl.read_text()
    assert 'SnP:COIL g0 g1 g2 g3 g4 g5 g6 g7 g8 g9 g10 g11 NumPorts=12 File="c.s12p"' in t
    assert (tmp_path / "w" / "c.s12p").exists()
    for i in iso:
        assert f"C:C{i} g{i} e{i} C=Ci" in t
    assert "C:C3 g3 e3 C=Cd2" in t and "C:CM2 a2 p2 C=Cm2" in t and "Cd2 = Cd1" in t
    assert "tune Ct min=1.7 pF max=5.1 pF scale=lin opt=1" in t
    assert "tune Cm1 min=0.05 pF max=300 pF scale=log opt=1" in t  # the series capacitor keeps its absolute range
    assert t.count("goal ") == 3 and "optimize algorithm=minimax" in t
    caps = pt.rung_caps(vals, 12, drive, iso)
    assert caps[0] == pytest.approx(1 / (1 / 3.5e-12 + 1 / 30e-12)) and caps[5] == 3e-12


# ------------------------------------------------------------------------------------------- resonance algebra
def _circulant_coil(N=12, f=np.linspace(100e6, 160e6, 31), r=500.0):
    """Lossless synthetic coil: gap admittance Y = G / (j w) with a symmetric circulant G (1/H), as a ladder of
    inductors would give. Harmonic k has G eigenvalue g_k, so with C on every gap k resonates at w^2 = g_k / C."""
    n = np.arange(N)
    # g_k = 2.4e7 - 1.8e7 cos(2 pi k / N) + 0.05e7 (-1)^k, plus 5e7 on k = 0 (no closed rung loop): k = 1 is lowest
    col = (np.where(n == 0, 2.4e7, 0.0) - np.where((n == 1) | (n == N - 1), 0.9e7, 0.0) + np.where(n == N // 2, 0.05e7, 0.0)
           + 5e7 / N)
    G = np.array([np.roll(col, i) for i in range(N)])
    I = np.eye(N)
    S = np.array([np.linalg.solve(I + r * G / (1j * 2 * np.pi * fk), I - r * G / (1j * 2 * np.pi * fk)) for fk in f])
    g = np.real(np.fft.fft(col))
    return f, S, r, g


def test_resonance_and_c_for_f0_on_synthetic_12_port():
    f, S, r, g = _circulant_coil()
    f0 = 3.0 * pt.GAMMA_BAR
    w0 = 2 * np.pi * f0
    c_exact = g[1] / w0 ** 2
    assert pt.c_k1_estimate(f, S, r, f0) == pytest.approx(c_exact, rel=1e-6)
    c, fr = pt.c_for_frequency(f, S, r, f0)
    assert c == pytest.approx(c_exact, rel=1e-5) and fr == pytest.approx(f0, rel=1e-6)
    k1, rs = pt.k1_frequency(f, S, r, 1.1 * c_exact)
    assert k1 == pytest.approx(f0 / np.sqrt(1.1), rel=1e-5)
    # a lossless N-port: S S^H = I
    assert np.abs(np.einsum("fij,fkj->fik", S, S.conj()) - np.eye(12)).max() < 1e-10
    # numpy termination: with every gap at c_exact and the drive gaps fed, the network is solvable at f0
    vals = {"Ct": c_exact, "Ci": c_exact, "Cd1": 1.2 * c_exact, "Cd2": 1.2 * c_exact, "Cp1": 6 * c_exact, "Cp2": 6 * c_exact,
            "Cm1": 30e-12, "Cm2": 30e-12}
    drive, iso = pt.layout(12)
    S2 = pt.python_two_port(np.array([f0]), (f, S), r, vals, 0.68, drive, iso)
    assert S2.shape == (1, 2, 2) and np.all(np.abs(S2) <= 1 + 1e-9)


def test_portcase_sweep_scales_with_b0(tmp_path):
    pytest.importorskip("gmsh")
    base = tmp_path / "b"
    cli.main(["build", "--legs", "8", "--diameter", "0.2", "--length", "0.16", "--shield", "0.28", "--b0", "1.5",
              "--h-coil", "0.02", "--h-gap", "0.01", "--h-max", "0.06", "--order", "1", "--modes", "2",
              "--symmetry", "eighth", "--sym-bc", "x=pec,y=pmc,z=pec", "--out", str(base)])
    cli.main(["portcase", "--case", str(base), "--sym-bc", "x=pmc,y=pec,z=pec", "--out", str(tmp_path / "p_me")])
    cfg = json.loads((tmp_path / "p_me" / "config.json").read_text())
    smp = cfg["Solver"]["Driven"]["Samples"][0]
    assert smp["MinFreq"] == pytest.approx(0.05) and smp["MaxFreq"] == pytest.approx(0.08) and smp["FreqStep"] == pytest.approx(0.001)
    assert cfg["Solver"]["Linear"]["Type"] == "MUMPS" and cfg["Solver"]["Linear"]["MGMaxLevels"] == 1
    assert len(cfg["Boundaries"]["LumpedPort"]) == 2  # 8 rungs / 4 per quarter


# ------------------------------------------------------------------------------------------- the web builder's package
PAGE_JS = r"""
const fs = require('fs'), vm = require('vm');
const html = fs.readFileSync(process.argv[process.argv.length - 1], 'utf8');
const core = html.match(/<script id="core">([\s\S]*?)<\/script>/)[1];
const ctx = {}; vm.createContext(ctx);
vm.runInContext(core + '\nthis.Ladder = Ladder; this.Spec = Spec; this.Pkg = Pkg; this.Cost = Cost;', ctx);
const base = Object.assign({}, ctx.Spec.DEFAULTS, { legs: 12, diameter: 0.24, mesh: 'coarse', h_coil: 0.008, h_gap: 0.004, h_max: 0.04 });
const out = { reasons: {} };
for (const [k, d] of Object.entries({ ok: {}, highpass: { type: 'highpass' }, legs10: { legs: 10 }, noshield: { shield_mode: 'without' },
    both: { shield_mode: 'both' }, phantom: { phantom: true }, open: { shield_ends: 'absorbing' } }))
  out.reasons[k] = ctx.Spec.portReason(Object.assign({}, base, d));
for (const tm of ['eigen', 'port']) {
  const S = Object.assign({}, base, { tune_method: tm }); S.shielded = true;
  const lad = ctx.Ladder.design(ctx.Spec.toDesign(S)), est = ctx.Cost.estimate(S);
  const f = ctx.Pkg.files(S, lad, est, '2026-10-23');
  out[tm] = { files: Object.keys(f), runsh: f['run.sh'], readme: f['README.md'], spec: JSON.parse(f['spec.json']),
    commands: ctx.Pkg.commands(S, 'birdcage'), dir: ctx.Pkg.dirName(S), cost: ctx.Cost.total(S) };
}
// a port choice on an ineligible design falls back to the eigenmode loop
const H = Object.assign({}, base, { type: 'highpass', tune_method: 'port' }); H.shielded = true;
out.hp_port_mode = ctx.Spec.portMode(H);
process.stdout.write(JSON.stringify(out));
"""


@pytest.fixture(scope="module")
def page():
    import shutil
    import subprocess
    from pathlib import Path
    node = shutil.which("node")
    if node is None:
        pytest.skip("node not installed")
    html = Path(__file__).resolve().parents[1] / "web" / "birdcage-builder" / "index.html"
    r = subprocess.run([node, "-e", PAGE_JS, str(html)], capture_output=True, text=True, check=True)
    return json.loads(r.stdout)


def test_page_port_option_only_for_eligible_designs(page):
    rs = page["reasons"]
    assert rs["ok"] == ""
    for k in ("highpass", "legs10", "noshield", "both", "phantom", "open"):
        assert rs[k], k
    assert "low-pass" in rs["highpass"] and "multiple of 4" in rs["legs10"]
    assert page["hp_port_mode"] is False
    assert page["eigen"]["dir"] == "birdcage_lp12_d24_3T_coarse" and page["port"]["dir"] == "birdcage_lp12_d24_3T_coarse_pt"


def test_page_port_package_commands_parse(page):
    """Every birdcage command the port-tuning package writes is a valid CLI call, in the guide's order."""
    p = page["port"]
    cmds = p["commands"]
    subs = [c.split()[1] + (" " + c.split()[2] if c.split()[1] == "porttune" else "") for c in cmds]
    assert subs == ["build", "portcase", "portcase", "portcase", "portcase", "solve", "solve", "solve", "solve",
                    "porttune assemble", "porttune resonance", "porttune tune", "build", "variant",
                    "solve", "solve", "porttune verify"]
    bcs = [c.split("--sym-bc ")[1].split()[0] for c in cmds if " portcase " in c]
    assert bcs == ["x=pec,y=pmc,z=pec", "x=pmc,y=pec,z=pec", "x=pmc,y=pmc,z=pec", "x=pec,y=pec,z=pec"]
    for c in cmds:
        for part in c.split(" && "):
            argv = part.replace("$NP", "4").split()[1:]
            a = cli.parse_args(argv)
            assert a.fn
            if argv[0] == "build":
                assert a.symmetry in ("eighth", "zhalf") and a.legs == 12
    assert "birdcage12.s12p" in " ".join(cmds)
    rs = p["runsh"]
    for s in ("porttune assemble", "porttune tune", "INSTALL_CIRCUITRF", "potatobeanradio/circuitRF/releases/download/1.0.9",
              "sha256sum -c", "STOP_AFTER", "SKIP_VERIFY", "circuitrf-tuning.md"):
        assert s in rs, s
    assert "MIT" in p["readme"] and "--help" in p["readme"]
    assert p["spec"]["tune_method"] == "circuitrf-port" and p["spec"]["build"]["symmetry"] == "eighth"
    assert page["eigen"]["spec"]["tune_method"] == "eigenmode-loop" and "porttune" not in page["eigen"]["runsh"]
    assert 0 < p["cost"]["lo"] < p["cost"]["hi"] < 3600
