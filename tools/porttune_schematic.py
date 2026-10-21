"""Schematic of the circuitRF port-tuning network (runs/pt_e/circuitrf/tune_unloaded_after.cnl), drawn with schemdraw.

Values are parsed from the netlist itself; the optimiser goal results come from results/lp16_porttune.json.

    python tools/porttune_schematic.py <netlist.cnl> <porttune.json> <out.png>
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import schemdraw  # noqa: E402
import schemdraw.elements as elm  # noqa: E402

INK, INK2, BLUE, ORANGE = "#0b0b0b", "#52514e", "#2a78d6", "#eb6834"


def parse(cnl: Path) -> tuple[dict, dict]:
    vals, caps = {}, {}
    for ln in cnl.read_text().splitlines():
        m = re.match(r"^(\w+)\s*=\s*([\d.eE+-]+)\s*(pF|Ohm)", ln)
        if m:
            vals[m.group(1)] = float(m.group(2))
        m = re.match(r"^C:C(\d+)\s+g\d+\s+\S+\s+C=(\w+)", ln)
        if m:
            caps[int(m.group(1))] = m.group(2)
    return vals, caps


def main(cnl: Path, pj: Path, out: Path):
    vals, caps = parse(cnl)
    goals = {g["name"]: g for g in json.loads(pj.read_text())["unloaded"]["raw"]["goals"]}
    esr = f"{vals['Resr']:g} Ω"
    fig = plt.figure(figsize=(16, 8.4), dpi=150, facecolor="white")
    ax = fig.add_axes([0.0, 0.065, 1.0, 0.79])
    ax.axis("off")
    d = schemdraw.Drawing(canvas=ax, show=False)
    d.config(fontsize=10.5, unit=2.0, lw=1.3)
    top = [9, 10, 11, 12, 13, 14, 15]
    bot = [1, 2, 3, 5, 6, 7, 8]
    sp = 2.3
    W = sp * 7 + 0.6
    pins = [elm.IcPin(name="0", anchorname="g0", side="left"), elm.IcPin(name="4", anchorname="g4", side="right")]
    pins += [elm.IcPin(name=str(n), anchorname=f"g{n}", side="bottom", slot=f"{i + 1}/7") for i, n in enumerate(bot)]
    pins += [elm.IcPin(name=str(n), anchorname=f"g{n}", side="top", slot=f"{i + 1}/7") for i, n in enumerate(top)]
    ic = d.add(elm.Ic(pins=pins, size=(W, 3.2), pinspacing=sp, edgepadW=0.3, plblsize=9, lofst=0.12)
               .at((0, 0)).color(BLUE))
    d.add(elm.Label().at((W / 2, 1.6)).label("COIL: Palace 16-port S-matrix\n(bare low-pass coil, one port per gap)",
                                               fontsize=11, color=BLUE))

    def branch(anchor, n, up):
        c = caps[n]
        s = 1 if up else -1
        d.add(elm.Line().at(anchor).toy(anchor[1] + s * 0.5))
        lab = f"C{n}: {c}\n{vals[c]:.3f} pF"
        col = ORANGE if c == "Ci" else INK
        e = elm.Capacitor().up() if up else elm.Capacitor().down()
        d.add(e.label(lab, loc="top" if up else "bottom", fontsize=8.6, ofst=0.12, color=col).color(col))
        r = elm.Resistor().up() if up else elm.Resistor().down()
        d.add(r.label(esr, loc="top" if up else "bottom", fontsize=8.2, ofst=0.12, color=INK2))
        d.add(elm.Ground().flip() if up else elm.Ground())

    for n in bot:
        branch(ic.absanchors[f"g{n}"], n, up=False)
    for n in top:
        branch(ic.absanchors[f"g{n}"], n, up=True)

    def drive(anchor, k, gap, left):
        go = (lambda e: e.left()) if left else (lambda e: e.right())
        d.add(go(elm.Line()).at(anchor).length(0.6))
        cd = f"Cd{k}"
        d.add(go(elm.Capacitor()).label(f"C{gap}: {cd}\n{vals[cd]:.3f} pF", loc="top", fontsize=8.6, ofst=0.15))
        d.add(go(elm.Resistor()).label(esr, loc="top", fontsize=8.2, color=INK2, ofst=0.15))
        node = d.add(elm.Dot())
        a = node.absanchors["center"]
        cp, cm = f"Cp{k}", f"Cm{k}"
        d.add(elm.Capacitor().down().at(a).label(f"{cp}\n{vals[cp]:.2f} pF", loc="bottom", fontsize=8.6, ofst=0.15,
                                                  color=BLUE))
        d.add(elm.Resistor().down().label(esr, loc="bottom", fontsize=8.2, color=INK2, ofst=0.15))
        d.add(elm.Ground())
        d.add(go(elm.Capacitor()).at(a).label(f"{cm}\n{vals[cm]:.2f} pF", loc="top", fontsize=8.6, ofst=0.15,
                                              color=BLUE))
        d.add(go(elm.Line()).length(0.6))
        d.add(elm.SourceSin().down().label(f"P{k}\n50 Ω", loc="top" if left else "bottom", fontsize=10, ofst=0.2))
        d.add(elm.Ground())

    drive(ic.absanchors["g0"], 1, 0, left=True)
    drive(ic.absanchors["g4"], 2, 4, left=False)
    d.draw(show=False)
    ax.set_aspect("equal")
    ax.autoscale_view()
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    ax.set_xlim(x0 - 0.3, x1 + 0.3)
    ax.set_ylim(y0 - 0.3, y1 + 0.3)

    fig.text(0.03, 0.955, "circuitRF tuning network around the Palace 16-port coil", fontsize=19, weight="bold",
             color=INK)
    fig.text(0.03, 0.915, "16-rung low-pass birdcage, 3 T, unloaded. Every capacitor gap is a port of the Palace "
             "S-matrix; circuitRF closes the gaps with capacitors (ESR 0.68 Ω each) and adds the drive network.",
             fontsize=10.5, color=INK2)
    fig.text(0.03, 0.885, "Gaps 0 and 4 (90° apart) are fed from 50 Ω ports P1 and P2 through a capacitive divider "
             "(Cp shunt, Cm series). Gaps 2 and 10 (orange) get their own value Ci for P1-P2 isolation.",
             fontsize=10.5, color=INK2)
    m1, m2, iso = goals["match1"], goals["match2"], goals["isol"]
    fig.text(0.03, 0.035, f"Tuned (circuitRF optimiser, all goals met), worst case over 127.73 MHz ± 2 kHz:  "
             f"|S11| {m1['value']:.1f} dB,  |S22| {m2['value']:.1f} dB,  |S21| {iso['value']:.1f} dB.  "
             f"Netlist: {cnl.name}", fontsize=10, color=INK2)
    fig.text(0.97, 0.035, "simulation only", fontsize=10, color=INK2, style="italic", ha="right")
    fig.savefig(out, dpi=150, facecolor="white")
    print("wrote", out)


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]))
