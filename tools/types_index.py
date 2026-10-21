"""Contact sheet (index.html) for the low-pass / high-pass / band-pass birdcage visuals.

    python tools/types_index.py --out <dir> --case lowpass=runs/lp16_final:runs/lp16_final/tune_history_part1.json ...
        [--status "text"]

Every number is read from the cases' modes.json / tune history; nothing is typed in by hand.
"""
from __future__ import annotations

import argparse
import html
import json
from pathlib import Path

NAMES = {"lowpass": ("Low-pass", "C in the rungs"), "highpass": ("High-pass", "C in the end rings"),
         "bandpass": ("Band-pass", "C in rungs and end rings")}


def row(kind: str, case: Path, first: Path | None) -> dict:
    res = json.loads((case / "modes.json").read_text())
    des = res["design"]
    hp = res["homogeneous_pair"]
    st = res.get("solve_stats") or {}
    k1 = [r for r in res["comparison"] if r["k"] == 1][0]
    fg = None
    if first and first.exists():
        h = json.loads(first.read_text())[0]
        fg = {"c": h.get("c_leg") or h.get("c_ring"), "f": h["f_k1_mhz"], "src": first.parent.name}
    notes = des.get("notes", {})
    ladder_c = notes.get("c_leg_first_guess") or notes.get("c_ring_first_guess")
    return {"kind": kind, "case": case.as_posix(), "c_leg": des.get("c_leg"), "c_ring": des.get("c_ring"),
            "f": hp["f_mean_mhz"], "split": hp.get("split_mhz"), "cv": hp["roi_bt_cv"], "b1p": hp.get("b1plus_roi_cv"),
            "bm": hp.get("b1minus_over_b1plus_roi"), "ladder_at_c": k1["f_ladder_mhz"],
            "mesh": f"{'production' if des['h_coil'] <= 0.0051 and des['order'] == 2 else 'coarse'} "
                    f"({des['h_coil'] * 1e3:g}/{des['h_gap'] * 1e3:g}/{des['h_max'] * 1e3:g} mm, order {des['order']}, "
                    f"{st.get('unknowns', 0) / 1e6:.2f} M unknowns)",
            "fg": fg, "ladder_c": ladder_c, "order": [m for m in res["modes"]]}


def pf(c):
    return "-" if not c else f"{c * 1e12:.2f} pF"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--case", action="append", default=[], help="kind=case_dir[:first_guess_history.json]")
    ap.add_argument("--first-guess-c", action="append", default=[], help="kind=F (ladder first guess, when the case "
                    "was built with an explicit C)")
    ap.add_argument("--status", default="")
    ap.add_argument("--extra", default="", help="JSON file with extra caption lines per kind")
    a = ap.parse_args()
    out = Path(a.out)
    fgc = {k: float(v) for k, v in (x.split("=") for x in a.first_guess_c)}
    rows = []
    for spec in a.case:
        kind, rest = spec.split("=", 1)
        parts = rest.split(":")
        case = Path(parts[0] + (":" + parts[1] if len(parts) > 2 or (len(parts) == 2 and len(parts[0]) == 1) else ""))
        first = Path(parts[-1]) if len(parts) > 1 and parts[-1].endswith(".json") else None
        r = row(kind, case, first)
        if kind in fgc:
            r["ladder_c"] = fgc[kind]
        rows.append(r)
    extra = json.loads(Path(a.extra).read_text()) if a.extra else {}

    def tr(r):
        name, where = NAMES[r["kind"]]
        fg = r["fg"]
        fgtxt = (f"{pf(fg['c'])} &rarr; {fg['f']:.2f} MHz ({(fg['f'] / 127.732434 - 1) * 100:+.1f} %)"
                 f"<br><span class=m>{html.escape(fg['src'])}</span>") if fg else "-"
        cc = " / ".join(x for x in (f"rung {pf(r['c_leg'])}" if r["c_leg"] else "",
                                    f"ring {pf(r['c_ring'])}" if r["c_ring"] else "") if x)
        return (f"<tr><td><b>{name}</b><br><span class=m>{where}</span></td><td>{cc}</td>"
                f"<td><b>{r['f']:.2f} MHz</b><br><span class=m>split {r['split'] * 1e3:.0f} kHz</span></td>"
                f"<td>{pf(r['ladder_c'])}</td><td>{fgtxt}</td>"
                f"<td>{r['ladder_at_c']:.1f} MHz ({(r['f'] / r['ladder_at_c'] - 1) * 100:+.1f} %)</td>"
                f"<td>{r['cv'][0] * 100:.1f} / {r['cv'][1] * 100:.1f} %</td><td>{r['b1p'] * 100:.1f} %</td>"
                f"<td>{r['bm'] * 100:.1f} %</td><td class=m>{html.escape(r['mesh'])}</td></tr>"
                + (f"<tr><td></td><td colspan=9 class=note>{extra[r['kind']]}</td></tr>" if r["kind"] in extra else ""))

    cards = [
        ("fields_types.gif", "Hero: k = 1 pair in quadrature for low-pass, high-pass and band-pass (1080 x 1080). "
         "Top: currents; bottom: |B| and B direction in the axial plane. MP4: fields_types.mp4."),
        ("fields_lowpass.gif", "Low-pass alone (1080 x 1350, 4:5). MP4: fields_lowpass.mp4."),
        ("fields_highpass.gif", "High-pass alone (1080 x 1350, 4:5). MP4: fields_highpass.mp4."),
        ("fields_bandpass.gif", "Band-pass alone (1080 x 1350, 4:5). MP4: fields_bandpass.mp4."),
        ("coils_types.png", "(a) The three coils, capacitor gaps in yellow, with the tuned values."),
        ("spectra_types.png", "(b) Mode spectra, Palace vs ladder model at the same C, with k labels and the 127.7 MHz "
         "line. Low-pass: k = 1 lowest; high-pass and band-pass: k = 1 highest of the rung modes."),
        ("b1_types.png", "(c) Quadrature |B1+| in the central axial plane, same colour scale for all three."),
        ("porttune_schematic.png", "(d) circuitRF tuning network around the Palace 16-port coil (drawn with schemdraw "
         "from tune_unloaded_after.cnl; values parsed from the netlist)."),
        ("porttune_schematic_circuitrf_native.png", "(d, alt) The same netlist drawn by circuitRF itself "
         "(circuitrf netlist --to-schematic, then circuitrf render). Complete but hard to read at post size."),
        ("porttune_netlist_text.png", "(d, backup) The .cnl netlist as text."),
        ("porttune_sparams.png", "Existing candidate: port-tuning S-parameters (from birdcage-symmetry)."),
    ]
    cards = [c for c in cards if (out / c[0]).exists()]
    body = "\n".join(f'<figure><a href="{f}"><img src="{f}" loading="lazy"></a><figcaption><b>{f}</b> &middot; '
                     f'{html.escape(c)}</figcaption></figure>' for f, c in cards)
    status = f'<p class=status>{a.status}</p>' if a.status else ""
    page = f"""<!doctype html><html><head><meta charset="utf-8"><meta name=viewport content="width=device-width,initial-scale=1">
<title>Birdcage types</title><style>
:root{{--bg:#fff;--ink:#0b0b0b;--ink2:#52514e;--line:#e4e3df;--note:#fff8e1}}
@media (prefers-color-scheme:dark){{:root:not([data-theme=light]){{--bg:#151515;--ink:#eee;--ink2:#aaa;--line:#333;--note:#2b2716}}}}
:root[data-theme=dark]{{--bg:#151515;--ink:#eee;--ink2:#aaa;--line:#333;--note:#2b2716}}
body{{background:var(--bg);color:var(--ink);font:15px/1.45 system-ui,sans-serif;margin:0 auto;max-width:1300px;padding:20px 16px}}
h1{{font-size:24px;margin:0 0 4px}} .sub,.m{{color:var(--ink2)}} .m{{font-size:12px}}
.status{{background:var(--note);border:1px solid var(--line);padding:10px 12px;border-radius:6px}}
.tw{{overflow-x:auto}} table{{border-collapse:collapse;font-size:13.5px;margin:14px 0}} td,th{{border-bottom:1px solid var(--line);padding:6px 8px;text-align:left;vertical-align:top}}
th{{font-size:12px;color:var(--ink2);font-weight:600}} .note{{font-size:12.5px;color:var(--ink2)}}
.grid{{display:grid;grid-template-columns:repeat(auto-fill,minmax(380px,1fr));gap:18px}}
figure{{margin:0;border:1px solid var(--line);border-radius:6px;padding:8px;background:#fff}} img{{width:100%;display:block}}
figcaption{{font-size:13px;color:#333;margin-top:6px}}
</style></head><body>
<h1>Birdcage types: low-pass, high-pass, band-pass</h1>
<p class=sub>16-rung birdcage, 28 cm diameter, 25 cm long, 34 cm closed PEC shield, unloaded, target k = 1 at 127.7 MHz (1H at 3 T).
Palace FEM eigenmode solves, simulation only. Repo: birdcage-palace (tools/types_anim.py, types_figs.py, porttune_schematic.py, types_index.py).</p>
{status}
<div class=tw><table><tr><th>type</th><th>C (tuned)</th><th>k = 1 (Palace)</th><th>ladder first-guess C</th><th>Palace k = 1 at first-guess C</th>
<th>ladder k = 1 at tuned C (Palace vs ladder)</th><th>ROI std/mean |B_t| (each member)</th><th>quadrature B1+ ROI std/mean</th><th>|B1-|/|B1+|</th><th>mesh of the tuned result</th></tr>
{''.join(tr(r) for r in rows)}</table></div>
<p class=m>ROI: r &lt; 0.6 R, |z| &lt; L/6 (3D probe set). "Palace k = 1 at first-guess C" comes from the first solve of each tune
(mesh named under the value).</p>
<div class=grid>{body}</div></body></html>"""
    (out / "index.html").write_text(page, encoding="utf-8")
    print("wrote", out / "index.html")


if __name__ == "__main__":
    main()
