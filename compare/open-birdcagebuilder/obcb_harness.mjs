// Run Open BirdcageBuilder's capacitor calculation (github.com/dbnyu/OpenBirdcageBuilder, AGPL-3.0) from Node.
//
// Nothing from Open BirdcageBuilder is copied here: this script loads its unmodified source files from a local
// clone given on the command line, in the same order as its index.html, and makes the same calls as its
// get_gui_args() + calculate() (BBGuiHandler.js), with the GUI's unit conversions done here (cm -> m, MHz -> Hz,
// pF -> F).
//
//   node obcb_harness.mjs <path to OpenBirdcageBuilder clone> <cases.json>  > obcb_raw.json
//
// cases.json: [{"name", "config": lowpass|highpass|bandpass_leg|bandpass_er, "n_legs", "f_mhz", "r_coil_cm",
//               "r_shield_cm" (0 = no shield), "leg_len_cm", "leg_width_cm", "er_width_cm", "predcap_pf"?,
//               "bp_equal": true?}]
// bp_equal: band-pass with C_leg = C_ring (ratio 1). Open BirdcageBuilder takes one fixed capacitor and returns the
// other, so the harness solves predcap = returned C by bisection (each step is one unmodified calc_capacitor() call).
import fs from "node:fs";
import path from "node:path";
import vm from "node:vm";

const [root, casesFile] = process.argv.slice(2);
const files = ["debug.js", "helpers.js", "BCLegs.js", "BCEndRings.js", "BirdcageBuilder.js"];  // index.html order
const ctx = vm.createContext({ Math, Number, Array, console: { log() {}, warn() {}, assert() {} } });
for (const f of files) vm.runInContext(fs.readFileSync(path.join(root, f), "utf8"), ctx, { filename: f });

// Same sequence as BBGuiHandler.get_gui_args() for rectangular legs and end rings, then calc_capacitor().
vm.runInContext(`
function obcb_run(c, predcap_pf) {
  const bb = new BirdcageBuilder(c.n_legs, mhz2hz(c.f_mhz), cm2meters(c.r_coil_cm), cm2meters(c.r_shield_cm), c.config);
  bb.legs.set_legs_rect(cm2meters(c.leg_len_cm), cm2meters(c.leg_width_cm));
  bb.endrings.set_endrings_rect(cm2meters(c.er_width_cm));
  bb.endrings.calc_er_mutual_inductance(bb.legs);
  if (c.config === "bandpass_leg" || c.config === "bandpass_er") bb.predcap = predcap_pf * 1e-12;
  const cap = bb.calc_capacitor();
  return { cap_pf: 1e12 * cap,
           leg_self_nH: 1e9 * bb.legs.leg_self_inductance, er_self_nH: 1e9 * bb.endrings.er_self_inductance,
           leg_mutual_nH: 1e9 * bb.legs.leg_mutual_inductance,
           er_mutual_nH: 1e9 * bb.endrings.er_mutual_inductance[bb.n_legs / 4] };  // the 4 values the GUI shows
}`, ctx);

// Self-check against the repo's own unit test (qunit_test_all_math.html, "App Defaults Tests"): 205.0515 pF.
const chk = ctx.obcb_run({ n_legs: 8, f_mhz: 21.28, r_coil_cm: 10, r_shield_cm: 0, config: "lowpass",
                           leg_len_cm: 10, leg_width_cm: 2, er_width_cm: 1 });
const out = { obcb_selfcheck: { case: "repo unit test 'App Defaults Tests'", expected_pf: 205.0515, got_pf: chk.cap_pf },
              cases: [] };

for (const c of JSON.parse(fs.readFileSync(casesFile, "utf8"))) {
  if (c.bp_equal) {
    // g(C) = calc_capacitor(predcap = C) - C. The returned C has a pole (it changes sign) at the fixed C where
    // the other capacitor alone would tune the coil, so scan a log grid for a sign change of g with a positive
    // returned C on both sides of it, then bisect.
    const g = (x) => ctx.obcb_run(c, x).cap_pf - x;
    const grid = Array.from({ length: 3001 }, (_, i) => Math.pow(10, -1 + 4 * i / 3000));  // 0.1 .. 1000 pF
    const roots = [];
    for (let i = 0; i + 1 < grid.length; i++) {
      const a = grid[i], b = grid[i + 1];
      if (ctx.obcb_run(c, a).cap_pf > 0 && ctx.obcb_run(c, b).cap_pf > 0 && g(a) * g(b) <= 0) roots.push([a, b]);
    }
    if (roots.length !== 1) { out.cases.push({ ...c, error: `found ${roots.length} positive fixed points` }); continue; }
    let [lo, hi] = roots[0];
    for (let i = 0; i < 200; i++) { const m = 0.5 * (lo + hi); (g(lo) * g(m) <= 0) ? (hi = m) : (lo = m); }
    const x = 0.5 * (lo + hi);
    out.cases.push({ ...c, predcap_pf: x, result: ctx.obcb_run(c, x), fixed_point_residual_pf: g(x) });
  } else {
    out.cases.push({ ...c, result: ctx.obcb_run(c, c.predcap_pf ?? -1) });
  }
}
console.log(JSON.stringify(out, null, 1));
