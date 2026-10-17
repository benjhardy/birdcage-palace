# birdcage-palace

**Live / project page: https://benjhardy.github.io/birdcage-palace/** *(placeholder, not yet published)*

An open-source, scriptable builder for MRI birdcage RF coils. It sizes the capacitors with a ladder-network model,
builds the coil in [Gmsh](https://gmsh.info), solves the full-wave eigenmodes with the
[Palace](https://github.com/awslabs/palace) finite-element solver, and then picks out the one mode that
matters for imaging, the homogeneous k = 1 pair, from all the others.

![16-rung low-pass birdcage in its RF shield](docs/img/coil.png)

The default design is a 16-rung low-pass head coil (rung diameter 28 cm, 25 cm between end-ring centrelines,
RF shield 34 cm) tuned to the proton Larmor frequency at 3 T: 42.577 MHz/T x 3 T = 127.7 MHz.

> Not to be confused with **BirdcageBuilder** (Chin, Collins, Li, Dardzinski and Smith, Penn State, 2002), the
> well-known ladder-network design tool this project's analytic step follows. birdcage-palace adds a 3D
> full-wave eigenmode solve and automatic mode identification on top of that kind of model.

## Why the mode matters

A birdcage with N rungs is a periodic ladder network, so its resonances are discrete Fourier modes. In mode k the
rung currents go as I_n ∝ cos(2πkn/N + φ), k = 0 ... N/2 (Hayes et al. 1985; Tropp 1989; Leifer 1997):

* **k = 1** puts one full period of current around the coil, the 2D analogue of a cos θ surface current. It
  makes a uniform transverse B1 in the bore. It is a degenerate pair (cos and sin patterns, 90 degrees apart),
  so two drive ports 90 degrees apart in space and phase give a circularly polarised B1+ (quadrature drive).
* **k ≥ 2** modes have k current periods. Their field is a multipole that vanishes at the centre and is
  confined near the conductors. Every k between 1 and N/2 − 1 is a degenerate pair; k = N/2 is single.
* **k = 0** has no net rung current pattern (low-pass: rungs all in phase, closing through the shield; high-pass:
  counter- and co-rotating end-ring currents).
* Where the capacitors sit sets the ordering. **Low-pass** (C in the rungs): k = 1 is the lowest resonance above
  DC and k = N/2 the highest. **High-pass** (C in the end-ring segments): the order flips, the end-ring k = 0 modes
  are on top and k = 1 is the highest rung mode below them. **Band-pass** has both.

So a full-wave solve returns a crowded spectrum (20 modes between 111 and 466 MHz for the default coil), and only one
pair images. This tool finds that pair automatically.

## What it does

1. **Analytic first guess** (`birdcage_palace/ladder.py`). Neumann double integrals over the rung and ring
   centre lines, regularised with the geometric-mean distance of a thin strip, give the self and mutual
   inductances. The shield is represented by image conductors. Their discrete Fourier transforms Lleg_k and
   Lring_k give every mode frequency in closed form:
   low-pass ω_k² = 1 / (C (Lleg_k + Lring_k / (2 sin²(πk/N)))), high-pass ω_k² = 1 / (C (Lring_k + 2 sin²(πk/N) Lleg_k)),
   band-pass a combination of the two. The capacitor that puts k = 1 on the target frequency follows directly.
2. **Geometry and mesh** (`geometry.py`). The rungs and end rings are zero-thickness PEC strips on the coil
   cylinder. Each capacitor is a small patch across a 5 mm gap (in the rungs for low-pass, in the ring segments for
   high-pass, both for band-pass). The domain is the inside of the RF shield, closed by PEC end caps 10 cm beyond the
   rings. The mesh is graded to 5 mm on the strip edges and 2.5 mm at the gaps, and grows to 35 mm in the bore. An
   optional saline phantom (εr 78, σ 0.6 S/m) can be added.
3. **Palace eigenmode solve** (`palace.py`). Each capacitor is a Palace `LumpedPort` with only `C` set, so the
   capacitors enter the generalised eigenproblem. The solver is shift-and-invert Krylov-Schur (SLEPc) above a target
   of 0.55 f0, with 20 modes, order-2 Nédélec elements, and a direct MUMPS factorisation (`MGMaxLevels = 1`).
   The closed shield keeps the problem real and lossless, so unloaded Q is infinite by construction; a phantom or
   an absorbing boundary (`--shield 0` or `--shield-ends absorbing`) makes it complex and Palace reports Q.
4. **Mode identification** (`modes.py`). For every eigenmode:
   * the 16 capacitor currents (`port-I.csv`) give the rung-current pattern; its spatial DFT gives k and a purity;
   * about 4,000 field probes on the axial plane, the coronal plane and a 3D ROI (r < 0.6 R, |z| < L/6) give the
     ROI std/mean of |B_transverse| and the transverse fraction of |B|²;
   * modes whose capacitors carry almost no current are flagged as non-ladder resonances (end-ring / rung-half
     modes set by stray capacitance).

   The homogeneous pair is the k = 1 mode with a transverse field and the lowest ROI std/mean, plus its partner.
   The pair is then combined 90 degrees apart to give the quadrature B1+ and B1− maps.
5. **Retuning** (`birdcage tune`). The mesh is reused, only C changes. Because 1/f² = L_eff (C + C_stray) is
   linear in C, a secant step on 1/f² lands on target in one or two solves.

## Quick start

```bash
pip install -e .            # numpy, matplotlib, gmsh (+ pyvista for the 3D renders)
birdcage build --legs 16 --diameter 0.28 --length 0.25 --type lowpass --b0 3.0 --shield 0.34 --out runs/lp16
birdcage solve --case runs/lp16 --np 4      # Palace (on Windows it runs inside WSL; see below)
birdcage modes --case runs/lp16             # modes.json + figures in runs/lp16/figs
birdcage tune  --case runs/lp16 --tol 0.005 # rescale C until the k=1 pair is within 0.5 % of 127.7 MHz
```

`birdcage build --help` lists every option: rung and ring widths, gap, high-pass and band-pass, shield size and end
caps, phantom, mesh sizes and FE order. Palace is called as `run_palace config.json <np>` (a wrapper around the
Palace Apptainer image `palace -np <np> config.json`). On Windows the command runs inside WSL2.

## Results: default 16-rung low-pass head coil at 3 T

Simulation only. Coil: 16 rungs, 12 mm wide, on a 28 cm diameter, 25 cm between end-ring centrelines, 15 mm end
rings, one 5 mm capacitor gap at the middle of each rung, closed PEC shield 34 cm in diameter and 45 cm long, no load.
Numbers are in [`results/lp16_3T.json`](results/lp16_3T.json).

**The homogeneous mode**

| | |
|---|---|
| k = 1 pair | **127.73 MHz** (127.726 and 127.741 MHz; split 15 kHz, 0.012 %) |
| capacitor per rung | **1.83 pF** (ladder first guess 2.81 pF) |
| ROI std/mean of \|B_t\|, one linear mode | 4.0 % (ROI: r < 0.6 R, \|z\| < L/6) |
| quadrature B1+ (pair driven 90 degrees apart) | ROI std/mean 2.3 %, \|B1−\| / \|B1+\| = 4.0 % |
| transverse fraction of \|B\|² in the ROI | 0.997 |

![Mode spectrum: ladder model vs Palace, before and after retuning](docs/img/spectrum_compare.png)

**Ladder model vs Palace** (same mesh, order 2, 1.12 M unknowns)

| k | ladder, C = 2.81 pF | Palace, C = 2.81 pF | ladder, C = 1.83 pF | Palace, C = 1.83 pF |
|---|---|---|---|---|
| 1 | **127.7** | 111.1 | 158.3 | **127.7** |
| 2 | 189.9 | 173.7 | 235.3 | 202.8 |
| 3 | 220.9 | 206.2 | 273.8 | 242.2 |
| 4 | 239.5 | 225.4 | 296.7 | 265.2 |
| 5 | 251.6 | 237.8 | 311.8 | 280.0 |
| 6 | 259.5 | 245.8 | 321.6 | 289.5 |
| 7 | 264.1 | 250.4 | 327.3 | 294.9 |
| 8 | 265.6 | 251.9 | 329.1 | 296.6 |

Frequencies in MHz; Palace pairs are averaged (the members agree to within 0.04 MHz). Palace also finds modes the
ladder model does not have: k = 0 with all rungs in phase, closing through the shield (465.7 MHz at 2.81 pF, 483.3 MHz
at 1.83 pF), and two pairs at 222.6 and 393.3 MHz whose capacitors carry about 0.1 % of the current of the ladder
modes. Those two pairs do not move when C changes, so they are resonances of the end rings and rung halves with
their stray capacitance, not birdcage modes.

The ladder model puts k = 1 13 % too high at the first-guess capacitor. Fitting 1/f² = L_eff (C + C_stray) to the
three fine-mesh solves gives L_eff = 511 nH, against 553 nH from the ladder model (the shield end caps, which the
model leaves out, lower the ring inductance), plus C_stray = 1.21 pF in parallel with each lumped capacitor, which is
two-thirds of its final value. At 3 T a 16-rung low-pass head coil needs only a couple of picofarads per rung, so the
stray capacitance to the shield and between the strips matters as much as the lumped part. This is the usual argument
for high-pass or band-pass designs at 3 T and above.

![|B| maps per mode](docs/img/b1_maps.png)

The maps show why the mode choice matters. Only k = 1 fills the ROI with a uniform transverse field (std/mean 4 %).
k = 2 has a null at the centre (35 %), k = 8 keeps its field at the conductors (126 %), and the in-phase k = 0 mode
gives an azimuthal field that vanishes on the axis (33 %).

![Rung-current patterns](docs/img/currents.png)

**Mesh convergence of the k = 1 frequency** (C = 2.81 pF)

| mesh (strip edge / gap / bore, mm) | tets | order | unknowns | k = 1 (MHz) | change | ranks | wall | peak RAM |
|---|---|---|---|---|---|---|---|---|
| 8 / 4 / 40 | 72 k | 1 | 92 k | 112.95 | | 4 | 3.0 min | 1.4 GB |
| 8 / 4 / 40 | 72 k | 2 | 484 k | 111.54 | −1.25 % | 4 | 6.6 min | 6.6 GB |
| 5 / 2.5 / 35 | 169 k | 2 | 1.12 M | 111.08 | −0.41 % | 4 | 20.4 min | 16.2 GB |

The last step moves k = 1 by 0.4 % and k = 8 by 0.7 %. Expect the production numbers to carry about that much
discretisation error; the next refinement step would likely need more RAM than the 32 GB machine used here has. The
tuning steps on the production mesh took 14.2 and 13.6 min each (4 MPI ranks, a second Palace job sharing the 8-thread
machine). Palace estimates the peak memory at 16.2 GB. Retuning took two solves: the first step used the 1/f² slope
from a pair of coarse solves and landed at +0.24 %, and a secant step then landed at +0.001 %.

![Mesh: y = 0 slice through the volume mesh and the surface mesh at the capacitor gaps](docs/img/mesh.png)

LOADED_PLACEHOLDER

## Limits

* Simulation only; nothing here is measured. Conductors are lossless PEC sheets and the capacitors are ideal, so
  the unloaded coil has no Q. Copper and capacitor ESR, cables, baluns, matching and the drive ports are not modelled.
* The shield is a closed PEC can (end caps 10 cm beyond the rings), not the open-ended shield of a real scanner. That
  lowers the ring inductance a little; `--shield-ends absorbing` is available but turns the eigenproblem complex.
* The quadrature B1+ map is a superposition of the two k = 1 eigenmodes, which is what an ideal, decoupled
  quadrature drive excites at resonance. It is not a driven two-port solve.
* The ladder model is quasi-static and ignores stray capacitance. It is a starting point, not the answer.
* Only the low-pass coil has been taken through the full tune loop here. High-pass and band-pass geometry, capacitor
  placement and ladder formulas are implemented, but their Palace results are not yet verified.

## References

1. C. E. Hayes, W. A. Edelstein, J. F. Schenck, O. M. Mueller, M. Eash, "An efficient, highly homogeneous
   radiofrequency coil for whole-body NMR imaging at 1.5 T," *J. Magn. Reson.* 63, 622–628 (1985).
2. J. Tropp, "The theory of the bird-cage resonator," *J. Magn. Reson.* 82, 51–62 (1989).
3. M. C. Leifer, "Resonant modes of the birdcage coil," *J. Magn. Reson.* 124, 51–60 (1997).
4. C.-L. Chin, C. M. Collins, S. Li, B. J. Dardzinski, M. B. Smith, "BirdcageBuilder: design of specified-geometry
   birdcage coils with desired current pattern and resonant frequency," *Concepts Magn. Reson.* 15(2), 156–163 (2002),
   doi:10.1002/cmr.10030.
5. J. Jin, *Electromagnetic Analysis and Design in Magnetic Resonance Imaging*, CRC Press (1998).
6. F. W. Grover, *Inductance Calculations*, Van Nostrand (1946) (strip self-inductance, geometric mean distance).
7. Palace: 3D finite element solver for computational electromagnetics, AWS Labs, https://github.com/awslabs/palace (version 0.16 used here).
8. C. Geuzaine, J.-F. Remacle, "Gmsh: a 3-D finite element mesh generator with built-in pre- and post-processing
   facilities," *Int. J. Numer. Meth. Eng.* 79(11), 1309–1331 (2009).

## License

MIT. Palace (Apache-2.0) and Gmsh (GPL) are used as external tools and are not redistributed here.
