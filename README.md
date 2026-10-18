# birdcage-palace

Project page: https://benjhardy.github.io/projects/birdcage-palace/ · part of [benjhardy.github.io](https://benjhardy.github.io)

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
* **k = 0** has the same current in every rung. In the low-pass coil that mode closes through the shield; in the
  high-pass coil the k = 0 modes are end-ring modes with no rung current (counter- and co-rotating ring currents).
* Where the capacitors sit sets the ordering. **Low-pass** (C in the rungs): k = 1 is the lowest resonance above
  DC and k = N/2 the highest. **High-pass** (C in the end-ring segments): the order flips, the end-ring k = 0 modes
  are on top and k = 1 is the highest rung mode below them. **Band-pass** has both.

So a full-wave solve returns a crowded spectrum (20 modes between 128 and 483 MHz for the tuned default coil), and
only one pair images. This tool finds that pair automatically.

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

`tune` starts from the plain (f/f0)² scaling, or from `--slope` (d(1/f²)/dC measured on a cheaper mesh), then takes
secant steps on 1/f². `tools/make_video.py`, `tools/export_docs.py` and `tools/hero.py` rebuild the clip, the
figures in `docs/img` and the hero image from a solved case.

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

At the first-guess capacitor Palace puts k = 1 13 % below the ladder value. Fitting 1/f² = L_eff (C + C_stray) to the
three fine-mesh solves gives L_eff = 511 nH, against 553 nH from the ladder model (most likely because of the shield
end caps, which the model leaves out and which lower the ring inductance), plus C_stray = 1.21 pF in parallel with each lumped capacitor, which is
two-thirds of its final value. At 3 T a 16-rung low-pass head coil needs only a couple of picofarads per rung, so the
distributed capacitance of the strips (to the shield and to each other) matters as much as the lumped part. This is the usual argument
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
tuning steps on the production mesh took 14.2 and 13.6 min each (4 MPI ranks on an 8-thread machine that was at
times shared with another Palace job). Palace estimates the peak memory at 16.2 GB. Retuning took two solves: the first step used the 1/f² slope
from a pair of coarse solves and landed at +0.24 %, and a secant step then landed at +0.001 %.

![Mesh: y = 0 slice through the volume mesh and the surface mesh at the capacitor gaps](docs/img/mesh.png)

## Loaded check: saline phantom (coarse mesh)

A cylindrical phantom (radius 8 cm, length 20 cm, εr 78, σ 0.6 S/m) inside the tuned coil (C = 1.83 pF), on the
coarse order-2 mesh, compared with the same mesh unloaded:

| | unloaded | loaded |
|---|---|---|
| k = 1 pair | 128.09 MHz | 127.33 MHz (−0.6 %) |
| Q (phantom loss only; conductors are lossless) | ∞ | 66 |
| quadrature B1+, ROI std/mean | 2.4 % | 18 % |
| \|B1−\| / \|B1+\| in the ROI | 4.0 % | 28 % |
| solve | 484 k unknowns, 3.8 min on 4 ranks, 6.6 GB | 498 k unknowns, 5.4 min on 6 ranks, 7.9 GB |

The phantom makes the eigenproblem complex (Palace reports Q directly). At 128 MHz the wavelength in the phantom is
about 27 cm, comparable to its 16 cm diameter, so the quadrature field brightens at the centre and B1+ and B1−
separate. That is the familiar 3 T head-imaging inhomogeneity, which the empty-coil maps above do not show. The ROI
(r < 8.4 cm) extends 4 mm past the phantom radius. The loaded solve requested only the k = 1 pair (2 modes,
tolerance 1e-6); with 12 modes the quadratic eigenproblem converged too slowly to be practical here.

![Quadrature B1+ unloaded vs phantom-loaded](docs/img/loaded_b1plus.png)

## High-pass check (coarse mesh, not tuned)

`birdcage build --type highpass` moves the capacitors into the 32 end-ring segments (ladder first guess 36.9 pF each).
One coarse solve (order 1, 100 k unknowns, 2.2 min on 3 ranks) shows the inverted spectrum the theory predicts:
k = 8 is now the lowest mode (58.5 MHz) and k = 1 the highest rung mode (135.2 MHz, ladder 127.7 MHz, ROI std/mean
3.6 %). Above it sit the two end-ring modes: co-rotating at 163.4 MHz (ladder 161.8 MHz; ring currents in the same
sense, axial field, transverse fraction 0.10) and counter-rotating at 163.7 MHz (ladder k = 0 at 155.7 MHz). The rung
modes and the counter-rotating mode are 5 to 13 % above the ladder values (the co-rotating mode 1 %). The shift has
the opposite sign to the low-pass coil: with 37 pF per gap the stray capacitance hardly matters, and the shield end caps, which the ladder model leaves out, probably lower the
inductance. This case is a functional
check on a coarse mesh only; it has not been retuned or refined.

![High-pass mode spectrum](docs/img/highpass_spectrum.png)


## Symmetry planes

Simulation only. The coil has the mirror symmetries of a square (16 rungs), so one linear member of the k = 1 pair can
be solved on a quarter of the domain, with a PEC wall on one cut plane and a PMC wall on the other. Same mesh density
as the production full model (5 / 2.5 / 35 mm, order 2), same capacitor (1.828 pF). Numbers are in
[`results/lp16_symmetry.json`](results/lp16_symmetry.json).

![Quarter model: x = 0 plane PEC (orange), y = 0 plane PMC (green)](docs/img/symmetry_quarter_model.png)

```bash
birdcage build ... --symmetry quarter --sym-bc x=pec,y=pmc --out runs/q_em   # B1 along y
birdcage variant --case runs/q_em --sym-bc x=pmc,y=pec --out runs/q_me        # same mesh, the partner (B1 along x)
birdcage solve --case runs/q_em && birdcage modes --case runs/q_em             # fields and rung currents unfolded to the full coil
```

**Which plane is which.** Take the member with rung currents I_n ∝ cos φ_n. A +z current at φ = 0 and a −z current at
φ = π both give B along −y at the centre, so this member's B1 is along y. Mirror it in the plane x = 0 (φ → π − φ):
cos changes sign, so the rung currents, which run parallel to that plane, are antisymmetric across it. The image of a
tangential current in a perfect electric conductor is reversed, so x = 0 is a **PEC plane** (tangential E = 0; B1
lies in it). Mirror it in y = 0 (φ → −φ): cos is unchanged, the currents are symmetric, and y = 0 is a **PMC plane**
(tangential H = 0; B1 crosses it at right angles). In short, the PEC plane contains B1 and the PMC plane is
perpendicular to it. The x-polarised member (sin φ_n) has the two swapped.

The full-model field confirms this. Rotating the solved pair so that one member has B1 along y and mirroring its ROI
field (938 point pairs per plane):

| member | plane | residual if PEC, \|B(Mr) − M B(r)\| / \|B\| | residual if PMC, \|B(Mr) + M B(r)\| / \|B\| | plane is |
|---|---|---|---|---|
| B1 along y | x = 0 | 0.66 % | 200 % | PEC |
| B1 along y | y = 0 | 200 % | 0.67 % | PMC |
| B1 along x | x = 0 | 200 % | 0.66 % | PMC |
| B1 along x | y = 0 | 0.65 % | 200 % | PEC |

Its rung currents fit −cos φ_n with a residual of 0.03 % (B along +y needs I ∝ −cos φ_n with this sign convention). In
Palace the PEC plane goes into `Boundaries.PEC`; the PMC plane is simply left without a boundary condition, because
PMC is the natural boundary of the curl-curl formulation (Palace says so in its log: "PMC/ZeroCharge condition is
assumed"). Fields and capacitor currents are unfolded to the full coil with the mirror rules B(Mr) = +M B(r) across a
PEC plane and −M B(r) across a PMC plane (rung current odd and even, respectively).

**Where the planes go.** With rungs at φ_n = 2πn/N the planes x = 0 and y = 0 pass through four rung centres, and
the plane cuts those rungs and their capacitors lengthwise. Palace defines a lumped element by its circuit value and
converts it to a per-square value with the element's own width and length, Cs = C · l / w. Half the width at the same
voltage carries half the current, so the half element needs **C / 2** (Cs is then unchanged; the Palace log shows
0.7616 pF/sq on every element). On a PEC plane the half rung is shorted by the wall and carries no current, which is
what the cos φ member has there anyway. Both versions were solved (quarter, x = 0 PEC, y = 0 PMC):

| quarter model | k = 1 (MHz) |
|---|---|
| planes between rungs (cage rotated by π/16) | 127.61 |
| planes through rungs, half element C / 2 | 127.70 |
| planes through rungs, full C on the half element (wrong) | 117.93 (−7.6 %) |
| full model, the pair | 127.73 (127.726 / 127.741) |

The default is **between rungs** (rotate the cage by π/N): no element is cut, each quarter holds four whole rungs with
the same C, and the z = 0 rule below is the only special case. It needs N divisible by 4.

**The z = 0 midplane (eighth model).** The capacitor gaps sit on it. In every ladder mode the rung current crosses
z = 0 and continues in the same direction, and the two end rings carry opposite currents, so the tangential E is odd
in z: z = 0 is a PEC plane. A PEC wall forces only the tangential E to zero; the gap field E_z is normal to the wall
and is allowed, so the capacitor is not shorted. The wall cuts the element across its length: each half sees half the
voltage at the full current, so the half element needs **2C** (again the same Cs). This is implemented
(`--symmetry eighth`). The z-even family, which a PMC wall at z = 0 would keep, has no voltage across the gaps; it
contains the 222 MHz pair whose capacitors carry no current, which the eighth model indeed no longer returns.

**Symmetry filters the spectrum.** A quarter model only supports modes with its own parities. For the rotated cage
(rungs at (2n+1)π/16) the four combinations split the 16 ladder modes into four families of four:

| x = 0 | y = 0 | rung currents | ladder modes kept | k = 1 | modes found below 490 MHz (MHz) |
|---|---|---|---|---|---|
| PEC | PMC | cos kφ, k odd | k = 1, 3, 5, 7 | B1 along y | 127.61, 222.46 (no cap. current), 242.07, 279.94, 294.84 |
| PMC | PEC | sin kφ, k odd | k = 1, 3, 5, 7 | B1 along x | 127.66, 222.51 (no cap. current), 242.07, 279.94, 294.70 |
| PMC | PMC | cos kφ, k even | k = 0, 2, 4, 6 | none | 202.72, 265.18, 289.46, 393.19 (no cap. current), 483.51 |
| PEC | PEC | sin kφ, k even | k = 2, 4, 6, 8 | none | 202.70, 265.13, 289.36, 296.51, 393.06 (no cap. current) |

(k = N/2 = 8 is a sine pattern on the rotated cage, cos 8φ_n = 0, so it appears only in PEC/PEC; k = 0 is a cosine
pattern.) The union of the four runs is exactly the full model's 20 modes: every one is matched by type and k, the
largest frequency difference is 0.11 %, and no run returns an extra mode below 490 MHz. A degenerate pair is split
across two runs, so the pair splitting itself (15 kHz in the full model, a mesh effect) needs both runs; on one
quarter mesh the two k = 1 members differ by 56 kHz (0.04 %), because that mesh is not symmetric about the diagonal.
For finding the imaging mode this filtering helps: in PEC/PMC the homogeneous mode is the lowest eigenvalue and only
four other modes lie below 300 MHz.

![Which modes each PEC/PMC combination keeps](docs/img/symmetry_spectrum.png)

**Full vs quarter vs eighth** (same mesh density, order 2, C = 1.828 pF; MUMPS, 4 ranks):

| | full | quarter (x = 0 PEC, y = 0 PMC) | eighth (+ z = 0 PEC) |
|---|---|---|---|
| k = 1 (MHz) | 127.73 (pair mean) | 127.61 (−0.10 %) | 127.68 (−0.05 %) |
| ROI std/mean of \|B_t\|, linear member | 4.03 % | 4.04 % | 4.08 % |
| quadrature B1+, ROI std/mean | 2.35 % | 2.32 % (two quarter runs) | |
| tetrahedra | 169 k | 44 k | 23 k |
| unknowns | 1.12 M | 0.30 M (3.8×) | 0.16 M (7.2×) |
| peak memory (Palace estimate) | 16.4 GB | 3.7 GB (4.4×) | 2.1 GB (7.8×) |
| Palace time, whole run | 12.8 min (20 modes) | 1.2 min (10 modes) | 4.0 min (10 modes), CPU shared |
| factorisation (setup + preconditioner) | 359 s | 25 s (14×) | 77 s, CPU shared |

The frequency differences are within the 0.4 % discretisation error of this mesh density (see the mesh-convergence
table above) and are of the size one expects between two different meshes at that density. The quarter run asked for
10 modes and the full run for 20; the factorisation row compares the part that does not depend on the mode count.
On a free machine one quarter run took 71 s, so all four (the whole spectrum) take under 5 min against 12.8 min for one
full run. All runs used 4 MPI ranks. Another job (8 ranks) shared the 8-thread machine during the eighth-model run
and during the other three quarter runs, so their times are not comparable and are not used here; the unknowns and
memory are unaffected.

![k = 1 field, full model vs mirrored quarter model](docs/img/symmetry_b1_compare.png)

![Cost of one eigen solve: full, quarter, eighth](docs/img/symmetry_costs.png)

## Limits

* Simulation only; nothing here is measured. Conductors are lossless PEC sheets and the capacitors are ideal, so
  the unloaded coil has no Q, and the loaded Q counts phantom loss only. The phantom case uses the coarse mesh. Copper and capacitor ESR, cables, baluns, matching and the drive ports are not modelled.
* The shield is a closed PEC can (end caps 10 cm beyond the rings), not the open-ended shield of a real scanner. It
  probably lowers the ring inductance a little. `--shield-ends absorbing` and `--shield 0` (absorbing outer boundary)
  are implemented but untested, and they turn the eigenproblem complex.
* The quadrature B1+ map is a superposition of the two k = 1 eigenmodes, which is what an ideal, decoupled
  quadrature drive excites at resonance. It is not a driven two-port solve.
* The ladder model is quasi-static and ignores stray capacitance. It is a starting point, not the answer.
* Only the low-pass coil has been taken through mesh refinement and the tune loop. The high-pass coil has one coarse
  solve (above). Band-pass geometry and ladder formulas are implemented but have not been solved in Palace yet.

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
