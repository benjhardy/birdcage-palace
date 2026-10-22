"""Coil design parameters (SI units throughout) and their JSON round trip."""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

GAMMA_BAR_1H = 42.577478e6  # proton gyromagnetic ratio / 2 pi, Hz/T (CODATA 2018)


@dataclass
class Phantom:
    radius: float = 0.08          # m
    length: float = 0.20          # m
    eps_r: float = 78.0           # saline-like at 128 MHz
    sigma: float = 0.6            # S/m
    # How Palace sees the loss: "conductivity" (sigma as a conductivity: damping term, quadratic eigenproblem) or
    # "losstan" (the same loss as a loss tangent tan d = sigma / (2 pi f_ref eps0 eps_r), exact at f_ref: the
    # eigenproblem stays linear, so the shift-and-invert solve is a plain direct solve).
    loss_model: str = "conductivity"
    f_ref: float | None = None    # Hz, frequency of the loss-tangent conversion; None = the Larmor frequency


@dataclass
class CoilDesign:
    legs: int = 16                # number of rungs (even)
    diameter: float = 0.28        # m, rung centre diameter
    length: float = 0.25          # m, rung length between end-ring centrelines
    leg_width: float = 0.012      # m, copper strip width of the rungs
    ring_width: float = 0.015     # m, copper strip width of the end rings
    gap: float = 0.005            # m, capacitor gap length (along the current)
    coil_type: str = "lowpass"    # lowpass | highpass | bandpass
    b0: float = 3.0               # T
    shield_diameter: float | None = 0.34   # m; None = unshielded (absorbing outer boundary)
    shield_length: float | None = None     # m; default length + 0.20
    shield_ends: str = "pec"      # pec (closed can) | absorbing
    # unshielded only: air cylinder closed by a second-order absorbing boundary. Default (None): a gap of
    # max(0.15 m, lambda0 / 4) around the coil, radially and past each end-ring centreline (geometry.air_clearance)
    air_radius: float | None = None
    air_length: float | None = None
    c_leg: float | None = None    # F, capacitor in each rung (lowpass/bandpass); None = analytic first guess
    c_ring: float | None = None   # F, capacitor in each end-ring segment (highpass/bandpass)
    bandpass_ratio: float = 1.0   # bandpass only: c_ring / c_leg for the first guess
    phantom: Phantom | None = None
    # mesh
    h_coil: float = 0.005         # m, target edge length on the conductors
    h_gap: float = 0.0025         # m, edge length at the capacitor gaps
    h_max: float = 0.035          # m, far-field edge length
    order: int = 2                # FE order in Palace
    # symmetry (low-pass only): model 1/2, 1/4 or 1/8 of the coil, cut by the planes x=0 (half), x=0 and y=0
    # (quarter), plus z=0 (eighth). Each plane is PEC or PMC (sym_bc); PMC is Palace's natural boundary.
    symmetry: str = "none"        # none | half | quarter | eighth | zhalf (z >= 0 only: keeps all rungs)
    sym_cut: str = "between"      # between (planes between rungs: cage rotated by pi/N) | through (planes through rungs)
    sym_bc: dict = field(default_factory=lambda: {"x": "pec", "y": "pmc", "z": "pec"})
    rung_offset: float | None = None  # rad, angle of rung 0; None = 0, or pi/N when sym_cut == "between" and symmetric
    c_leg_rungs: list | None = None   # F per rung (full-coil rung index), overrides c_leg; for asymmetric tunings
    # coil losses (None = lossless): copper surface impedance on the strips and the shield (Palace Conductivity
    # boundary, frequency dependent: nonlinear eigensolve), and a series ESR per capacitor (entered in Palace as the
    # equivalent parallel R at the Larmor frequency, R_p = ESR (1 + Q_c^2), Q_c = 1 / (w0 C ESR))
    copper_sigma: float | None = None  # S/m, e.g. 5.8e7
    esr: float | None = None           # ohm, e.g. 0.68
    notes: dict = field(default_factory=dict)

    @property
    def radius(self) -> float:
        return self.diameter / 2

    @property
    def f_larmor(self) -> float:
        return GAMMA_BAR_1H * self.b0

    @property
    def planes(self) -> list[str]:
        return {"none": [], "half": ["x"], "quarter": ["x", "y"], "eighth": ["x", "y", "z"], "zhalf": ["z"]}[self.symmetry]

    @property
    def phi0(self) -> float:
        if self.rung_offset is not None:
            return self.rung_offset
        return 3.141592653589793 / self.legs if (self.symmetry != "none" and self.sym_cut == "between") else 0.0

    @property
    def shield_len(self) -> float | None:
        if self.shield_diameter is None:
            return None
        return self.shield_length if self.shield_length else self.length + 0.20

    def validate(self) -> None:
        if self.legs < 4 or self.legs % 2:
            raise ValueError("legs must be an even number >= 4")
        if self.coil_type not in ("lowpass", "highpass", "bandpass"):
            raise ValueError("coil_type must be lowpass, highpass or bandpass")
        if self.shield_diameter is not None and self.shield_diameter <= self.diameter + 2 * 0.005:
            raise ValueError("shield must clear the coil")
        if self.leg_width * self.legs >= 0.8 * 3.14159 * self.diameter:
            raise ValueError("rungs overlap")
        if self.phantom is not None and self.phantom.loss_model not in ("conductivity", "losstan"):
            raise ValueError("phantom loss_model must be conductivity or losstan")
        if self.symmetry not in ("none", "half", "quarter", "eighth", "zhalf"):
            raise ValueError("symmetry must be none, half, quarter, eighth or zhalf")
        if self.symmetry != "none":
            if self.coil_type != "lowpass":
                raise ValueError("symmetry planes are implemented for the low-pass coil only")
            if self.legs % 4:
                raise ValueError("symmetry planes need a multiple of 4 rungs")
            if self.sym_cut not in ("between", "through"):
                raise ValueError("sym_cut must be between or through")
            for a in self.planes:
                if self.sym_bc.get(a) not in ("pec", "pmc"):
                    raise ValueError(f"sym_bc[{a}] must be pec or pmc")
            if "z" in self.planes and self.sym_bc["z"] != "pec":
                raise ValueError("the z=0 plane must be PEC (rung currents pass straight through it)")

    def to_json(self, path: Path) -> None:
        d = asdict(self)
        Path(path).write_text(json.dumps(d, indent=2))

    @classmethod
    def from_json(cls, path: Path) -> "CoilDesign":
        d = json.loads(Path(path).read_text())
        ph = d.pop("phantom", None)
        des = cls(**d)
        des.phantom = Phantom(**ph) if ph else None
        return des
