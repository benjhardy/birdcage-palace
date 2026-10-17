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
    c_leg: float | None = None    # F, capacitor in each rung (lowpass/bandpass); None = analytic first guess
    c_ring: float | None = None   # F, capacitor in each end-ring segment (highpass/bandpass)
    bandpass_ratio: float = 1.0   # bandpass only: c_ring / c_leg for the first guess
    phantom: Phantom | None = None
    # mesh
    h_coil: float = 0.005         # m, target edge length on the conductors
    h_gap: float = 0.0025         # m, edge length at the capacitor gaps
    h_max: float = 0.035          # m, far-field edge length
    order: int = 2                # FE order in Palace
    notes: dict = field(default_factory=dict)

    @property
    def radius(self) -> float:
        return self.diameter / 2

    @property
    def f_larmor(self) -> float:
        return GAMMA_BAR_1H * self.b0

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
