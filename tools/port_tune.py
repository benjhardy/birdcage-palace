"""Port tuning with circuitRF. The code now lives in birdcage_palace/porttune.py (installed with the package, as
`birdcage porttune ...`); this script keeps the old entry point:

    python tools/port_tune.py assemble  --cases runs/pt_e --out runs/pt_e/birdcage16.s16p --rotation-average
    python tools/port_tune.py resonance --snp runs/pt_e/birdcage16.s16p
    python tools/port_tune.py tune      --snp runs/pt_e/birdcage16.s16p --name unloaded --esr 0.68 [--crf PATH] [--wsl-distro NAME]
    python tools/port_tune.py verify    --snp ... --ref-case ... --tuned-case ...
    python tools/port_tune.py check     --cases runs/ptc_e --full runs/ptc_full

Step-by-step guide: docs/circuitrf-tuning.md.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from birdcage_palace.porttune import *  # noqa: E402,F401,F403
from birdcage_palace.porttune import main  # noqa: E402

if __name__ == "__main__":
    main()
