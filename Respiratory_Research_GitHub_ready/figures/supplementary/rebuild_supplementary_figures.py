#!/usr/bin/env python3
"""One-command builder for Respiratory Research Supplementary Figures S1-S25."""
from pathlib import Path
import subprocess

script = Path(__file__).parent / "code/rebuild_supplementary_figures.R"
subprocess.run(["Rscript", str(script)], check=True)
