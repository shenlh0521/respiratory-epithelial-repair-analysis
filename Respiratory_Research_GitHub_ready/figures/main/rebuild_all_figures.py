#!/usr/bin/env python3
"""Rebuild Figures 1-8 from the Excel-grounded source workbook."""
from pathlib import Path
import subprocess, sys

here = Path(__file__).parent
subprocess.run([sys.executable, str(here / "code/rebuild_figures_from_excel.py")], check=True)
