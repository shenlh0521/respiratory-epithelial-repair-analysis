#!/usr/bin/env python3
"""One-command entry point for the Excel-grounded R figure builder."""
from pathlib import Path
import subprocess
import pandas as pd

script = Path(__file__).with_suffix(".R")
base = script.parent.parent
workbook = base / "source_data/Figure_1-8_Source_Data.xlsx"
extract = base / "qa/extracted_from_workbook_for_plot"
extract.mkdir(parents=True, exist_ok=True)
for sheet in ["Fig1_design", "Fig2_discovery", "Fig3_severity", "Fig4_persistence", "Fig5_scRNA", "Fig6_ALI_trajectory", "Fig7_maturation", "Fig8_exposure"]:
    pd.read_excel(workbook, sheet_name=sheet).to_csv(extract / f"{sheet}.csv", index=False)
subprocess.run(["Rscript", str(script)], check=True)
