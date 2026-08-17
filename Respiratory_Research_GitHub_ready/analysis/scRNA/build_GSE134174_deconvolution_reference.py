#!/usr/bin/env python3
"""Build a donor-aware airway epithelial reference from existing GSE134174 files.

No external data are downloaded. The output is intentionally marker-restricted and
is used only for exploratory cross-platform composition projection.
"""
import os
from pathlib import Path
import gzip
import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("PROJECT_ROOT", ".")).resolve()
BASE = ROOT / "02_scRNA/GSE134174"
DATA = BASE / "data"
OUT = BASE / "third_stage_results/deconvolution_reference"
OUT.mkdir(parents=True, exist_ok=True)

MARKERS = {
    "Basal_repair": ["KRT5", "TP63", "KRT14", "KRT17", "KRT15", "NGFR"],
    "Ciliated": ["FOXJ1", "PIFO", "TPPP3", "CAPS", "RSPH1", "TUBB4B", "DNAH5"],
    "Secretory": ["SCGB1A1", "SCGB3A1", "MUC5AC", "MUC5B", "SPDEF", "AGR2", "BPIFB1", "WFDC2"],
    "Transitional": ["KRT4", "KRT13", "KRT8", "CLDN4", "SFN", "LGALS3", "KRT19"],
}
TARGET = set(sum(MARKERS.values(), []))

def broad_state(x):
    if x == "Ciliated": return "Ciliated"
    if x in {"Differentiating.basal", "Proliferating.basal", "Proteasomal.basal", "SMG.basal"}: return "Basal_repair"
    if x == "KRT8.high": return "Transitional"
    if x in {"Mucus.secretory", "SMG.secretory"}: return "Secretory"
    return None

md = pd.read_csv(DATA / "GSE134174_Processed_invivo_metadata.txt.gz", sep="\t")
md["broad_state"] = md["cluster_ident"].map(broad_state)
md = md[(md["Smoke_status"].str.lower() == "never") & md["broad_state"].notna()].copy()
md["group"] = md["Donor"].astype(str) + "|" + md["broad_state"]

raw_path = DATA / "GSE134174_Processed_invivo_raw.txt.gz"
with gzip.open(raw_path, "rt") as fh:
    cells = fh.readline().rstrip("\n").split("\t")
cell_to_group = dict(zip(md["Cell"], md["group"]))
groups = sorted(md["group"].unique())
gidx = {g:i for i,g in enumerate(groups)}
selected = [(i, gidx[cell_to_group[c]]) for i,c in enumerate(cells) if c in cell_to_group]

counts = {}
with gzip.open(raw_path, "rt") as fh:
    fh.readline()
    for line in fh:
        fields = line.rstrip("\n").split("\t")
        gene = fields[0].upper()
        if gene not in TARGET:
            continue
        vals = np.zeros(len(groups), dtype=float)
        for ci, gi in selected:
            vals[gi] += float(fields[ci + 1])
        counts[gene] = vals

mat = pd.DataFrame.from_dict(counts, orient="index", columns=groups).sort_index()
libs = mat.sum(axis=0).replace(0, np.nan)
logcpm = np.log2(mat.divide(libs, axis=1) * 1e6 + 1)
long = logcpm.T.reset_index(names="group")
long[["donor", "state"]] = long["group"].str.split("|", expand=True)
reference = long.drop(columns=["group", "donor"]).groupby("state").mean().T
reference = reference.reindex(columns=["Basal_repair", "Ciliated", "Secretory", "Transitional"])

cell_counts = md.groupby(["Donor", "broad_state"], observed=True).size().rename("n_cells").reset_index()
coverage = []
for state, genes in MARKERS.items():
    observed = sorted(set(genes) & set(reference.index))
    coverage.append({"state": state, "requested_n": len(genes), "observed_n": len(observed),
                     "observed_genes": ";".join(observed),
                     "missing_genes": ";".join(sorted(set(genes)-set(observed)))})

reference.reset_index(names="gene").to_csv(OUT / "GSE134174_never_smoker_reference_logCPM.csv", index=False)
logcpm.T.reset_index(names="donor_state").to_csv(OUT / "GSE134174_donor_state_marker_logCPM.csv", index=False)
cell_counts.to_csv(OUT / "GSE134174_reference_cell_counts.csv", index=False)
pd.DataFrame(coverage).to_csv(OUT / "GSE134174_reference_marker_coverage.csv", index=False)

with open(OUT / "README.md", "w") as fh:
    fh.write("# GSE134174 deconvolution reference\n\n")
    fh.write(f"Reference restricted to never-smoker donors and {len(reference)} predefined epithelial marker genes. ")
    fh.write("Rare PNEC/ionocyte/tuft states were excluded. Broad states are donor-pseudobulk averages; ")
    fh.write("the reference is for exploratory cross-platform projection, not absolute cell-fraction estimation.\n")
print(f"Wrote {len(reference)} genes across {len(groups)} donor-state pseudobulks from {len(md)} cells")
