#!/usr/bin/env python3
"""Create a non-overwriting ALI donor-day pseudobulk with third-stage targets."""
from __future__ import annotations
import os
import csv, gzip
from pathlib import Path
import numpy as np

PROJECT = Path(os.environ.get("PROJECT_ROOT", ".")).resolve()
ROOT = PROJECT / "02_scRNA" / "GSE134174"
DATA = ROOT / "data"
OUT = ROOT / "derived" / "invitro_target_pseudobulk_third_stage.csv"
TARGETS = {
    "H2AX","TP53BP1","ATM","ATR","CHEK1","CHEK2","KRT8","KRT17","KRT19","CLDN4","SFN","LGALS3","KRT14",
    "FOXJ1","MCIDAS","GMNC","MYB","TP73","RFX2","RFX3","CCNO","CDC20B","DEUP1","TUBB4B","DNAH5","DNAH9",
    "DNAH11","DNAI1","DNAI2","CCDC39","CCDC40","SPEF2","HYDIN","PIFO","CFAP43","CFAP44","RSPH1","RSPH4A",
    "RSPH9","CEP78","CETN2","PLK4","STIL","TPPP3","CAPS","KRT4","KRT13"
}

meta = {}
with gzip.open(DATA / "GSE134174_Processed_invitro_metadata.txt.gz", "rt", newline="") as fh:
    reader = csv.reader(fh, delimiter="\t"); next(reader)
    for fields in reader:
        if len(fields) == 4:
            cell, day, donor, cluster = fields
            meta[cell] = (donor, day, cluster)

with gzip.open(DATA / "GSE134174_Processed_invitro_raw.txt.gz", "rt") as fh:
    cells = fh.readline().rstrip("\n").split("\t")
    groups = sorted(set(meta[c] for c in cells)); gid = {g:i for i,g in enumerate(groups)}
    cell_group = np.fromiter((gid[meta[c]] for c in cells), dtype=np.int32)
    n_cells = np.bincount(cell_group, minlength=len(groups)).astype(np.int64)
    library = np.zeros(len(groups), dtype=np.float64); retained = {}; detected = {}; n_genes = 0
    for line in fh:
        gene, sep, values_text = line.partition("\t")
        if not sep: continue
        values = np.fromstring(values_text, sep="\t", dtype=np.float64)
        sums = np.bincount(cell_group, weights=values, minlength=len(groups)); library += sums
        symbol = gene.upper()
        if symbol in TARGETS:
            retained[symbol] = retained.get(symbol, 0) + sums
            detected[symbol] = detected.get(symbol, 0) + np.bincount(cell_group, weights=(values > 0), minlength=len(groups))
        n_genes += 1

with OUT.open("w", newline="") as fh:
    writer = csv.writer(fh)
    writer.writerow(["donor","day","cluster_ident","n_cells","library_size","gene","count","detected_cells","log2_cpm"])
    for gene in sorted(retained):
        for i, group in enumerate(groups):
            count = float(retained[gene][i]); log2_cpm = np.log2((count + .5)/(library[i] + 1)*1_000_000)
            writer.writerow(list(group)+[int(n_cells[i]),int(library[i]),gene,int(count),int(detected[gene][i]),float(log2_cpm)])

qc = OUT.with_suffix(".qc.txt")
qc.write_text(f"cells={len(cells)}\ngroups={len(groups)}\ngenes_streamed={n_genes}\ntargets_detected={len(retained)}\nmissing_targets={';'.join(sorted(TARGETS-set(retained)))}\n")
print(f"wrote {OUT} with {len(retained)} targets across {len(groups)} donor-day-cluster groups")
