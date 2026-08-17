#!/usr/bin/env python3
"""Stream GSE134174 ALI counts into a compact ciliated-lineage matrix."""
import os
from pathlib import Path
import csv
import gzip
import numpy as np
import pandas as pd

PROJECT = Path(os.environ.get("PROJECT_ROOT", ".")).resolve()
DATA = PROJECT / "02_scRNA/GSE134174/data"
OUT = PROJECT / "12_targeted_manuscript_upgrades/02_ALI_trajectory"
LINEAGE = {
    "basal.colonies", "basal.subconfluent", "basal.confluent", "p.basal",
    "d.basal", "ciliating.early", "ciliating.late", "ciliated",
}
SIGNATURES = {
    "ABNORMAL_REPAIR": ["KRT8", "KRT17", "KRT19", "CLDN4", "SFN", "LGALS3", "KRT14"],
    "KRT14_KRT17_REPAIR": ["KRT14", "KRT17"],
    "CILIA_CONSENSUS": ["FOXJ1", "MCIDAS", "GMNC", "MYB", "TP73", "RFX2", "RFX3", "CCNO", "CDC20B", "DEUP1", "TUBB4B", "DNAH5", "DNAH9", "DNAH11", "DNAI1", "DNAI2", "CCDC39", "CCDC40", "SPEF2", "HYDIN", "PIFO", "CFAP43", "CFAP44", "RSPH1", "RSPH4A", "RSPH9"],
    "MULTICILIOGENESIS": ["FOXJ1", "MCIDAS", "GMNC", "MYB", "TP73", "RFX2", "RFX3", "CCNO", "CDC20B", "DEUP1", "CEP78", "CETN2", "PLK4", "STIL"],
    "MATURE_CILIATED": ["FOXJ1", "PIFO", "TPPP3", "CAPS", "RSPH1"],
    "KRT4_KRT13_TRANSITIONAL_FATE": ["KRT4", "KRT13"],
}
CORE = set(sum(SIGNATURES.values(), [])) | {"FOXJ1", "MCIDAS", "GMNC"}

def metadata():
    rows = []
    with gzip.open(DATA / "GSE134174_Processed_invitro_metadata.txt.gz", "rt", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t"); next(reader)
        for fields in reader:
            if len(fields) == 4:
                cell, day, donor, cluster = fields
                if cluster in LINEAGE: rows.append((cell, day, donor, cluster))
    return pd.DataFrame(rows, columns=["cell", "day", "donor", "cluster"])

def main():
    md = metadata()
    matrix_path = DATA / "GSE134174_Processed_invitro_raw.txt.gz"
    with gzip.open(matrix_path, "rt") as handle:
        all_cells = handle.readline().rstrip("\n").split("\t")
        idx = [i for i, x in enumerate(all_cells) if x in set(md.cell)]
        cells = [all_cells[i] for i in idx]
        library = np.zeros(len(cells))
        for line in handle:
            f = line.rstrip("\n").split("\t")
            library += np.asarray([float(f[i + 1]) for i in idx])
    variances, detection = {}, {}
    with gzip.open(matrix_path, "rt") as handle:
        handle.readline()
        for line in handle:
            f = line.rstrip("\n").split("\t"); gene = f[0].upper()
            values = np.asarray([float(f[i + 1]) for i in idx])
            logcpm = np.log2(values / np.maximum(library, 1) * 1e6 + 1)
            variances[gene] = float(np.var(logcpm, ddof=1)); detection[gene] = int(np.sum(values > 0))
    candidates = [g for g in variances if detection[g] >= max(10, int(0.01 * len(cells)))]
    hvg = sorted(candidates, key=lambda g: (-variances[g], g))[:1000]
    keep = set(hvg) | CORE
    expression = {}
    with gzip.open(matrix_path, "rt") as handle:
        handle.readline()
        for line in handle:
            f = line.rstrip("\n").split("\t"); gene = f[0].upper()
            if gene not in keep: continue
            values = np.asarray([float(f[i + 1]) for i in idx])
            expression[gene] = np.log2(values / np.maximum(library, 1) * 1e6 + 1)
    expr = pd.DataFrame.from_dict(expression, orient="index", columns=cells)
    md = md.set_index("cell").loc[cells].reset_index()
    md["day_numeric"] = md.day.map(lambda x: -3 if x == "seed_day" else (-2 if x == "day_minus2" else int(x.replace("day_", ""))))
    md["library_size"] = library.astype(int)
    md.to_csv(OUT / "ALI_lineage_metadata.csv", index=False)
    expr.reset_index(names="gene").to_csv(OUT / "ALI_lineage_log2CPM.csv", index=False)
    pd.DataFrame({"gene": hvg, "variance": [variances[g] for g in hvg], "detected_cells": [detection[g] for g in hvg]}).to_csv(OUT / "ALI_HVG_1000.csv", index=False)
    pd.DataFrame([{"signature": k, "defined_genes": ";".join(v), "detected_genes": ";".join([g for g in v if g in expr.index]), "n_detected": sum(g in expr.index for g in v)} for k, v in SIGNATURES.items()]).to_csv(OUT / "trajectory_signature_coverage.csv", index=False)
    print(f"Prepared {expr.shape[0]} genes x {expr.shape[1]} cells across {md.donor.nunique()} donors")

if __name__ == "__main__": main()
