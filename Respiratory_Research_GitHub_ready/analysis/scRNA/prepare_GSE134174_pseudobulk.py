#!/usr/bin/env python3
"""Stream GEO dense count matrices into compact donor/state pseudobulks.

The full matrices are never materialized in memory. Library sizes are calculated
from every gene, while only prespecified project genes are retained.
"""
from __future__ import annotations
import os

import csv
import gzip
import argparse
from pathlib import Path
import numpy as np

PROJECT = Path(os.environ.get("PROJECT_ROOT", ".")).resolve()
ROOT = PROJECT / "02_scRNA" / "GSE134174"
DATA = ROOT / "data"
OUT = ROOT / "derived"
OUT.mkdir(parents=True, exist_ok=True)

GENE_SETS = {
    "DNA_DAMAGE_RESPONSE": ["H2AX", "TP53BP1", "ATM", "ATR", "CHEK1", "CHEK2"],
    "OXIDATIVE_STRESS": ["NFE2L2", "HMOX1", "NQO1", "GCLC", "GCLM", "SOD2", "TXNRD1"],
    "ABNORMAL_REPAIR": ["KRT8", "KRT17", "KRT19", "CLDN4", "SFN", "LGALS3", "KRT14"],
    "KRT14_KRT17_REPAIR": ["KRT14", "KRT17"],
    "CILIA_CONSENSUS": ["FOXJ1", "MCIDAS", "GMNC", "MYB", "TP73", "RFX2", "RFX3", "CCNO",
                         "CDC20B", "DEUP1", "TUBB4B", "DNAH5", "DNAH9", "DNAH11", "DNAI1",
                         "DNAI2", "CCDC39", "CCDC40", "SPEF2", "HYDIN", "PIFO", "CFAP43",
                         "CFAP44", "RSPH1", "RSPH4A", "RSPH9"],
    "MULTICILIOGENESIS": ["FOXJ1", "MCIDAS", "GMNC", "MYB", "TP73", "RFX2", "RFX3", "CCNO",
                           "CDC20B", "DEUP1", "CEP78", "CETN2", "PLK4", "STIL"],
    "MATURE_CILIATED": ["FOXJ1", "PIFO", "TPPP3", "CAPS", "RSPH1"],
    "BASAL_MARKERS": ["KRT5", "TP63", "KRT14", "KRT17"],
}
TARGETS = set(sum(GENE_SETS.values(), [])) | {"KRT14", "KRT17", "FOXJ1", "MCIDAS", "GMNC"}


def load_invivo_metadata(path: Path):
    with gzip.open(path, "rt", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return {r["Cell"]: r for r in rows}


def load_invitro_metadata(path: Path):
    result = {}
    with gzip.open(path, "rt", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader)
        for fields in reader:
            if len(fields) != 4:
                continue
            cell, day, donor, cluster = fields
            result[cell] = {"Cell": cell, "day": day, "donor": donor, "clust_ident": cluster}
    return result


def stream_pseudobulk(matrix_path: Path, metadata: dict, mode: str, out_path: Path):
    with gzip.open(matrix_path, "rt") as fh:
        cells = fh.readline().rstrip("\n").split("\t")
        missing = [c for c in cells if c not in metadata]
        if missing:
            raise RuntimeError(f"{len(missing)} matrix cells absent from metadata; first={missing[0]}")

        if mode == "invivo":
            group_tuples = [
                (metadata[c]["Donor"], metadata[c]["Smoke_status"],
                 metadata[c]["cluster_ident"], metadata[c]["subcluster_ident"])
                for c in cells
            ]
            group_cols = ["donor", "smoke_status", "cluster_ident", "subcluster_ident"]
        else:
            group_tuples = [
                (metadata[c]["donor"], metadata[c]["day"], metadata[c]["clust_ident"])
                for c in cells
            ]
            group_cols = ["donor", "day", "cluster_ident"]

        groups = sorted(set(group_tuples))
        group_id = {g: i for i, g in enumerate(groups)}
        cell_group = np.fromiter((group_id[g] for g in group_tuples), dtype=np.int32)
        n_groups = len(groups)
        n_cells = np.bincount(cell_group, minlength=n_groups).astype(np.int64)
        library = np.zeros(n_groups, dtype=np.float64)
        retained = {}
        detected = {}
        n_genes = 0

        for line in fh:
            gene, sep, values_text = line.partition("\t")
            if not sep:
                continue
            values = np.fromstring(values_text, sep="\t", dtype=np.float64)
            if values.size != len(cells):
                raise RuntimeError(f"Malformed row {gene}: {values.size} values for {len(cells)} cells")
            sums = np.bincount(cell_group, weights=values, minlength=n_groups)
            library += sums
            gene = gene.upper()
            if gene in TARGETS:
                if gene in retained:
                    retained[gene] += sums
                    detected[gene] += np.bincount(cell_group, weights=(values > 0), minlength=n_groups)
                else:
                    retained[gene] = sums
                    detected[gene] = np.bincount(cell_group, weights=(values > 0), minlength=n_groups)
            n_genes += 1

    header = group_cols + ["n_cells", "library_size", "gene", "count", "detected_cells", "log2_cpm"]
    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        for gene in sorted(retained):
            for i, group in enumerate(groups):
                count = float(retained[gene][i])
                log2_cpm = np.log2((count + 0.5) / (library[i] + 1.0) * 1_000_000.0)
                writer.writerow(list(group) + [int(n_cells[i]), int(library[i]), gene, int(count),
                                               int(detected[gene][i]), float(log2_cpm)])

    qc_path = out_path.with_suffix(".qc.txt")
    with qc_path.open("w") as fh:
        fh.write(f"mode={mode}\n")
        fh.write(f"matrix={matrix_path}\n")
        fh.write(f"cells={len(cells)}\n")
        fh.write(f"groups={n_groups}\n")
        fh.write(f"genes_streamed={n_genes}\n")
        fh.write(f"targets_detected={len(retained)}\n")
        fh.write("missing_targets=" + ";".join(sorted(TARGETS - set(retained))) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["all", "invivo", "invitro"], default="all")
    args = parser.parse_args()
    invivo_matrix = DATA / "GSE134174_Processed_invivo_raw.txt.gz"
    invitro_matrix = DATA / "GSE134174_Processed_invitro_raw.txt.gz"
    required_files = []
    if args.mode in ("all", "invivo"):
        required_files += [invivo_matrix, DATA / "GSE134174_Processed_invivo_metadata.txt.gz"]
    if args.mode in ("all", "invitro"):
        required_files += [invitro_matrix, DATA / "GSE134174_Processed_invitro_metadata.txt.gz"]
    for required in required_files:
        if not required.exists():
            raise FileNotFoundError(required)
    if args.mode in ("all", "invivo"):
        stream_pseudobulk(invivo_matrix,
                          load_invivo_metadata(DATA / "GSE134174_Processed_invivo_metadata.txt.gz"),
                          "invivo", OUT / "invivo_target_pseudobulk.csv")
    if args.mode in ("all", "invitro"):
        stream_pseudobulk(invitro_matrix,
                          load_invitro_metadata(DATA / "GSE134174_Processed_invitro_metadata.txt.gz"),
                          "invitro", OUT / "invitro_target_pseudobulk.csv")


if __name__ == "__main__":
    main()
