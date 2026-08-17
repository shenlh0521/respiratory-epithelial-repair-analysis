#!/usr/bin/env python3
import os
import gzip
import shutil
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

PROJECT = Path(os.environ.get("PROJECT_ROOT", ".")).resolve()
DATA = PROJECT / "02_scRNA/SECOND_VALIDATION/GSE233145/data"
DERIVED = PROJECT / "02_scRNA/SECOND_VALIDATION/GSE233145/derived"
WORK = PROJECT / "work"

GENE_SETS = {
    "DNA_DAMAGE_RESPONSE": ["H2AX", "TP53BP1", "ATM", "ATR", "CHEK1", "CHEK2"],
    "OXIDATIVE_STRESS": ["NFE2L2", "HMOX1", "NQO1", "GCLC", "GCLM", "SOD2", "TXNRD1"],
    "ABNORMAL_REPAIR": ["KRT8", "KRT17", "KRT19", "CLDN4", "SFN", "LGALS3", "KRT14"],
    "KRT14_KRT17_REPAIR": ["KRT14", "KRT17"],
    "CILIA_CONSENSUS": ["FOXJ1", "MCIDAS", "GMNC", "MYB", "TP73", "RFX2", "RFX3", "CCNO",
        "CDC20B", "DEUP1", "TUBB4B", "DNAH5", "DNAH9", "DNAH11", "DNAI1", "DNAI2", "CCDC39",
        "CCDC40", "SPEF2", "HYDIN", "PIFO", "CFAP43", "CFAP44", "RSPH1", "RSPH4A", "RSPH9"],
    "MULTICILIOGENESIS": ["FOXJ1", "MCIDAS", "GMNC", "MYB", "TP73", "RFX2", "RFX3", "CCNO",
        "CDC20B", "DEUP1", "CEP78", "CETN2", "PLK4", "STIL"],
    "MATURE_CILIATED": ["FOXJ1", "PIFO", "TPPP3", "CAPS", "RSPH1"],
}
TARGETS = sorted(set(sum(GENE_SETS.values(), []) + ["KRT5", "TP63"] + [f"RFX{i}" for i in range(1, 8)]))


def decode_array(values):
    values = np.asarray(values)
    if values.dtype.kind in "SO":
        return np.array([x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x) for x in values])
    return values


def read_element(node):
    if isinstance(node, h5py.Dataset):
        return decode_array(node[()])
    if isinstance(node, h5py.Group) and "codes" in node and "categories" in node:
        codes = np.asarray(node["codes"][()])
        categories = decode_array(node["categories"][()])
        return np.array([categories[x] if x >= 0 else None for x in codes], dtype=object)
    raise ValueError(f"Unsupported H5AD element: {node.name}")


def read_frame(group):
    index_key = group.attrs.get("_index", "_index")
    if isinstance(index_key, bytes):
        index_key = index_key.decode()
    order = group.attrs.get("column-order", [k for k in group.keys() if k != index_key])
    order = [x.decode() if isinstance(x, bytes) else str(x) for x in order]
    frame = pd.DataFrame({name: read_element(group[name]) for name in order if name in group})
    if index_key in group:
        frame.index = decode_array(group[index_key][()])
    return frame


def read_sparse(group):
    shape = tuple(group.attrs["shape"])
    data = group["data"][()]
    indices = group["indices"][()]
    indptr = group["indptr"][()]
    encoding = group.attrs.get("encoding-type", "csr_matrix")
    if isinstance(encoding, bytes):
        encoding = encoding.decode()
    cls = sparse.csc_matrix if "csc" in encoding else sparse.csr_matrix
    return cls((data, indices, indptr), shape=shape)


def aggregate_file(gz_path, label):
    temp = WORK / f"{gz_path.name[:-3]}"
    print(f"Decompressing {gz_path.name} to temporary H5AD")
    with gzip.open(gz_path, "rb") as source, temp.open("wb") as target:
        shutil.copyfileobj(source, target, length=8 * 1024 * 1024)
    try:
        with h5py.File(temp, "r") as handle:
            print("H5AD keys:", list(handle.keys()))
            obs = read_frame(handle["obs"])
            source_group = handle["raw"] if "raw" in handle and "X" in handle["raw"] else handle
            var = read_frame(source_group["var"])
            var_names = np.array([str(x).upper() for x in var.index])
            selected = [(gene, int(np.where(var_names == gene)[0][0])) for gene in TARGETS if np.any(var_names == gene)]
            genes = [x[0] for x in selected]
            indices = [x[1] for x in selected]
            x_node = source_group["X"]
            if isinstance(x_node, h5py.Group):
                matrix = read_sparse(x_node)
                target_matrix = matrix[:, indices].tocsr()
                library_size = np.asarray(matrix.sum(axis=1)).ravel()
                is_counts = np.allclose(target_matrix.data, np.round(target_matrix.data))
            else:
                target_matrix = np.asarray(x_node[:, indices])
                library_size = np.asarray(x_node[:]).sum(axis=1)
                is_counts = np.allclose(target_matrix[np.isfinite(target_matrix)], np.round(target_matrix[np.isfinite(target_matrix)]))

            aliases = {
                "patient": ["patient", "donor", "Patient"],
                "health_state": ["health_state", "condition", "Health_state"],
                "time_point": ["time_point", "day", "time"],
                "cell_type": ["cell_type", "Cell_type", "annotation"],
            }
            meta = {}
            for wanted, choices in aliases.items():
                found = next((c for c in choices if c in obs.columns), None)
                if found is None:
                    raise KeyError(f"Missing {wanted}; available obs columns: {list(obs.columns)}")
                meta[wanted] = obs[found].astype(str).to_numpy()
            # The combined author H5ADs encode time as d0/d28, whereas the
            # companion metadata and downstream tables use "day 0"/"day 28".
            meta["time_point"] = np.repeat("day 0" if label == "day0" else "day 28", len(obs))
            frame = pd.DataFrame(meta)
            frame["row"] = np.arange(len(frame))
            frame["source_object"] = label

            groups = []
            for keys, part in frame.groupby(["patient", "health_state", "time_point", "cell_type", "source_object"], dropna=False):
                rows = part["row"].to_numpy()
                if sparse.issparse(target_matrix):
                    sums = np.asarray(target_matrix[rows].sum(axis=0)).ravel()
                    means = np.asarray(target_matrix[rows].mean(axis=0)).ravel()
                else:
                    sums = np.nansum(target_matrix[rows], axis=0)
                    means = np.nanmean(target_matrix[rows], axis=0)
                total = float(np.nansum(library_size[rows]))
                values = np.log2((sums + 0.5) / (total + 1) * 1e6) if is_counts else means
                for gene, value in zip(genes, values):
                    groups.append(dict(zip(["patient", "health_state", "time_point", "cell_type", "source_object"], keys),
                                       n_cells=len(rows), gene=gene, value=float(value),
                                       value_scale="log2_CPM_from_raw" if is_counts else "mean_author_normalized_expression"))

            # Add combined basal and all-epithelial summaries without changing author annotations.
            for state_name, mask in {
                "Basal_combined": np.array([str(x).startswith("Basal") for x in frame["cell_type"]], dtype=bool),
                "All_epithelial": np.ones(len(frame), dtype=bool),
            }.items():
                subset = frame[mask]
                for keys, part in subset.groupby(["patient", "health_state", "time_point", "source_object"], dropna=False):
                    rows = part["row"].to_numpy()
                    if sparse.issparse(target_matrix):
                        sums = np.asarray(target_matrix[rows].sum(axis=0)).ravel()
                        means = np.asarray(target_matrix[rows].mean(axis=0)).ravel()
                    else:
                        sums = np.nansum(target_matrix[rows], axis=0)
                        means = np.nanmean(target_matrix[rows], axis=0)
                    total = float(np.nansum(library_size[rows]))
                    values = np.log2((sums + 0.5) / (total + 1) * 1e6) if is_counts else means
                    for gene, value in zip(genes, values):
                        groups.append({"patient": keys[0], "health_state": keys[1], "time_point": keys[2],
                                       "cell_type": state_name, "source_object": keys[3], "n_cells": len(rows),
                                       "gene": gene, "value": float(value),
                                       "value_scale": "log2_CPM_from_raw" if is_counts else "mean_author_normalized_expression"})
            return pd.DataFrame(groups), {"source": label, "n_cells": len(obs), "n_genes": len(var_names),
                                          "n_targets": len(genes), "targets": ";".join(genes), "value_scale": groups[0]["value_scale"]}
    finally:
        temp.unlink(missing_ok=True)


DERIVED.mkdir(parents=True, exist_ok=True)
all_frames, qc = [], []
for filename, label in [
    ("GSE233145_basal_cells_nonCLD_COPD_day0.h5ad.gz", "day0"),
    ("GSE233145_basal_cells_nonCLD_COPD_day28.h5ad.gz", "day28"),
]:
    frame, info = aggregate_file(DATA / filename, label)
    all_frames.append(frame)
    qc.append(info)

pd.concat(all_frames, ignore_index=True).to_csv(DERIVED / "target_donor_state_pseudobulk.csv", index=False)
pd.DataFrame(qc).to_csv(DERIVED / "target_donor_state_pseudobulk_qc.csv", index=False)
print("Wrote", DERIVED / "target_donor_state_pseudobulk.csv")
