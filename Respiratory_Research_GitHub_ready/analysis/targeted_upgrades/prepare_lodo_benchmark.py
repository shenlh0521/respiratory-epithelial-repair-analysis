#!/usr/bin/env python3
"""Prepare frozen-marker LODO mixtures and GSE320381 marker expression."""
import os
from pathlib import Path
import gzip
import numpy as np
import pandas as pd

PROJECT = Path(os.environ.get("PROJECT_ROOT", ".")).resolve()
OUT = PROJECT / "12_targeted_manuscript_upgrades/01_deconvolution_robustness"
SC = PROJECT / "02_scRNA/GSE134174/data"
SEED = 20260813
STATES = ["Basal_repair", "Ciliated", "Secretory", "Transitional"]
MARKERS = {
    "Basal_repair": ["KRT5", "TP63", "KRT14", "KRT17", "KRT15", "NGFR"],
    "Ciliated": ["FOXJ1", "PIFO", "TPPP3", "CAPS", "RSPH1", "TUBB4B", "DNAH5"],
    "Secretory": ["SCGB1A1", "SCGB3A1", "MUC5AC", "MUC5B", "SPDEF", "AGR2", "BPIFB1", "WFDC2"],
    "Transitional": ["KRT4", "KRT13", "KRT8", "CLDN4", "SFN", "LGALS3", "KRT19"],
}
TARGET = sorted(set(sum(MARKERS.values(), [])))

def broad_state(value):
    if value == "Ciliated": return "Ciliated"
    if value in {"Differentiating.basal", "Proliferating.basal", "Proteasomal.basal", "SMG.basal"}: return "Basal_repair"
    if value == "KRT8.high": return "Transitional"
    if value in {"Mucus.secretory", "SMG.secretory"}: return "Secretory"
    return None

def marker_matrix(path, wanted_cells):
    with gzip.open(path, "rt") as handle:
        all_cells = handle.readline().rstrip("\n").split("\t")
        selected_idx = [i for i, cell in enumerate(all_cells) if cell in wanted_cells]
        selected_cells = [all_cells[i] for i in selected_idx]
        rows = {}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            gene = fields[0].upper()
            if gene in TARGET and gene not in rows:
                rows[gene] = np.asarray([float(fields[i + 1]) for i in selected_idx])
    return selected_cells, pd.DataFrame.from_dict(rows, orient="index", columns=selected_cells).reindex(TARGET).fillna(0.0)

def cpm(vector):
    total = float(vector.sum())
    return vector / total * 1e6 if total > 0 else np.zeros_like(vector)

def prepare_lodo():
    md = pd.read_csv(SC / "GSE134174_Processed_invivo_metadata.txt.gz", sep="\t")
    md["broad_state"] = md["cluster_ident"].map(broad_state)
    md = md[(md["Smoke_status"].str.lower() == "never") & md["broad_state"].notna()].copy()
    cells, counts = marker_matrix(SC / "GSE134174_Processed_invivo_raw.txt.gz", set(md["Cell"]))
    md = md.set_index("Cell").loc[cells].reset_index()
    cell_index = {cell: i for i, cell in enumerate(cells)}
    donors = sorted(md["Donor"].unique())
    rng = np.random.default_rng(SEED)
    signatures, mixtures, truth, count_rows = [], [], [], []
    for donor in donors:
        held_in = [d for d in donors if d != donor]
        profiles = {}
        for ref_donor in held_in:
            for state in STATES:
                ids = md.loc[(md.Donor == ref_donor) & (md.broad_state == state), "Cell"]
                idx = [cell_index[x] for x in ids]
                profiles[(ref_donor, state)] = cpm(counts.iloc[:, idx].sum(axis=1).to_numpy(float))
        signature = np.column_stack([np.mean([profiles[(d, state)] for d in held_in], axis=0) for state in STATES])
        for gi, gene in enumerate(counts.index):
            for si, state in enumerate(STATES):
                value = signature[gi, si]
                signatures.append({"held_out_donor": donor, "gene": gene, "state": state, "cpm": value, "log2_cpm": np.log2(value + 1)})
        state_cells = {state: [cell_index[x] for x in md.loc[(md.Donor == donor) & (md.broad_state == state), "Cell"]] for state in STATES}
        for state in STATES: count_rows.append({"donor": donor, "state": state, "n_cells": len(state_cells[state])})
        for mix_number in range(1, 41):
            allocated = rng.multinomial(500, rng.dirichlet(np.repeat(1.5, 4)))
            for i in np.where(allocated == 0)[0]:
                j = int(np.argmax(allocated)); allocated[j] -= 1; allocated[i] = 1
            sampled = []
            for state, n in zip(STATES, allocated): sampled.extend(rng.choice(state_cells[state], size=int(n), replace=True))
            mix_cpm = cpm(counts.iloc[:, sampled].sum(axis=1).to_numpy(float))
            mix_id = f"{donor}_M{mix_number:02d}"
            for state, n in zip(STATES, allocated): truth.append({"mixture_id": mix_id, "held_out_donor": donor, "state": state, "observed": n / 500})
            for gi, gene in enumerate(counts.index): mixtures.append({"mixture_id": mix_id, "held_out_donor": donor, "gene": gene, "cpm": mix_cpm[gi], "log2_cpm": np.log2(mix_cpm[gi] + 1)})
    pd.DataFrame(signatures).to_csv(OUT / "lodo_reference_signatures.csv", index=False)
    pd.DataFrame(mixtures).to_csv(OUT / "lodo_mixture_expression.csv", index=False)
    pd.DataFrame(truth).to_csv(OUT / "lodo_mixture_truth.csv", index=False)
    pd.DataFrame(count_rows).to_csv(OUT / "reference_donor_state_counts.csv", index=False)

def prepare_gse320381():
    ref = pd.read_csv(PROJECT / "02_scRNA/GSE134174/third_stage_results/deconvolution_reference/GSE134174_never_smoker_reference_logCPM.csv")
    mapping = pd.read_csv(PROJECT / "files-mentioned-by-the-user-users/outputs/lung_bulla_cilia_cancer/data/metadata/target_symbol_ensembl.csv")
    extra = pd.DataFrame({"symbol": ["KRT4", "KRT13"], "ensembl_id": ["ENSG00000170477", "ENSG00000171401"]})
    mapping = pd.concat([mapping, extra]).drop_duplicates("ensembl_id")
    id2sym = dict(zip(mapping.ensembl_id, mapping.symbol.str.upper()))
    rows = {}
    with gzip.open(PROJECT / "01_PML_extended/GSE320381/raw/GSE320381_BX_counts.tsv.gz", "rt") as handle:
        sample_ids = handle.readline().rstrip("\n").split("\t")[1:]
        for line in handle:
            f = line.rstrip("\n").split("\t"); symbol = id2sym.get(f[0].split(".")[0])
            if symbol in set(ref.gene) and symbol not in rows: rows[symbol] = np.asarray(f[1:], dtype=float)
    bulk = pd.DataFrame.from_dict(rows, orient="index", columns=sample_ids).sort_index()
    bulk_cpm = bulk.divide(bulk.sum(axis=0), axis=1) * 1e6
    long = bulk_cpm.T.reset_index(names="sample_id").melt(id_vars="sample_id", var_name="gene", value_name="cpm")
    long["log2_cpm"] = np.log2(long.cpm + 1)
    long.to_csv(OUT / "GSE320381_marker_expression.csv", index=False)

if __name__ == "__main__":
    prepare_lodo(); prepare_gse320381()
    print("Prepared 240 LODO mixtures and GSE320381 marker expression")
