#!/usr/bin/env python3
"""Exploratory scRNA-reference projection and composition-adjusted PML severity."""
import os
from pathlib import Path
import gzip
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.stats import t as student_t
import matplotlib.pyplot as plt

ROOT = Path(os.environ.get("PROJECT_ROOT", ".")).resolve()
BASE = ROOT / "01_PML_extended/GSE320381"
OUT = BASE / "composition_adjusted_results"
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(exist_ok=True)
REFPATH = ROOT / "02_scRNA/GSE134174/third_stage_results/deconvolution_reference/GSE134174_never_smoker_reference_logCPM.csv"
COUNT = BASE / "raw/GSE320381_BX_counts.tsv.gz"
SCORES = BASE / "results/biopsy_metadata_scores.csv"
MAP = ROOT / "files-mentioned-by-the-user-users/outputs/lung_bulla_cilia_cancer/data/metadata/target_symbol_ensembl.csv"

ref = pd.read_csv(REFPATH).set_index("gene")
states = list(ref.columns)
mp = pd.read_csv(MAP)
extra = pd.DataFrame({"symbol":["KRT4","KRT13"], "ensembl_id":["ENSG00000170477","ENSG00000171401"]})
mp = pd.concat([mp, extra], ignore_index=True).drop_duplicates("ensembl_id")
id2sym = dict(zip(mp["ensembl_id"], mp["symbol"].str.upper()))

bulk_rows = {}
with gzip.open(COUNT, "rt") as fh:
    header = fh.readline().rstrip("\n").split("\t")
    sample_ids = header[1:]
    for line in fh:
        f = line.rstrip("\n").split("\t")
        ens = f[0].split(".")[0]
        sym = id2sym.get(ens)
        if sym in ref.index and sym not in bulk_rows:
            bulk_rows[sym] = np.asarray(f[1:], dtype=float)
bulk = pd.DataFrame.from_dict(bulk_rows, orient="index", columns=sample_ids)
genes = sorted(set(ref.index) & set(bulk.index))
ref = ref.loc[genes]
bulk = bulk.loc[genes]
bulk_log = np.log2(bulk.divide(bulk.sum(axis=0), axis=1) * 1e6 + 1)

# Gene-wise 0-1 scaling preserves within-gene state/sample contrasts while reducing
# cross-platform intensity differences. Fractions are consequently relative estimates.
def row_minmax(x):
    lo = x.min(axis=1); hi = x.max(axis=1); den = (hi-lo).replace(0, np.nan)
    return x.sub(lo, axis=0).divide(den, axis=0).fillna(0)
A = row_minmax(ref).to_numpy()
Y = row_minmax(bulk_log).to_numpy()
fractions = []
for j, sid in enumerate(sample_ids):
    coef, resid = nnls(A, Y[:, j])
    coef = coef / coef.sum() if coef.sum() else np.repeat(np.nan, len(coef))
    row = {"sample_id":sid, "nnls_residual":resid, "n_shared_markers":len(genes)}
    row.update(dict(zip(states, coef)))
    fractions.append(row)
fractions = pd.DataFrame(fractions)
fractions.to_csv(OUT / "GSE320381_relative_epithelial_composition.csv", index=False)

d = pd.read_csv(SCORES).merge(fractions, on="sample_id", how="inner")
d["imbalance_z"] = (d["Repair_Cilia_Imbalance"]-d["Repair_Cilia_Imbalance"].mean())/d["Repair_Cilia_Imbalance"].std(ddof=1)

def cluster_model(label, covariates):
    z = d.dropna(subset=["imbalance_z", "grade_ordinal", "patient_id"]+covariates).copy()
    cols = ["intercept", "grade_ordinal"] + covariates
    X = np.column_stack([np.ones(len(z)), z[["grade_ordinal"]+covariates].to_numpy(float)])
    y = z["imbalance_z"].to_numpy(float)
    bread = np.linalg.inv(X.T @ X)
    beta = bread @ X.T @ y
    resid = y - X @ beta
    meat = np.zeros((X.shape[1], X.shape[1]))
    groups = z["patient_id"].astype(str).to_numpy()
    for group in np.unique(groups):
        idx = np.where(groups == group)[0]
        u = X[idx].T @ resid[idx]
        meat += np.outer(u, u)
    G, N, P = len(np.unique(groups)), len(z), X.shape[1]
    vcov = (G/(G-1))*((N-1)/(N-P)) * bread @ meat @ bread
    est = beta[1]; se = np.sqrt(vcov[1,1]); stat = est/se
    pval = 2*student_t.sf(abs(stat), df=G-1)
    crit = student_t.ppf(.975, df=G-1)
    return {"model":label, "term":"grade_ordinal", "effect_size_SD_per_grade":est,
            "std_error":se, "CI95_lower":est-crit*se, "CI95_upper":est+crit*se,
            "P":pval, "n_samples":len(z), "n_subjects":G,
            "analysis_method":"OLS with patient-cluster robust SE"}

# Transitional is the omitted compositional reference category.
unadj = cluster_model("unadjusted", [])
adj = cluster_model("composition_adjusted", ["Basal_repair", "Ciliated", "Secretory"])
models = pd.DataFrame([unadj, adj])
# The two rows are planned sensitivity models of one estimand; retain BH explicitly.
order = np.argsort(models["P"].to_numpy())
q = np.empty(len(models)); ranked = models["P"].to_numpy()[order] * len(models) / np.arange(1, len(models)+1)
ranked = np.minimum.accumulate(ranked[::-1])[::-1]
q[order] = np.minimum(ranked, 1)
models["FDR"] = q
models["attenuation_percent_vs_unadjusted"] = 100*(1-models["effect_size_SD_per_grade"]/unadj["effect_size_SD_per_grade"])
models.to_csv(OUT / "composition_adjusted_severity_models.csv", index=False)

summ = d.groupby("grade_ordinal")[states].agg(["mean","count"])
summ.to_csv(OUT / "composition_by_histology_grade.csv")
pd.DataFrame({"gene":genes}).to_csv(OUT / "shared_deconvolution_markers.csv", index=False)

fig, ax = plt.subplots(figsize=(7.2,4.4))
y = np.arange(2)
ax.errorbar(models["effect_size_SD_per_grade"], y,
            xerr=[models["effect_size_SD_per_grade"]-models["CI95_lower"], models["CI95_upper"]-models["effect_size_SD_per_grade"]],
            fmt="o", color="#263238", capsize=3)
ax.axvline(0, color="grey", ls="--", lw=1)
ax.set_yticks(y, ["Unadjusted", "Composition-adjusted"])
ax.set_xlabel("Repair–Cilia Imbalance, SD per histology grade (95% CI)")
ax.set_title("GSE320381: exploratory composition adjustment")
fig.tight_layout(); fig.savefig(OUT / "figures/composition_adjusted_severity.png", dpi=260); plt.close(fig)

atten = adj["effect_size_SD_per_grade"]/unadj["effect_size_SD_per_grade"]*100
with open(OUT / "analysis_summary.md", "w") as fh:
    fh.write("# GSE320381 composition-adjusted analysis\n\n")
    fh.write(f"The GSE134174 never-smoker epithelial reference used {len(genes)} shared predefined markers. ")
    fh.write(f"The unadjusted severity coefficient was {unadj['effect_size_SD_per_grade']:.3f} SD/grade ")
    fh.write(f"(95% CI {unadj['CI95_lower']:.3f} to {unadj['CI95_upper']:.3f}, P={unadj['P']:.3g}); ")
    fh.write(f"after relative composition adjustment it was {adj['effect_size_SD_per_grade']:.3f} ")
    fh.write(f"(95% CI {adj['CI95_lower']:.3f} to {adj['CI95_upper']:.3f}, P={adj['P']:.3g}), ")
    fh.write(f"or {atten:.1f}% of the unadjusted coefficient.\n\n")
    fh.write("This is an exploratory cross-platform marker projection, not absolute deconvolution. ")
    fh.write("Persistence of a coefficient after adjustment is consistent with, but does not prove, a within-cell-state component.\n")
print(models.to_string(index=False))
