# Required external frozen inputs

The supplementary-figure script is archival code and is **not fully self-contained** in this repository. It reads frozen result tables from the original analysis project. Set `PROJECT_ROOT` to the root of that project before running it.

The script references these project-relative files:

- `01_PML_extended/GSE109743/pml_extended_scores.csv`
- `01_PML_extended/GSE109743/pml_grade_trends.csv`
- `01_PML_extended/GSE109743/pml_subtype_effects.csv`
- `01_PML_extended/GSE109743/repair_cilia_correlations.csv`
- `01_PML_extended/GSE114489/results/persistent_vs_regressive_effects.csv`
- `01_PML_extended/GSE79210/results/baseline_persistent_vs_regressive.csv`
- `02_scRNA/GSE134174/results/cell_state_signature_means.csv`
- `02_scRNA/GSE134174/results/smoking_cell_state_composition_effects.csv`
- `02_scRNA/GSE134174/third_stage_results/ALI_full_timecourse_scores.csv`
- `02_scRNA/SECOND_VALIDATION/GSE233145/results/cell_state_signature_profiles.csv`
- `05_stretch/GSE59128_unified/paired_signature_effects.csv`
- `11_BMC_Pulmonary_Medicine_submission/10_source_archive/provenance/final_result_provenance.csv`
- `12_targeted_manuscript_upgrades/01_deconvolution_robustness/benchmark_method_status.csv`
- `12_targeted_manuscript_upgrades/01_deconvolution_robustness/deconvolution_method_comparison.csv`
- `12_targeted_manuscript_upgrades/02_ALI_trajectory/slingshot_cell_pseudotime.csv`
- `12_targeted_manuscript_upgrades/03_ALI_maturation_projection/ALI_maturation_signature.csv`
- `12_targeted_manuscript_upgrades/03_ALI_maturation_projection/GSE114489_ALI_maturation_scores.csv`
- `12_targeted_manuscript_upgrades/03_ALI_maturation_projection/PML_maturation_projection_results.csv`
- `12_targeted_manuscript_upgrades/03_ALI_maturation_projection/platform_coverage.csv`
- `12_targeted_manuscript_upgrades/04_DECAMP_validation/GSE300258_frozen_scores.csv`
