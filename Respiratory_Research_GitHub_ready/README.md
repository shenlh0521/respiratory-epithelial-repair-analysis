# Respiratory Research analysis code and figure-reproduction archive

This repository is a cleaned public-code package prepared from the frozen submission archive. It contains the final analysis scripts, the self-contained main-figure source data and rebuild code, supplementary-figure rebuild code, and compact provenance/inventory files.

## Repository structure

```text
analysis/                     Final R/Python analysis scripts
  PML/
  scRNA/
  PM25/
  smoke/
  stretch/
  targeted_upgrades/
figures/
  main/                       Main Figures 1-8 rebuild code + frozen source data
  supplementary/              Supplementary Figures S1-S25 rebuild code
metadata/                     Code/data inventories and analysis-freeze note
requirements.txt              Python dependencies (unversioned; original versions unavailable)
R_PACKAGES.md                 R/Bioconductor dependency list
```

## Data availability

The analyses use public GEO datasets listed in `metadata/Data_Availability_Inventory.csv`. Large raw matrices and downloaded public data are intentionally not included in this repository. The archived submission states that processed public data were used preferentially and that no participant-level restricted data or large raw sequencing matrices are included.

## Path configuration for analysis scripts

Hard-coded author-machine paths were removed. Analysis scripts now obtain the frozen analysis-project root from the environment variable `PROJECT_ROOT`. If it is not set, they default to the current working directory.

Example:

```bash
export PROJECT_ROOT=/path/to/frozen_analysis_project
Rscript analysis/PML/analyze_GSE109743_extended.R
python analysis/scRNA/prepare_GSE134174_pseudobulk.py
```

The scripts retain the original project-relative folder layout beneath `PROJECT_ROOT` so the analysis logic is unchanged. Because the public repository does not contain all large downloaded GEO files or all intermediate frozen outputs, not every analysis script is one-command runnable from a fresh clone alone.

## Rebuild Main Figures 1-8

The main-figure archive is self-contained with respect to the frozen plotting source tables. It uses `figures/main/source_data/Figure_1-8_Source_Data.xlsx`.

```bash
python figures/main/rebuild_all_figures.py
```

This extracts the workbook sheets and invokes the R plotting script. Outputs are written to `figures/main/outputs/`. Requirements include Python/pandas/openpyxl, R, `ggplot2`, `patchwork`, `openxlsx`, and a Cairo-capable graphics setup. The plotting code was authored using Arial; font substitution may be needed on systems where Arial is unavailable.

## Rebuild Supplementary Figures S1-S25

```bash
export PROJECT_ROOT=/path/to/frozen_analysis_project
python figures/supplementary/rebuild_supplementary_figures.py
```

The supplementary builder also uses frozen result files from the original analysis project. See `figures/supplementary/REQUIRED_EXTERNAL_INPUTS.md` for the referenced inputs. Two Figure 6 QC/source tables needed by the supplementary builder are bundled under `figures/main/source_data/`.

## Reproducibility notes

- The frozen submission described the analysis results as fixed; this packaging changes paths and repository organization only, not scientific calculations. See `metadata/ANALYSIS_FREEZE.md`.
- Python and R package versions were not present in the supplied archive, so exact version pins are not fabricated here. Add the original environment/session information if available.
- Generated figure files are omitted from version control by default; the frozen source tables and rebuild scripts are retained.
- Before public release, add an explicit software license and replace any manuscript placeholder repository URL with the final GitHub URL.

## Suggested citation / code availability wording

After the repository is public, the manuscript can state: “Analysis and figure-reproduction code are available at the project GitHub repository [URL]. Public datasets are identified by GEO accession in the repository data-availability inventory.” For the final accepted version, consider archiving a tagged release in Zenodo and citing its DOI.
