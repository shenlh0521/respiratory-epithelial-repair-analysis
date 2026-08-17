# GitHub cleanup report

This public-code package was prepared from the supplied frozen Respiratory Research submission archive.

## Changes made

- Copied the final analysis scripts into `analysis/`.
- Removed hard-coded author-machine `/Users/...` project roots and replaced them with the `PROJECT_ROOT` environment variable.
- Copied the Main Figures 1-8 rebuild code and frozen source workbook/CSV tables into `figures/main/`.
- Redirected regenerated main figures to `figures/main/outputs/` instead of the original submission-package figure directory.
- Copied Supplementary Figures S1-S25 rebuild code into `figures/supplementary/` and made its external frozen-project root configurable through `PROJECT_ROOT`.
- Redirected regenerated supplementary figures to `figures/supplementary/outputs/`.
- Omitted already-generated PDF/PNG/SVG figure outputs from the repository package.
- Omitted macOS metadata files and added `.gitignore` rules for generated outputs, environments, and large raw-data formats.
- Added `README.md`, `requirements.txt`, `R_PACKAGES.md`, metadata inventories, and a license note.

## Static checks performed

- No remaining `/Users/...`, `/home/...`, or `/mnt/...` absolute paths were found in the cleaned package.
- No obvious API-key/password/private-key patterns were found in the cleaned package.
- All Python files passed Python bytecode compilation (`compileall`).
- The archived main-figure workbook is readable and contains the expected figure source sheets.

## Known limitations

- Exact Python/R package versions were not recorded in the supplied submission archive, so versions were not fabricated or pinned.
- Large downloaded GEO files and many original frozen intermediate result tables are not included. Therefore the analysis scripts are transparent/archival but are not all fully runnable from a fresh clone alone.
- The main-figure rebuild archive is self-contained with respect to its frozen plotting source data.
- The supplementary-figure rebuild script still requires original frozen project outputs listed in `figures/supplementary/REQUIRED_EXTERNAL_INPUTS.md`.
- R syntax/runtime execution was not performed in the packaging environment because `Rscript` was not available there; the R modifications were limited to path/output configuration lines.
