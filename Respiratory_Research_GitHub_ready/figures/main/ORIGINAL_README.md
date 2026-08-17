# Figure 1-8 source archive

`source_data/Figure_1-8_Source_Data.xlsx` is the Excel-grounded source for Figures 1-8. It contains one named sheet per figure plus a source manifest. Critical statistics are read from frozen project outputs and are not typed into plotting code.

Run `python rebuild_all_figures.py` to extract the workbook sheets and recreate PDF, SVG, and 600-dpi PNG outputs. Figure 6 uses presentation-only equal-cell pseudotime bins generated under the prespecified minimum of 150 cells and at least two donors; the Slingshot pseudotime and frozen donor-aware statistics are unchanged.

The `outputs` directory records the rebuilt outputs used in the live figure directory. The `code` directory contains the workbook builder and plotting scripts.
