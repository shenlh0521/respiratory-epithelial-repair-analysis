#!/usr/bin/env Rscript
project <- Sys.getenv("PROJECT_ROOT", unset = normalizePath(".", mustWork = FALSE))
legacy_lib <- file.path(project, "files-mentioned-by-the-user-users/outputs/lung_bulla_cilia_cancer/work/Rlib")
.libPaths(c(legacy_lib, .libPaths()))
suppressPackageStartupMessages({library(data.table); library(AnnotationDbi); library(illuminaHumanv4.db)})

legacy <- file.path(project, "files-mentioned-by-the-user-users/outputs/lung_bulla_cilia_cancer")
out <- file.path(project, "05_stretch/GSE59128_unified")
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)

d <- fread(file.path(legacy, "results/GSE59128/metadata_clean.csv"))
d[, cell_line_id := sub("^cell line id: ", "", characteristic_3)]

# Recover FOXJ1 from the already-downloaded processed matrix; no data are downloaded.
lines <- readLines(gzfile(file.path(legacy, "data/raw/GSE59128_series_matrix.txt.gz")), warn = FALSE)
b <- grep("^!series_matrix_table_begin", lines); e <- grep("^!series_matrix_table_end", lines)
tab <- fread(text = paste(lines[(b + 1):(e - 1)], collapse = "\n"), check.names = FALSE, data.table = FALSE)
probe <- tab[[1]]; mat <- as.matrix(tab[, -1, drop = FALSE]); storage.mode(mat) <- "double"; rownames(mat) <- probe
ann <- AnnotationDbi::select(illuminaHumanv4.db::illuminaHumanv4.db, keys = probe,
                             columns = "SYMBOL", keytype = "PROBEID")
fox_probes <- ann$PROBEID[toupper(ann$SYMBOL) == "FOXJ1"]
fox_probes <- intersect(fox_probes, rownames(mat))
if (!length(fox_probes)) stop("FOXJ1 not mapped on GPL10558")
fox_probe <- fox_probes[which.max(rowMeans(mat[fox_probes, , drop = FALSE], na.rm = TRUE))]
fox <- as.numeric(scale(mat[fox_probe, ])); names(fox) <- colnames(mat)
d[, FOXJ1 := fox[sample]]
rm(mat, tab, lines)

d[, ABNORMAL_REPAIR_Z := as.numeric(scale(ABNORMAL_REPAIR))]
d[, CILIA_CONSENSUS_Z := as.numeric(scale(CILIA_CONSENSUS))]
d[, Repair_Cilia_Imbalance := ABNORMAL_REPAIR_Z - CILIA_CONSENSUS_Z]

signature_map <- c(
  DNA_DAMAGE_RESPONSE = "DNA_DAMAGE_SENTINELS",
  OXIDATIVE_STRESS = "OXIDATIVE_STRESS",
  ABNORMAL_REPAIR = "ABNORMAL_REPAIR",
  KRT14_KRT17_REPAIR = "INJURY_BASAL",
  CILIA_CONSENSUS = "CILIA_CONSENSUS",
  FOXJ1 = "FOXJ1",
  Repair_Cilia_Imbalance = "Repair_Cilia_Imbalance",
  YAP_TAZ_MECHANOTRANSDUCTION = "YAP_TAZ_MECHANOTRANSDUCTION",
  PIEZO_MECHANICAL = "PIEZO_MECHANICAL",
  NOTCH = "NOTCH"
)

results <- list(); paired_rows <- list()
for (tt in c(8, 24)) {
  for (nm in names(signature_map)) {
    col <- signature_map[[nm]]
    a <- d[cell_state == "AEC" & treatment == "stretch" & time_h == tt,
           .(cell_line_id, stretch = get(col))]
    b0 <- d[cell_state == "AEC" & treatment == "sham" & time_h == tt,
            .(cell_line_id, sham = get(col))]
    p <- merge(a, b0, by = "cell_line_id")
    p[, delta := stretch - sham]
    ttst <- t.test(p$delta)
    results[[length(results) + 1]] <- data.frame(
      time_h = tt, signature = nm, n_pairs = nrow(p),
      mean_delta_stretch_minus_sham = mean(p$delta),
      CI95_lower = ttst$conf.int[1], CI95_upper = ttst$conf.int[2],
      paired_effect_dz = mean(p$delta) / sd(p$delta),
      statistic = unname(ttst$statistic), p_value = ttst$p.value,
      analysis_method = "matched cell-line paired t-test", stringsAsFactors = FALSE)
    p[, `:=`(time_h = tt, signature = nm)]
    paired_rows[[length(paired_rows) + 1]] <- p
  }
}
results <- rbindlist(results)
results[, FDR := p.adjust(p_value, "BH"), by = time_h]
fwrite(results, file.path(out, "paired_signature_effects.csv"))
fwrite(rbindlist(paired_rows), file.path(out, "paired_signature_differences.csv"))
fwrite(d[, c("sample", "cell_line_id", "cell_state", "treatment", "time_h",
             unname(signature_map)), with = FALSE], file.path(out, "unified_scores.csv"))

png(file.path(out, "figures/paired_stretch_effects.png"), width = 2000, height = 1400, res = 240)
par(mar = c(5, 15, 3, 2)); z <- results[order(time_h, signature)]; y <- seq_len(nrow(z))
plot(z$mean_delta_stretch_minus_sham, y,
     xlim = range(c(z$CI95_lower, z$CI95_upper)), yaxt = "n",
     pch = ifelse(z$FDR < .05, 19, 1), xlab = "Stretch minus matched sham (95% CI)", ylab = "",
     main = "GSE59128 unified paired stretch analysis")
segments(z$CI95_lower, y, z$CI95_upper, y, lwd = 2)
axis(2, y, paste0(z$time_h, "h | ", z$signature), las = 2, cex.axis = .72)
abline(v = 0, lty = 2, col = "grey40"); dev.off()

focus <- results[signature %in% c("DNA_DAMAGE_RESPONSE", "ABNORMAL_REPAIR", "KRT14_KRT17_REPAIR",
                                  "CILIA_CONSENSUS", "FOXJ1", "Repair_Cilia_Imbalance",
                                  "YAP_TAZ_MECHANOTRANSDUCTION", "PIEZO_MECHANICAL", "NOTCH")]
txt <- c("# GSE59128 unified mechanical-stretch analysis", "",
         "The processed matrix and canonical project signatures were reused; no data were downloaded.",
         "Cell-line IDs in GEO characteristics permit four matched stretch-versus-sham pairs at each time point. This supersedes the earlier unpaired 4-versus-8 screen for inference while leaving prior results untouched.", "",
         "## Effects", "")
for (i in seq_len(nrow(focus))) {
  z <- focus[i]
  txt <- c(txt, sprintf("- %dh | %s: mean paired delta %.3f (95%% CI %.3f to %.3f), dz=%.3f, P=%s, FDR=%s.",
    z$time_h, z$signature, z$mean_delta_stretch_minus_sham, z$CI95_lower, z$CI95_upper,
    z$paired_effect_dz, format(z$p_value, digits = 3), format(z$FDR, digits = 3)))
}
txt <- c(txt, "", "Interpretation emphasizes the most stable downstream stress response. YAP/PIEZO/NOTCH are evaluated without assuming sustained activation of any single upstream pathway.")
writeLines(txt, file.path(out, "analysis_summary.md"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
