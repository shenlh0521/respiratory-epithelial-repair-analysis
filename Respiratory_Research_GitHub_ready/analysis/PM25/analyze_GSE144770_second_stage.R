#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
})

set.seed(20260811)
project <- Sys.getenv("PROJECT_ROOT", unset = normalizePath(".", mustWork = FALSE))
root <- file.path(project, "03_PM25/SECOND_VALIDATION/GSE144770")
data_dir <- file.path(root, "data")
out <- file.path(root, "results")
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)

gene_sets <- list(
  DNA_DAMAGE_RESPONSE = c("H2AX", "TP53BP1", "ATM", "ATR", "CHEK1", "CHEK2"),
  OXIDATIVE_STRESS = c("NFE2L2", "HMOX1", "NQO1", "GCLC", "GCLM", "SOD2", "TXNRD1"),
  ABNORMAL_REPAIR = c("KRT8", "KRT17", "KRT19", "CLDN4", "SFN", "LGALS3", "KRT14"),
  KRT14_KRT17_REPAIR = c("KRT14", "KRT17"),
  CILIA_CONSENSUS = c("FOXJ1", "MCIDAS", "GMNC", "MYB", "TP73", "RFX2", "RFX3", "CCNO",
                      "CDC20B", "DEUP1", "TUBB4B", "DNAH5", "DNAH9", "DNAH11", "DNAI1",
                      "DNAI2", "CCDC39", "CCDC40", "SPEF2", "HYDIN", "PIFO", "CFAP43",
                      "CFAP44", "RSPH1", "RSPH4A", "RSPH9"),
  MULTICILIOGENESIS = c("FOXJ1", "MCIDAS", "GMNC", "MYB", "TP73", "RFX2", "RFX3", "CCNO",
                        "CDC20B", "DEUP1", "CEP78", "CETN2", "PLK4", "STIL"),
  MATURE_CILIATED = c("FOXJ1", "PIFO", "TPPP3", "CAPS", "RSPH1")
)
focus_genes <- c("FOXJ1", "MCIDAS", "GMNC", paste0("RFX", 1:7), "KRT14", "KRT17")

paired_effect <- function(x1, x0, label1, label0) {
  ok <- is.finite(x1) & is.finite(x0)
  d <- x1[ok] - x0[ok]
  n <- length(d)
  est <- mean(d); se <- sd(d) / sqrt(n); df <- n - 1
  stat <- est / se
  p <- if (is.finite(stat)) 2 * pt(abs(stat), df, lower.tail = FALSE) else NA_real_
  crit <- qt(.975, df)
  data.frame(group_1 = label1, group_0 = label0, n_donors = n, estimate = est,
             std_error = se, CI95_lower = est - crit * se, CI95_upper = est + crit * se,
             statistic = stat, df = df, p_value = p)
}

one_sample_slope <- function(values, dose_steps) {
  ok <- is.finite(values) & is.finite(dose_steps)
  unname(coef(lm(values[ok] ~ dose_steps[ok]))[2])
}

counts <- fread(cmd = sprintf("gzip -cd %s", shQuote(file.path(data_dir, "GSE144770_AirPol_raw_counts.txt.gz"))),
                check.names = FALSE)
setnames(counts, 1, "gene")
counts[, gene := toupper(gene)]
counts[, (setdiff(names(counts), "gene")) := lapply(.SD, function(x) suppressWarnings(as.numeric(x))),
       .SDcols = setdiff(names(counts), "gene")]
counts <- counts[complete.cases(counts)]  # source matrix has three malformed non-signature cells
if (anyDuplicated(counts$gene)) counts <- counts[, lapply(.SD, sum), by = gene]
metadata <- fread(file.path(data_dir, "GSE144770_sample_metadata.csv"))
stopifnot(setequal(names(counts)[-1], metadata$matrix_id))
metadata <- metadata[match(names(counts)[-1], matrix_id)]

count_mat <- as.matrix(counts[, -1]); storage.mode(count_mat) <- "integer"; rownames(count_mat) <- counts$gene
dge <- DGEList(counts = count_mat)
dge <- calcNormFactors(dge)
logcpm <- cpm(dge, log = TRUE, prior.count = 1)
zexpr <- t(scale(t(logcpm)))
zexpr[!is.finite(zexpr)] <- NA_real_

score <- function(genes) {
  use <- intersect(genes, rownames(zexpr))
  if (!length(use)) return(rep(NA_real_, ncol(zexpr)))
  colMeans(zexpr[use, , drop = FALSE], na.rm = TRUE)
}
scores <- copy(metadata)
for (nm in names(gene_sets)) scores[[nm]] <- score(gene_sets[[nm]])
for (gene in intersect(focus_genes, rownames(zexpr))) scores[[gene]] <- zexpr[gene, ]
scores[, ABNORMAL_REPAIR_Z := as.numeric(scale(ABNORMAL_REPAIR))]
scores[, CILIA_CONSENSUS_Z := as.numeric(scale(CILIA_CONSENSUS))]
scores[, Repair_Cilia_Imbalance := ABNORMAL_REPAIR_Z - CILIA_CONSENSUS_Z]
fwrite(scores, file.path(out, "GSE144770_sample_scores.csv"))

coverage <- rbindlist(lapply(names(gene_sets), function(nm) data.frame(
  signature = nm, n_defined = length(gene_sets[[nm]]),
  n_detected = length(intersect(gene_sets[[nm]], rownames(zexpr))),
  detected_genes = paste(intersect(gene_sets[[nm]], rownames(zexpr)), collapse = ";"),
  missing_genes = paste(setdiff(gene_sets[[nm]], rownames(zexpr)), collapse = ";")
)))
fwrite(coverage, file.path(out, "signature_coverage.csv"))

outcomes <- unique(c(names(gene_sets), "Repair_Cilia_Imbalance", focus_genes))
outcomes <- intersect(outcomes, names(scores))

# Primary independent validation: all 12 donors, moderate organic extract versus matched control.
primary <- rbindlist(lapply(outcomes, function(nm) {
  wide <- dcast(scores[treatment_code %in% c("OECtrl", "7.5OE")], donor ~ treatment_code, value.var = nm)
  ans <- paired_effect(wide[["7.5OE"]], wide[["OECtrl"]], "moderate_OE", "OE_control")
  ans$outcome <- nm; ans
}))
primary[, FDR := p.adjust(p_value, "BH")]
primary[, comparison := "moderate organic PM2.5 extract vs matched vehicle control"]
primary[, analysis_method := "paired donor-level t test on prespecified mean gene-z scores"]
fwrite(primary, file.path(out, "moderate_vs_control_effects.csv"))

# Healthy-donor four-level dose response; estimate one slope per donor, then test donor slopes.
dose_data <- scores[asthma_status == "Healthy" & treatment_code %in% c("OECtrl", "0.75OE", "7.5OE", "75OE")]
complete_dose_donors <- dose_data[, .(n_levels = uniqueN(treatment_code)), by = donor][n_levels == 4, donor]
dose_data <- dose_data[donor %in% complete_dose_donors]
dose_effects <- rbindlist(lapply(outcomes, function(nm) {
  donor_slopes <- dose_data[, .(slope = one_sample_slope(get(nm), dose_step)), by = donor]$slope
  n <- length(donor_slopes); est <- mean(donor_slopes); se <- sd(donor_slopes) / sqrt(n); df <- n - 1
  stat <- est / se; p <- 2 * pt(abs(stat), df, lower.tail = FALSE); crit <- qt(.975, df)
  data.frame(outcome = nm, n_donors = n, estimate_per_dose_step = est, std_error = se,
             CI95_lower = est - crit * se, CI95_upper = est + crit * se,
             statistic = stat, df = df, p_value = p)
}))
dose_effects[, FDR := p.adjust(p_value, "BH")]
dose_effects[, comparison := "ordinal OE dose response among 5 healthy donors (control/low/moderate/high)"]
dose_effects[, analysis_method := "per-donor linear dose slope followed by one-sample t test"]
fwrite(dose_effects, file.path(out, "dose_response_effects.csv"))

# High-dose endpoint in the same five healthy donors.
high <- rbindlist(lapply(outcomes, function(nm) {
  wide <- dcast(dose_data[treatment_code %in% c("OECtrl", "75OE")], donor ~ treatment_code, value.var = nm)
  ans <- paired_effect(wide[["75OE"]], wide[["OECtrl"]], "high_OE", "OE_control")
  ans$outcome <- nm; ans
}))
high[, FDR := p.adjust(p_value, "BH")]
high[, comparison := "high organic PM2.5 extract vs matched vehicle control"]
high[, analysis_method := "paired donor-level t test"]
fwrite(high, file.path(out, "high_vs_control_effects.csv"))

primary_focus <- primary[outcome %in% c(names(gene_sets), "Repair_Cilia_Imbalance")]
primary_focus[, label := factor(outcome, levels = rev(outcome))]
p <- ggplot(primary_focus, aes(estimate, label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = CI95_lower, xmax = CI95_upper), height = .18) +
  geom_point(aes(fill = FDR < .05), shape = 21, size = 3) +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "white"), guide = "none") +
  labs(title = "GSE144770: paired moderate PM2.5 organic extract effect",
       subtitle = "12 primary human mucociliary airway epithelial donors",
       x = "Paired effect on mean gene-wise z score (95% CI)", y = NULL) +
  theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank())
ggsave(file.path(out, "figures", "GSE144770_PM25_validation.png"), p, width = 8.5, height = 5.5, dpi = 240)

dose_focus <- dose_effects[outcome %in% c(names(gene_sets), "Repair_Cilia_Imbalance")]
dose_focus[, label := factor(outcome, levels = rev(outcome))]
p_dose <- ggplot(dose_focus, aes(estimate_per_dose_step, label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = CI95_lower, xmax = CI95_upper), height = .18) +
  geom_point(aes(fill = FDR < .05), shape = 21, size = 3) +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "white"), guide = "none") +
  labs(title = "GSE144770: PM2.5 organic-extract dose response",
       subtitle = "Five healthy donors with matched control, low, moderate and high doses",
       x = "Per-donor effect per ordinal dose step (95% CI)", y = NULL) +
  theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank())
ggsave(file.path(out, "figures", "GSE144770_PM25_dose_response.png"), p_dose, width = 8.5, height = 5.5, dpi = 240)

key <- primary[outcome %in% c("DNA_DAMAGE_RESPONSE", "OXIDATIVE_STRESS", "ABNORMAL_REPAIR",
                               "KRT14_KRT17_REPAIR", "CILIA_CONSENSUS", "MULTICILIOGENESIS",
                               "MATURE_CILIATED", "Repair_Cilia_Imbalance")]
lines <- c("# Second-stage PM2.5 validation: GSE144770", "",
  "Dataset: 12 donor-matched primary human nasal mucociliary ALI cultures, with organic PM2.5 extract dose series and processed RNA-seq counts.",
  "Primary comparison is moderate organic extract versus matched vehicle across all 12 donors. Signatures are unchanged from the existing project.", "",
  "## Primary paired effects", "")
for (i in seq_len(nrow(key))) {
  z <- key[i]
  lines <- c(lines, sprintf("- %s: effect %.3f (95%% CI %.3f to %.3f), P=%s, FDR=%s.",
    z$outcome, z$estimate, z$CI95_lower, z$CI95_upper,
    format(z$p_value, digits = 3), format(z$FDR, digits = 3)))
}
lines <- c(lines, "", "Dose-response results are saved separately and use five healthy donors with all four OE levels.",
  sprintf("Dose response: OXIDATIVE_STRESS %.3f per dose step (FDR %s); CILIA_CONSENSUS %.3f (FDR %s); MULTICILIOGENESIS %.3f (FDR %s); MATURE_CILIATED %.3f (FDR %s); imbalance %.3f (FDR %s).",
    dose_effects[outcome == "OXIDATIVE_STRESS"]$estimate_per_dose_step, format(dose_effects[outcome == "OXIDATIVE_STRESS"]$FDR, digits = 3),
    dose_effects[outcome == "CILIA_CONSENSUS"]$estimate_per_dose_step, format(dose_effects[outcome == "CILIA_CONSENSUS"]$FDR, digits = 3),
    dose_effects[outcome == "MULTICILIOGENESIS"]$estimate_per_dose_step, format(dose_effects[outcome == "MULTICILIOGENESIS"]$FDR, digits = 3),
    dose_effects[outcome == "MATURE_CILIATED"]$estimate_per_dose_step, format(dose_effects[outcome == "MATURE_CILIATED"]$FDR, digits = 3),
    dose_effects[outcome == "Repair_Cilia_Imbalance"]$estimate_per_dose_step, format(dose_effects[outcome == "Repair_Cilia_Imbalance"]$FDR, digits = 3)),
  sprintf("High vs control: FOXJ1 %.3f (FDR %s), CILIA_CONSENSUS %.3f (FDR %s), MULTICILIOGENESIS %.3f (FDR %s), MATURE_CILIATED %.3f (FDR %s).",
    high[outcome == "FOXJ1"]$estimate, format(high[outcome == "FOXJ1"]$FDR, digits = 3),
    high[outcome == "CILIA_CONSENSUS"]$estimate, format(high[outcome == "CILIA_CONSENSUS"]$FDR, digits = 3),
    high[outcome == "MULTICILIOGENESIS"]$estimate, format(high[outcome == "MULTICILIOGENESIS"]$FDR, digits = 3),
    high[outcome == "MATURE_CILIATED"]$estimate, format(high[outcome == "MATURE_CILIATED"]$FDR, digits = 3)),
  "The system is differentiated nasal rather than bronchial epithelium, and the exposure is an extract rather than intact airborne PM2.5; these are explicit limitations.")
writeLines(lines, file.path(project, "03_PM25/SECOND_VALIDATION/second_stage_PM25_summary.md"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
