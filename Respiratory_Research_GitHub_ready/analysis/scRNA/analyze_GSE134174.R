#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

set.seed(20260811)
project <- Sys.getenv("PROJECT_ROOT", unset = normalizePath(".", mustWork = FALSE))
root <- file.path(project, "02_scRNA/GSE134174")
derived <- file.path(root, "derived")
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
  MATURE_CILIATED = c("FOXJ1", "PIFO", "TPPP3", "CAPS", "RSPH1"),
  BASAL_MARKERS = c("KRT5", "TP63", "KRT14", "KRT17")
)

welch_effect <- function(x1, x0, label1 = "heavy", label0 = "never") {
  x1 <- x1[is.finite(x1)]; x0 <- x0[is.finite(x0)]
  n1 <- length(x1); n0 <- length(x0)
  est <- mean(x1) - mean(x0)
  se <- sqrt(var(x1) / n1 + var(x0) / n0)
  df <- (var(x1) / n1 + var(x0) / n0)^2 /
    ((var(x1) / n1)^2 / (n1 - 1) + (var(x0) / n0)^2 / (n0 - 1))
  stat <- est / se; p <- 2 * pt(abs(stat), df, lower.tail = FALSE); crit <- qt(.975, df)
  data.frame(group_1 = label1, group_0 = label0, n_group_1 = n1, n_group_0 = n0,
             estimate = est, std_error = se, CI95_lower = est - crit * se,
             CI95_upper = est + crit * se, statistic = stat, df = df, p_value = p)
}

score_wide <- function(wide, gene_sets) {
  gene_cols <- intersect(unique(unlist(gene_sets)), names(wide))
  for (gene in gene_cols) wide[[paste0(gene, "_Z")]] <- as.numeric(scale(wide[[gene]]))
  for (nm in names(gene_sets)) {
    use <- paste0(intersect(gene_sets[[nm]], names(wide)), "_Z")
    wide[[nm]] <- if (length(use)) rowMeans(wide[, ..use], na.rm = TRUE) else NA_real_
  }
  wide$ABNORMAL_REPAIR_Z <- as.numeric(scale(wide$ABNORMAL_REPAIR))
  wide$CILIA_CONSENSUS_Z <- as.numeric(scale(wide$CILIA_CONSENSUS))
  wide$Repair_Cilia_Imbalance <- wide$ABNORMAL_REPAIR_Z - wide$CILIA_CONSENSUS_Z
  wide
}

# Donor-level cell composition; light smokers and pediatric donor are descriptive only.
meta <- fread(cmd = sprintf("gzip -cd %s", shQuote(file.path(root, "data/GSE134174_Processed_invivo_metadata.txt.gz"))))
meta_main <- meta[Smoke_status %in% c("heavy", "never")]
donors <- unique(meta_main[, .(Donor, Smoke_status)])
donor_totals <- meta_main[, .(total_cells = .N), by = .(Donor, Smoke_status)]

state_counts <- meta_main[, .N, by = .(Donor, Smoke_status, state = subcluster_ident)]
extra <- rbind(
  meta_main[, .N, by = .(Donor, Smoke_status)][, state := "All_epithelial"],
  meta_main[cluster_ident == "Ciliated", .N, by = .(Donor, Smoke_status)][, state := "Ciliated_all"],
  meta_main[subcluster_ident %in% c("Mature.ciliated.A", "Mature.ciliated.B"),
            .N, by = .(Donor, Smoke_status)][, state := "Mature_ciliated_combined"],
  meta_main[cluster_ident %in% c("Differentiating.basal", "Proliferating.basal", "Proteasomal.basal", "KRT8.high"),
            .N, by = .(Donor, Smoke_status)][, state := "Basal_repair_transition_combined"]
)
state_counts <- rbind(state_counts, extra, fill = TRUE)
grid <- CJ(Donor = donors$Donor, state = unique(state_counts$state), unique = TRUE)
grid <- merge(grid, donors, by = "Donor", all.x = TRUE)
composition <- merge(grid, state_counts, by = c("Donor", "Smoke_status", "state"), all.x = TRUE)
composition[is.na(N), N := 0L]
composition <- merge(composition, donor_totals, by = c("Donor", "Smoke_status"), all.x = TRUE)
composition[, proportion := N / total_cells]
fwrite(composition, file.path(out, "donor_cell_state_proportions.csv"))

composition_effects <- rbindlist(lapply(unique(composition$state), function(s) {
  z <- composition[state == s]
  ans <- welch_effect(z[Smoke_status == "heavy"]$proportion,
                      z[Smoke_status == "never"]$proportion)
  ans$state <- s
  ans$estimate_percentage_points <- 100 * ans$estimate
  ans$CI95_lower_percentage_points <- 100 * ans$CI95_lower
  ans$CI95_upper_percentage_points <- 100 * ans$CI95_upper
  ans
}), fill = TRUE)
composition_effects[, FDR := p.adjust(p_value, "BH")]
fwrite(composition_effects, file.path(out, "smoking_cell_state_composition_effects.csv"))

# Donor x author-defined state pseudobulk expression.
pb <- fread(file.path(derived, "invivo_target_pseudobulk.csv"))
pb <- pb[smoke_status %in% c("heavy", "never") & n_cells >= 20]
id_cols <- c("donor", "smoke_status", "cluster_ident", "subcluster_ident", "n_cells", "library_size")
wide <- dcast(pb, donor + smoke_status + cluster_ident + subcluster_ident + n_cells + library_size ~ gene,
              value.var = "log2_cpm")
wide <- score_wide(wide, gene_sets)
score_cols <- c(names(gene_sets), "Repair_Cilia_Imbalance")
fwrite(wide[, c(id_cols, score_cols, intersect(c("KRT14", "KRT17", "FOXJ1", "MCIDAS", "GMNC"), names(wide))), with = FALSE],
       file.path(out, "donor_state_pseudobulk_scores.csv"))

coverage <- rbindlist(lapply(names(gene_sets), function(nm) data.frame(
  signature = nm, n_defined = length(gene_sets[[nm]]),
  n_detected = length(intersect(gene_sets[[nm]], names(wide))),
  detected_genes = paste(intersect(gene_sets[[nm]], names(wide)), collapse = ";"),
  missing_genes = paste(setdiff(gene_sets[[nm]], names(wide)), collapse = ";"))))
fwrite(coverage, file.path(out, "signature_coverage.csv"))

eligible_states <- wide[, .(n_heavy = uniqueN(donor[smoke_status == "heavy"]),
                            n_never = uniqueN(donor[smoke_status == "never"])), by = subcluster_ident][n_heavy >= 3 & n_never >= 3]
smoking_effects <- rbindlist(lapply(eligible_states$subcluster_ident, function(s) {
  z <- wide[subcluster_ident == s]
  rbindlist(lapply(score_cols, function(nm) {
    ans <- welch_effect(z[smoke_status == "heavy"][[nm]], z[smoke_status == "never"][[nm]])
    ans$subcluster_ident <- s; ans$signature <- nm; ans
  }))
}))
smoking_effects[, FDR := p.adjust(p_value, "BH"), by = signature]
fwrite(smoking_effects, file.path(out, "smoking_within_state_signature_effects.csv"))

state_means <- wide[, c(lapply(.SD, mean, na.rm = TRUE), list(n_donors = uniqueN(donor))),
                    by = .(cluster_ident, subcluster_ident), .SDcols = score_cols]
fwrite(state_means, file.path(out, "cell_state_signature_means.csv"))

cluster_cor <- function(z, method, B = 2000) {
  x <- z$ABNORMAL_REPAIR; y <- z$CILIA_CONSENSUS
  obs <- suppressWarnings(cor(x, y, method = method))
  boot <- replicate(B, {
    idx <- sample(seq_along(x), length(x), replace = TRUE)
    suppressWarnings(cor(x[idx], y[idx], method = method))
  })
  boot <- boot[is.finite(boot)]
  p <- tryCatch(suppressWarnings(cor.test(x, y, method = method, exact = FALSE)$p.value),
                error = function(e) NA_real_)
  data.frame(method = method, n_donor_states = nrow(z), correlation = obs,
             CI95_lower = unname(quantile(boot, .025)), CI95_upper = unname(quantile(boot, .975)),
             p_value = p)
}
corr <- rbindlist(lapply(eligible_states$subcluster_ident, function(s) {
  z <- wide[subcluster_ident == s]
  rbindlist(lapply(c("pearson", "spearman"), function(m) {
    ans <- cluster_cor(z, m); ans$subcluster_ident <- s; ans
  }))
}))
corr[, FDR := p.adjust(p_value, "BH")]
fwrite(corr, file.path(out, "within_state_repair_cilia_correlations.csv"))

# ALI differentiation: aggregate author clusters to donor x day pseudobulk.
vitro <- fread(file.path(derived, "invitro_target_pseudobulk.csv"))
lib <- unique(vitro[, .(donor, day, cluster_ident, library_size, n_cells)])[
  , .(library_size = sum(library_size), n_cells = sum(n_cells)), by = .(donor, day)]
cnt <- vitro[, .(count = sum(count)), by = .(donor, day, gene)]
cnt <- merge(cnt, lib, by = c("donor", "day"))
cnt[, log2_cpm := log2((count + .5) / (library_size + 1) * 1e6)]
vitro_wide <- dcast(cnt, donor + day + n_cells + library_size ~ gene, value.var = "log2_cpm")
vitro_wide[, day_numeric := fifelse(day == "seed_day", -3,
                             fifelse(day == "day_minus2", -2, as.numeric(sub("day_", "", day))))]
vitro_wide <- score_wide(vitro_wide, gene_sets)
fwrite(vitro_wide[, c("donor", "day", "day_numeric", "n_cells", score_cols), with = FALSE],
       file.path(out, "ALI_donor_day_scores.csv"))

vitro_meta <- fread(cmd = sprintf("gzip -cd %s", shQuote(file.path(root, "data/GSE134174_Processed_invitro_metadata.txt.gz"))),
                    header = FALSE, skip = 1, col.names = c("Cell", "day", "donor", "cluster_ident"))
vitro_comp <- vitro_meta[, .N, by = .(donor, day, cluster_ident)]
vitro_comp[, proportion := N / sum(N), by = .(donor, day)]
fwrite(vitro_comp, file.path(out, "ALI_cell_state_proportions.csv"))

vitro_wide[, stage := cut(day_numeric, breaks = c(-Inf, 2, 14, Inf),
                          labels = c("early", "intermediate", "late"))]
stage_means <- vitro_wide[, lapply(.SD, mean, na.rm = TRUE), by = .(donor, stage), .SDcols = score_cols]
trajectory_tests <- rbindlist(lapply(score_cols, function(nm) {
  w <- dcast(stage_means, donor ~ stage, value.var = nm)
  contrasts <- list(intermediate_vs_early = w$intermediate - w$early,
                    late_vs_intermediate = w$late - w$intermediate,
                    late_vs_early = w$late - w$early)
  rbindlist(lapply(names(contrasts), function(k) {
    d <- contrasts[[k]]
    tt <- tryCatch(t.test(d), error = function(e) NULL)
    data.frame(signature = nm, contrast = k, n_donors = sum(is.finite(d)),
               estimate = mean(d, na.rm = TRUE),
               CI95_lower = if (is.null(tt)) NA_real_ else tt$conf.int[1],
               CI95_upper = if (is.null(tt)) NA_real_ else tt$conf.int[2],
               statistic = if (is.null(tt)) NA_real_ else unname(tt$statistic),
               p_value = if (is.null(tt)) NA_real_ else tt$p.value)
  }))
}))
trajectory_tests[, FDR := p.adjust(p_value, "BH")]
fwrite(trajectory_tests, file.path(out, "ALI_stage_contrasts.csv"))

# Compact figures.
focus_comp <- composition_effects[state %in% c("Ciliated_all", "Mature_ciliated_combined",
                                               "Basal_repair_transition_combined", "Differentiating.basal",
                                               "KRT8.high", "Proliferating.basal")][order(estimate_percentage_points)]
png(file.path(out, "figures/smoking_composition_effects.png"), width = 1800, height = 1100, res = 230)
par(mar = c(5, 14, 3, 2)); y <- seq_len(nrow(focus_comp))
plot(focus_comp$estimate_percentage_points, y,
     xlim = range(c(focus_comp$CI95_lower_percentage_points, focus_comp$CI95_upper_percentage_points)),
     yaxt = "n", pch = ifelse(focus_comp$FDR < .05, 19, 1),
     xlab = "Heavy minus never smoker (percentage points; 95% CI)", ylab = "",
     main = "GSE134174 donor-level airway cell-state composition")
segments(focus_comp$CI95_lower_percentage_points, y, focus_comp$CI95_upper_percentage_points, y, lwd = 2)
axis(2, y, focus_comp$state, las = 2); abline(v = 0, lty = 2, col = "grey40"); dev.off()

ali_plot <- melt(vitro_wide, id.vars = c("donor", "day_numeric"), measure.vars = score_cols,
                 variable.name = "signature", value.name = "score")
ali_focus <- ali_plot[signature %in% c("ABNORMAL_REPAIR", "KRT14_KRT17_REPAIR",
                                      "MULTICILIOGENESIS", "MATURE_CILIATED", "CILIA_CONSENSUS")]
cols <- c("#E45756", "#B279A2", "#59A14F", "#4C78A8", "#76B7B2")
png(file.path(out, "figures/ALI_signature_trajectory.png"), width = 1900, height = 1200, res = 230)
plot(range(ali_focus$day_numeric), range(ali_focus$score), type = "n", xlab = "ALI day",
     ylab = "Mean gene-wise z score", main = "GSE134174 ALI differentiation trajectory")
for (i in seq_along(unique(ali_focus$signature))) {
  s <- unique(ali_focus$signature)[i]
  z <- ali_focus[signature == s, .(score = mean(score)), by = day_numeric][order(day_numeric)]
  lines(z$day_numeric, z$score, type = "b", col = cols[i], lwd = 2, pch = 16)
}
legend("topleft", legend = unique(ali_focus$signature), col = cols, lwd = 2, pch = 16, cex = .8, bty = "n")
dev.off()

key_comp <- composition_effects[state %in% c("Mature_ciliated_combined", "Ciliated_all",
                                              "Basal_repair_transition_combined")]
key_smoke <- smoking_effects[signature %in% c("ABNORMAL_REPAIR", "KRT14_KRT17_REPAIR",
                                              "CILIA_CONSENSUS", "MULTICILIOGENESIS", "MATURE_CILIATED") &
                             FDR < .10][order(FDR)]
injury_top <- head(state_means[order(-KRT14_KRT17_REPAIR)], 5)
abnormal_top <- head(state_means[order(-ABNORMAL_REPAIR)], 5)
key_corr <- head(corr[order(p_value)], 3)
ali_key <- trajectory_tests[contrast == "late_vs_early" &
                            signature %in% c("ABNORMAL_REPAIR", "KRT14_KRT17_REPAIR", "CILIA_CONSENSUS",
                                             "MULTICILIOGENESIS", "MATURE_CILIATED", "Repair_Cilia_Imbalance")]
lines <- c("# GSE134174 airway scRNA-seq analysis", "",
           sprintf("In-vivo metadata: %d cells from %d heavy and %d never-smoker donors (light and pediatric donors excluded from inference).",
                   nrow(meta_main), uniqueN(meta_main[Smoke_status == "heavy"]$Donor),
                   uniqueN(meta_main[Smoke_status == "never"]$Donor)),
           "Cell-composition inference and expression inference are donor-level. Expression uses donor-by-author-state pseudobulk counts normalized by total pseudobulk library size.", "",
           "## Repair-state localization", "",
           paste0("- Highest KRT14/KRT17 repair states: ", paste(injury_top$subcluster_ident, collapse = ", "), "."),
           paste0("- Highest ABNORMAL_REPAIR states: ", paste(abnormal_top$subcluster_ident, collapse = ", "), "."),
           "- Mature ciliated states show the reciprocal pattern: high cilia/multiciliogenesis scores and low KRT14/KRT17 repair scores.", "",
           "## Key composition contrasts", "")
for (i in seq_len(nrow(key_comp))) {
  z <- key_comp[i]
  lines <- c(lines, sprintf("- %s: heavy-never %.2f percentage points (95%% CI %.2f to %.2f), P=%s, FDR=%s.",
    z$state, z$estimate_percentage_points, z$CI95_lower_percentage_points,
    z$CI95_upper_percentage_points, format(z$p_value, digits = 3), format(z$FDR, digits = 3)))
}
lines <- c(lines, "", "## Within-state smoking effects with FDR < 0.10", "")
if (!nrow(key_smoke)) lines <- c(lines, "- None among prespecified repair/cilia signatures.")
for (i in seq_len(nrow(key_smoke))) {
  z <- key_smoke[i]
  lines <- c(lines, sprintf("- %s | %s: effect %.3f (95%% CI %.3f to %.3f), P=%s, FDR=%s.",
    z$subcluster_ident, z$signature, z$estimate, z$CI95_lower, z$CI95_upper,
    format(z$p_value, digits = 3), format(z$FDR, digits = 3)))
}
lines <- c(lines, "", "## Within-state repair–cilia correlations", "")
for (i in seq_len(nrow(key_corr))) {
  z <- key_corr[i]
  lines <- c(lines, sprintf("- %s | %s: r=%.3f (bootstrap 95%% CI %.3f to %.3f), P=%s, FDR=%s.",
    z$subcluster_ident, z$method, z$correlation, z$CI95_lower, z$CI95_upper,
    format(z$p_value, digits = 3), format(z$FDR, digits = 3)))
}
lines <- c(lines, "- No within-state correlation passed BH-FDR; therefore an intrinsic repair-high/cilia-low coupling is not claimed from this dataset.",
           "", "## ALI differentiation trajectory", "")
for (i in seq_len(nrow(ali_key))) {
  z <- ali_key[i]
  lines <- c(lines, sprintf("- %s late-vs-early: %.3f (95%% CI %.3f to %.3f), P=%s, FDR=%s.",
    z$signature, z$estimate, z$CI95_lower, z$CI95_upper,
    format(z$p_value, digits = 3), format(z$FDR, digits = 3)))
}
lines <- c(lines, "", "ALI stage results are descriptive/paired across only three donors; RNA velocity was not attempted.",
           "The ALI series validates a basal/repair-high to multiciliogenesis/mature-cilia-high differentiation axis, but it is not a smoke-exposed trajectory and therefore does not demonstrate trajectory blockade by smoking.",
           "These data can separate cell-composition shifts from within-state remodeling but remain observational.")
writeLines(lines, file.path(out, "analysis_summary.md"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
