#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

set.seed(20260811)
project <- Sys.getenv("PROJECT_ROOT", unset = normalizePath(".", mustWork = FALSE))
legacy <- file.path(project, "files-mentioned-by-the-user-users/outputs/lung_bulla_cilia_cancer")
out <- file.path(project, "01_PML_extended/GSE109743")
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)

md <- fread(file.path(legacy, "results/GSE109743/metadata_clean.csv"), data.table = FALSE)
md$patient_id <- factor(md$patient_id)
md$flow_cell_id <- factor(md$flow_cell_id)
md$molecular_subtype <- factor(md$molecular_subtype,
                               levels = c("Normal", "Inflammatory", "Proliferative", "Secretory"))
bio <- as.logical(md$passed_qc_flag) & grepl("biopsy", md$source, ignore.case = TRUE)
d <- md[bio, , drop = FALSE]
d$patient_id <- droplevels(d$patient_id)
d$flow_cell_id <- droplevels(d$flow_cell_id)

# Extract only prespecified genes from the existing residual matrix.
target_map <- fread(file.path(legacy, "data/metadata/target_symbol_ensembl.csv"), data.table = FALSE)
targets <- c("FOXJ1", "KRT14", "KRT17")
target_map <- target_map[target_map$symbol %in% targets, , drop = FALSE]
expr_path <- file.path(legacy, "data/raw/pending/GSE109743_residuals.txt.gz")
expr <- fread(cmd = sprintf("gzip -cd %s", shQuote(expr_path)), check.names = FALSE,
              data.table = FALSE)
names(expr)[1] <- "ensembl_id"
expr$ensembl_id <- sub("\\.[0-9]+$", "", expr$ensembl_id)
expr <- expr[expr$ensembl_id %in% target_map$ensembl_id, , drop = FALSE]
for (gene in targets) {
  ens <- target_map$ensembl_id[target_map$symbol == gene][1]
  hit <- expr[expr$ensembl_id == ens, , drop = FALSE]
  if (!nrow(hit)) {
    d[[gene]] <- NA_real_
  } else {
    values <- as.numeric(hit[1, -1, drop = TRUE])
    names(values) <- names(hit)[-1]
    d[[gene]] <- values[d$sample_id]
  }
}
rm(expr)

zscore <- function(x) as.numeric(scale(x))
rank_normal <- function(x) {
  ans <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  n <- sum(ok)
  ans[ok] <- qnorm((rank(x[ok], ties.method = "average") - 0.5) / n)
  ans
}
d$ABNORMAL_REPAIR_Z <- zscore(d$ABNORMAL_REPAIR)
d$CILIA_CONSENSUS_Z <- zscore(d$CILIA_CONSENSUS)
d$Repair_Cilia_Imbalance <- d$ABNORMAL_REPAIR_Z - d$CILIA_CONSENSUS_Z
d$ABNORMAL_REPAIR_RANKN <- rank_normal(d$ABNORMAL_REPAIR)
d$CILIA_CONSENSUS_RANKN <- rank_normal(d$CILIA_CONSENSUS)
d$Repair_Cilia_Imbalance_rank <- d$ABNORMAL_REPAIR_RANKN - d$CILIA_CONSENSUS_RANKN

cluster_fit <- function(data, outcome, rhs, cluster = "patient_id") {
  vars <- unique(c(outcome, all.vars(as.formula(paste("~", rhs))), cluster))
  z <- data[complete.cases(data[, vars, drop = FALSE]), , drop = FALSE]
  y <- z[[outcome]]
  X0 <- model.matrix(as.formula(paste("~", rhs)), data = z)
  q <- qr(X0)
  keep <- sort(q$pivot[seq_len(q$rank)])
  X <- X0[, keep, drop = FALSE]
  beta <- solve(crossprod(X), crossprod(X, y))
  resid <- as.numeric(y - X %*% beta)
  bread <- solve(crossprod(X))
  groups <- split(seq_len(nrow(z)), droplevels(factor(z[[cluster]])), drop = TRUE)
  meat <- matrix(0, ncol(X), ncol(X))
  for (idx in groups) {
    u <- as.numeric(crossprod(X[idx, , drop = FALSE], resid[idx]))
    meat <- meat + tcrossprod(u)
  }
  G <- length(groups); N <- nrow(z); P <- ncol(X)
  correction <- if (G > 1 && N > P) (G / (G - 1)) * ((N - 1) / (N - P)) else 1
  vc <- correction * bread %*% meat %*% bread
  list(beta = beta, vcov = vc, names = colnames(X), n = N, clusters = G,
       df = max(G - 1, 1), data = z)
}

contrast_result <- function(fit, weights, label, outcome, scale) {
  cvec <- setNames(rep(0, length(fit$beta)), fit$names)
  common <- intersect(names(weights), fit$names)
  cvec[common] <- weights[common]
  est <- sum(cvec * fit$beta)
  se <- sqrt(drop(t(cvec) %*% fit$vcov %*% cvec))
  stat <- est / se
  p <- 2 * pt(abs(stat), df = fit$df, lower.tail = FALSE)
  crit <- qt(0.975, df = fit$df)
  data.frame(outcome = outcome, contrast = label, scale = scale,
             n_samples = fit$n, n_patients = fit$clusters,
             estimate = est, std_error = se,
             CI95_lower = est - crit * se, CI95_upper = est + crit * se,
             statistic = stat, p_value = p, stringsAsFactors = FALSE)
}

subtypes <- levels(d$molecular_subtype)
pair_weights <- function(a, b) {
  # Returns a minus b under treatment coding with Normal reference.
  w <- numeric()
  if (a != "Normal") w[paste0("molecular_subtype", a)] <- 1
  if (b != "Normal") w[paste0("molecular_subtype", b)] <- -1
  w
}

outcomes <- c("CILIA_CONSENSUS", "ABNORMAL_REPAIR", "FOXJ1", "KRT14", "KRT17",
              "Repair_Cilia_Imbalance", "Repair_Cilia_Imbalance_rank")
scale_map <- c(CILIA_CONSENSUS = "legacy mean-gene z score",
               ABNORMAL_REPAIR = "legacy mean-gene z score",
               FOXJ1 = "study residual expression", KRT14 = "study residual expression",
               KRT17 = "study residual expression",
               Repair_Cilia_Imbalance = "biopsy-standardized z difference",
               Repair_Cilia_Imbalance_rank = "inverse-normal rank difference")
rhs_subtype <- "molecular_subtype + median_tin_num + flow_cell_id + patient_id"
subtype_results <- list()
for (outcome in outcomes) {
  fit <- cluster_fit(d, outcome, rhs_subtype)
  for (i in seq_len(length(subtypes) - 1)) for (j in (i + 1):length(subtypes)) {
    a <- subtypes[j]; b <- subtypes[i]
    subtype_results[[length(subtype_results) + 1]] <- contrast_result(
      fit, pair_weights(a, b), paste0(a, "_vs_", b), outcome, scale_map[[outcome]])
  }
}
subtype_results <- rbindlist(subtype_results, fill = TRUE)
subtype_results[, FDR := p.adjust(p_value, "BH"), by = outcome]
fwrite(subtype_results, file.path(out, "pml_subtype_effects.csv"))

# Histology trend models. Patient fixed effects absorb cohort; flow cell and TIN are retained.
grade_results <- list()
for (outcome in outcomes) {
  fit <- cluster_fit(d, outcome, "grade_ordinal + median_tin_num + flow_cell_id + patient_id")
  grade_results[[length(grade_results) + 1]] <- contrast_result(
    fit, c(grade_ordinal = 1), "per_one_grade_increase", outcome, scale_map[[outcome]])
}
grade_results <- rbindlist(grade_results)
grade_results[, FDR := p.adjust(p_value, "BH")]
fwrite(grade_results, file.path(out, "pml_grade_trends.csv"))

cluster_cor <- function(data, method, stratum, B = 2000) {
  z <- data[complete.cases(data[, c("ABNORMAL_REPAIR", "CILIA_CONSENSUS", "patient_id")]), ]
  obs <- cor(z$ABNORMAL_REPAIR, z$CILIA_CONSENSUS, method = method)
  ids <- unique(z$patient_id)
  boot <- replicate(B, {
    draw <- sample(ids, length(ids), replace = TRUE)
    idx <- unlist(lapply(draw, function(id) which(z$patient_id == id)), use.names = FALSE)
    suppressWarnings(cor(z$ABNORMAL_REPAIR[idx], z$CILIA_CONSENSUS[idx], method = method))
  })
  boot <- boot[is.finite(boot)]
  se <- sd(boot)
  stat <- obs / se
  p <- 2 * pnorm(abs(stat), lower.tail = FALSE)
  patient_mean <- aggregate(cbind(ABNORMAL_REPAIR, CILIA_CONSENSUS) ~ patient_id, z, mean)
  patient_mean_r <- if (nrow(patient_mean) >= 4) cor(patient_mean$ABNORMAL_REPAIR,
                                                    patient_mean$CILIA_CONSENSUS,
                                                    method = method) else NA_real_
  data.frame(stratum = stratum, method = method, n_samples = nrow(z),
             n_patients = length(ids), correlation = obs,
             CI95_lower = quantile(boot, 0.025), CI95_upper = quantile(boot, 0.975),
             bootstrap_se = se, statistic = stat, p_value = p,
             patient_mean_correlation_sensitivity = patient_mean_r,
             inference = "patient-cluster bootstrap", stringsAsFactors = FALSE)
}

strata <- list(overall = d,
               Proliferative = d[d$molecular_subtype == "Proliferative", , drop = FALSE],
               other_subtypes = d[d$molecular_subtype != "Proliferative", , drop = FALSE])
cor_results <- rbindlist(lapply(names(strata), function(s) rbind(
  cluster_cor(strata[[s]], "pearson", s),
  cluster_cor(strata[[s]], "spearman", s)
)))
cor_results[, FDR := p.adjust(p_value, "BH")]
fwrite(cor_results, file.path(out, "repair_cilia_correlations.csv"))

score_cols <- c("sample_id", "sample", "patient_id", "grade", "grade_ordinal",
                "molecular_subtype", "flow_cell_id", "median_tin_num",
                "CILIA_CONSENSUS", "ABNORMAL_REPAIR", "FOXJ1", "KRT14", "KRT17",
                "ABNORMAL_REPAIR_Z", "CILIA_CONSENSUS_Z", "Repair_Cilia_Imbalance",
                "ABNORMAL_REPAIR_RANKN", "CILIA_CONSENSUS_RANKN",
                "Repair_Cilia_Imbalance_rank")
fwrite(d[, score_cols], file.path(out, "pml_extended_scores.csv"))

availability <- data.frame(
  field = c("participant_id", "lesion_id", "pathology_grade", "molecular_subtype",
            "longitudinal_timepoint", "progressive_persistent_regressive", "age", "sex",
            "smoking_status", "pack_years", "batch_flow_cell", "anatomical_site"),
  available = c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE),
  source = c("GEO characteristics", "not present", "GEO characteristics", "GEO characteristics",
             "not present", "not present", "not present", "not present", "not present",
             "not present", "GEO characteristics", "not present"),
  notes = c("50 total; 49 represented among passed-QC biopsies",
            "sample_id is unique but is not a validated lesion identifier",
            "normal through severe dysplasia; 4 unknown among passed-QC biopsies",
            "Normal/Inflammatory/Proliferative/Secretory",
            "GSE109743 is cross-sectional in public metadata",
            "Requires longitudinal dataset such as GSE79210",
            rep("", 6)), stringsAsFactors = FALSE)
fwrite(availability, file.path(out, "outcome_covariate_availability.csv"))

# Effect-size figure: all non-reference subtypes versus Normal.
plotdat <- subtype_results[contrast %in% c("Inflammatory_vs_Normal", "Proliferative_vs_Normal",
                                          "Secretory_vs_Normal") &
                           outcome %in% c("CILIA_CONSENSUS", "ABNORMAL_REPAIR", "FOXJ1", "KRT14", "KRT17")]
png(file.path(out, "figures/subtype_effects_forest.png"), width = 2400, height = 1700, res = 260)
par(mar = c(5, 13, 3, 2))
ord <- order(plotdat$outcome, plotdat$contrast)
p <- plotdat[ord]
y <- seq_len(nrow(p))
xlim <- range(c(p$CI95_lower, p$CI95_upper), finite = TRUE)
plot(p$estimate, y, xlim = xlim, ylim = c(0.5, length(y) + 0.5), pch = 19,
     yaxt = "n", xlab = "Adjusted effect estimate (95% CI)", ylab = "",
     main = "GSE109743: subtype effects with patient-cluster robust inference")
segments(p$CI95_lower, y, p$CI95_upper, y, lwd = 2)
axis(2, at = y, labels = paste(p$outcome, p$contrast, sep = " | "), las = 2, cex.axis = 0.75)
abline(v = 0, lty = 2, col = "grey40")
dev.off()

png(file.path(out, "figures/repair_vs_cilia_scatter.png"), width = 2100, height = 700, res = 220)
par(mfrow = c(1, 3), mar = c(5, 5, 3, 1))
for (s in names(strata)) {
  z <- strata[[s]]
  plot(z$CILIA_CONSENSUS, z$ABNORMAL_REPAIR, pch = 16, col = "#4C78A866",
       xlab = "CILIA_CONSENSUS", ylab = "ABNORMAL_REPAIR", main = s)
  abline(lm(ABNORMAL_REPAIR ~ CILIA_CONSENSUS, data = z), col = "#E45756", lwd = 2)
}
dev.off()

prolif <- subtype_results[outcome %in% c("CILIA_CONSENSUS", "ABNORMAL_REPAIR", "FOXJ1", "KRT14", "KRT17",
                                         "Repair_Cilia_Imbalance", "Repair_Cilia_Imbalance_rank") &
                          contrast == "Proliferative_vs_Normal"]
grade_report <- grade_results[outcome %in% c("CILIA_CONSENSUS", "ABNORMAL_REPAIR", "FOXJ1", "KRT14", "KRT17",
                                               "Repair_Cilia_Imbalance", "Repair_Cilia_Imbalance_rank")]
report <- c(
  "# GSE109743 extended PML analysis", "",
  sprintf("Passed-QC endobronchial biopsies: %d; patients: %d.", nrow(d), length(unique(d$patient_id))),
  "Primary inference uses patient fixed effects plus flow-cell and median-TIN adjustment, with patient-cluster robust standard errors.",
  "Patient fixed effects absorb cohort because cohort is constant within patient.", "",
  "## Proliferative versus Normal", "",
  paste0("- ", prolif$outcome, ": estimate=", sprintf("%.3f", prolif$estimate),
         ", 95% CI ", sprintf("%.3f", prolif$CI95_lower), " to ", sprintf("%.3f", prolif$CI95_upper),
         ", P=", format(prolif$p_value, digits = 3), ", FDR=", format(prolif$FDR, digits = 3)), "",
  "## Per one-grade increase in dysplasia severity", "",
  paste0("- ", grade_report$outcome, ": estimate=", sprintf("%.3f", grade_report$estimate),
         ", 95% CI ", sprintf("%.3f", grade_report$CI95_lower), " to ", sprintf("%.3f", grade_report$CI95_upper),
         ", P=", format(grade_report$p_value, digits = 3), ", FDR=", format(grade_report$FDR, digits = 3)), "",
  "## Longitudinal outcome", "",
  "No progressive/persistent/regressive field or public time point exists in GSE109743. No longitudinal claim is made from this dataset; GSE79210 is analyzed separately.", "",
  "## Caveats", "",
  "Effects are associations in study-supplied residual expression. Bulk signatures mix epithelial composition and within-cell-state expression. Cluster-bootstrap correlation P values are normal approximations from bootstrap SEs; patient-mean correlations are included as sensitivity estimates."
)
writeLines(report, file.path(out, "analysis_summary.md"))

capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
cat(sprintf("Completed GSE109743 extended analysis: %d biopsies, %d patients\n",
            nrow(d), length(unique(d$patient_id))))
