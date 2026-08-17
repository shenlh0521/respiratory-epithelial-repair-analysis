#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(yaml)
})

set.seed(20260811)
project <- Sys.getenv("PROJECT_ROOT", unset = normalizePath(".", mustWork = FALSE))
legacy <- file.path(project, "files-mentioned-by-the-user-users/outputs/lung_bulla_cilia_cancer")
data_dir <- file.path(project, "01_PML_extended/GSE79210/data")
out <- file.path(project, "01_PML_extended/GSE79210/results")
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
gene_sets <- yaml.load_file(file.path(legacy, "config/gene_sets.yaml"))$gene_sets

parse_geo_metadata <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  meta <- lines[grepl("^!Sample_", lines)]
  get_rows <- function(key) {
    hit <- meta[grepl(paste0("^!Sample_", key, "\\t"), meta)]
    lapply(hit, function(x) {
      z <- strsplit(x, "\t", fixed = TRUE)[[1]][-1]
      gsub('^"|"$', '', z)
    })
  }
  ids <- get_rows("geo_accession")[[1]]
  ans <- data.frame(gsm = ids,
                    title = get_rows("title")[[1]],
                    sample_id = get_rows("description")[[1]],
                    stringsAsFactors = FALSE)
  chars <- get_rows("characteristics_ch1")
  for (z in chars) {
    key <- tolower(trimws(sub(":.*$", "", z[1])))
    key <- gsub("[^a-z0-9]+", "_", key)
    ans[[key]] <- trimws(sub("^[^:]+:", "", z))
  }
  ans
}

md <- parse_geo_metadata(file.path(data_dir, "GSE79210_series_matrix.txt.gz"))
counts <- fread(cmd = sprintf("gzip -cd %s", shQuote(file.path(data_dir, "GSE79210_P_counts.tsv.gz"))),
                check.names = FALSE, data.table = FALSE)
names(counts)[1] <- "ensembl_id"
counts$ensembl_id <- sub("\\.[0-9]+$", "", counts$ensembl_id)
mat <- as.matrix(counts[, -1, drop = FALSE])
storage.mode(mat) <- "numeric"
rownames(mat) <- counts$ensembl_id
stopifnot(all(md$sample_id %in% colnames(mat)))
mat <- mat[, md$sample_id, drop = FALSE]

map <- fread(file.path(legacy, "data/metadata/target_symbol_ensembl.csv"), data.table = FALSE)
map <- map[map$ensembl_id %in% rownames(mat), , drop = FALSE]
dge <- DGEList(mat)
dge <- calcNormFactors(dge)
all_logcpm <- cpm(dge, log = TRUE, prior.count = 1)
logcpm <- all_logcpm[map$ensembl_id, , drop = FALSE]
rownames(logcpm) <- map$symbol
if (anyDuplicated(rownames(logcpm))) {
  av <- rowMeans(logcpm)
  ord <- order(rownames(logcpm), -av)
  logcpm <- logcpm[ord, , drop = FALSE]
  logcpm <- logcpm[!duplicated(rownames(logcpm)), , drop = FALSE]
}
gene_z <- t(scale(t(logcpm)))
gene_z[!is.finite(gene_z)] <- NA_real_
score <- sapply(gene_sets, function(genes) {
  hit <- intersect(toupper(genes), rownames(gene_z))
  if (length(hit) < 2) rep(NA_real_, ncol(gene_z)) else colMeans(gene_z[hit, , drop = FALSE], na.rm = TRUE)
})
score <- as.data.frame(score, check.names = FALSE)
score$sample_id <- rownames(score)
md <- merge(md, score, by = "sample_id", all.x = TRUE, sort = FALSE)

for (gene in c("FOXJ1", "KRT14", "KRT17")) {
  md[[gene]] <- if (gene %in% rownames(logcpm)) as.numeric(logcpm[gene, md$sample_id]) else NA_real_
}
md$time_num <- as.integer(sub("[^0-9]", "", md$time_point))
md$age_num <- as.numeric(md$age)
md$pack_years_num <- as.numeric(md$pack_years)
md$sex <- factor(md$sex)
md$patient_id <- factor(md$patient_id)
md$smoking_status <- factor(md$smoking_status)
md$copd_status <- factor(md$copd_status)
grade_map <- c("Normal" = 0, "Hyperplasia" = 1, "Metaplasia" = 2,
               "Mild dysplasia" = 3, "Moderate dysplasia" = 4,
               "Severe dysplasia" = 5, "Carcinoma in situ" = 6)
md$grade_ordinal <- unname(grade_map[md$max_histology])
md$dysplasia_binary <- as.integer(md$dysplasia_status == "Dysplasia")

zscore <- function(x) as.numeric(scale(x))
rank_normal <- function(x) {
  ans <- rep(NA_real_, length(x)); ok <- is.finite(x); n <- sum(ok)
  ans[ok] <- qnorm((rank(x[ok], ties.method = "average") - 0.5) / n); ans
}
md$Repair_Cilia_Imbalance <- zscore(md$ABNORMAL_REPAIR) - zscore(md$CILIA_CONSENSUS)
md$Repair_Cilia_Imbalance_rank <- rank_normal(md$ABNORMAL_REPAIR) - rank_normal(md$CILIA_CONSENSUS)

# Derive subject-level future outcome from public repeated max-histology/dysplasia fields.
md <- md[order(md$patient_id, md$time_num), , drop = FALSE]
outcome <- rbindlist(lapply(split(md, md$patient_id, drop = TRUE), function(z) {
  b <- z[which.min(z$time_num), , drop = FALSE]
  l <- z[which.max(z$time_num), , drop = FALSE]
  category <- if (b$dysplasia_binary == 1 && l$dysplasia_binary == 0) "regressive" else
    if (b$dysplasia_binary == 1 && l$dysplasia_binary == 1) "persistent" else
      if (b$dysplasia_binary == 0 && l$dysplasia_binary == 1) "progressive" else "stable_normal"
  data.frame(patient_id = as.character(b$patient_id), baseline_sample_id = b$sample_id,
             baseline_gsm = b$gsm, baseline_time = b$time_num, last_time = l$time_num,
             baseline_histology = b$max_histology, last_histology = l$max_histology,
             baseline_dysplasia = b$dysplasia_binary, last_dysplasia = l$dysplasia_binary,
             grade_change = l$grade_ordinal - b$grade_ordinal,
             outcome = category, stringsAsFactors = FALSE)
}))
md <- merge(md, outcome[, c("patient_id", "outcome")], by = "patient_id", all.x = TRUE, sort = FALSE)

analysis_vars <- c("CILIA_CONSENSUS", "ABNORMAL_REPAIR", "FOXJ1", "KRT14", "KRT17",
                   "Repair_Cilia_Imbalance", "Repair_Cilia_Imbalance_rank")
baseline <- md[md$time_num == 1, , drop = FALSE]
baseline <- merge(baseline, outcome, by = "patient_id", suffixes = c("", "_outcome"), sort = FALSE)
primary <- baseline[baseline$baseline_dysplasia == 1 & baseline$outcome %in% c("persistent", "regressive"), ]
primary$outcome_binary <- as.integer(primary$outcome == "persistent")

group_effect <- function(data, variable) {
  x <- data[data$outcome == "persistent" & is.finite(data[[variable]]), variable]
  y <- data[data$outcome == "regressive" & is.finite(data[[variable]]), variable]
  if (length(x) < 2 || length(y) < 2) {
    return(data.frame(variable = variable, contrast = "persistent_vs_regressive",
                      n_persistent = length(x), n_regressive = length(y),
                      mean_persistent = NA_real_, mean_regressive = NA_real_,
                      mean_difference = NA_real_, CI95_lower = NA_real_, CI95_upper = NA_real_,
                      hedges_g = NA_real_, p_value = NA_real_, stringsAsFactors = FALSE))
  }
  tt <- t.test(x, y)
  nx <- length(x); ny <- length(y)
  sp <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  d <- (mean(x) - mean(y)) / sp
  J <- 1 - 3 / (4 * (nx + ny) - 9)
  data.frame(variable = variable, contrast = "persistent_vs_regressive",
             n_persistent = nx, n_regressive = ny,
             mean_persistent = mean(x), mean_regressive = mean(y),
             mean_difference = mean(x) - mean(y),
             CI95_lower = tt$conf.int[1], CI95_upper = tt$conf.int[2],
             hedges_g = J * d, p_value = tt$p.value, stringsAsFactors = FALSE)
}
group_results <- rbindlist(lapply(analysis_vars, function(v) group_effect(primary, v)))
group_results[, FDR := p.adjust(p_value, "BH")]
fwrite(group_results, file.path(out, "baseline_persistent_vs_regressive.csv"))

logistic_one <- function(data, variable, adjusted = FALSE) {
  z <- data[is.finite(data[[variable]]), , drop = FALSE]
  if (nrow(z) < 8 || length(unique(z$outcome_binary)) < 2) {
    return(data.frame(variable = variable,
                      model = if (adjusted) "adjusted_age_sex_packyears_exploratory" else "unadjusted_primary",
                      n_subjects = nrow(z), events_persistent = sum(z$outcome_binary),
                      log_OR_per_SD = NA_real_, OR_per_SD = NA_real_, CI95_lower = NA_real_,
                      CI95_upper = NA_real_, p_value = NA_real_, converged = FALSE,
                      stringsAsFactors = FALSE))
  }
  z$score_std <- zscore(z[[variable]])
  form <- if (adjusted) outcome_binary ~ score_std + age_num + sex + pack_years_num else outcome_binary ~ score_std
  fit <- glm(form, family = binomial(), data = z)
  co <- summary(fit)$coefficients["score_std", ]
  est <- unname(co["Estimate"]); se <- unname(co["Std. Error"])
  data.frame(variable = variable,
             model = if (adjusted) "adjusted_age_sex_packyears_exploratory" else "unadjusted_primary",
             n_subjects = nobs(fit), events_persistent = sum(z$outcome_binary),
             log_OR_per_SD = est, OR_per_SD = exp(est),
             CI95_lower = exp(est - 1.96 * se), CI95_upper = exp(est + 1.96 * se),
             p_value = unname(co["Pr(>|z|)"]), converged = fit$converged,
             stringsAsFactors = FALSE)
}
logistic_results <- rbindlist(lapply(analysis_vars, function(v) rbind(
  logistic_one(primary, v, FALSE), logistic_one(primary, v, TRUE)
)))
logistic_results[, FDR := p.adjust(p_value, "BH"), by = model]
fwrite(logistic_results, file.path(out, "baseline_future_outcome_logistic.csv"))

# Within-subject longitudinal association with concurrent max histology.
cluster_fit <- function(data, outcome_col) {
  vars <- c(outcome_col, "grade_ordinal", "time_num", "patient_id")
  z <- data[complete.cases(data[, vars]), , drop = FALSE]
  if (nrow(z) < 8 || length(unique(z$patient_id)) < 4) {
    return(data.frame(variable = outcome_col, contrast = "within_subject_per_grade_increase",
                      n_samples = nrow(z), n_subjects = length(unique(z$patient_id)),
                      estimate = NA_real_, std_error = NA_real_, CI95_lower = NA_real_,
                      CI95_upper = NA_real_, statistic = NA_real_, p_value = NA_real_,
                      stringsAsFactors = FALSE))
  }
  X0 <- model.matrix(~ grade_ordinal + time_num + patient_id, data = z)
  q <- qr(X0); keep <- sort(q$pivot[seq_len(q$rank)]); X <- X0[, keep, drop = FALSE]
  y <- z[[outcome_col]]; beta <- solve(crossprod(X), crossprod(X, y)); e <- as.numeric(y - X %*% beta)
  bread <- solve(crossprod(X)); groups <- split(seq_len(nrow(z)), droplevels(z$patient_id), drop = TRUE)
  meat <- matrix(0, ncol(X), ncol(X))
  for (idx in groups) { u <- as.numeric(crossprod(X[idx, , drop = FALSE], e[idx])); meat <- meat + tcrossprod(u) }
  G <- length(groups); N <- nrow(z); P <- ncol(X)
  vc <- (G/(G-1))*((N-1)/(N-P))*bread%*%meat%*%bread
  j <- which(colnames(X) == "grade_ordinal")
  est <- beta[j]; se <- sqrt(vc[j,j]); crit <- qt(.975, G-1); stat <- est/se
  data.frame(variable = outcome_col, contrast = "within_subject_per_grade_increase",
             n_samples = N, n_subjects = G, estimate = est, std_error = se,
             CI95_lower = est - crit*se, CI95_upper = est + crit*se,
             statistic = stat, p_value = 2*pt(abs(stat), G-1, lower.tail = FALSE),
             stringsAsFactors = FALSE)
}
longitudinal_results <- rbindlist(lapply(analysis_vars, function(v) cluster_fit(md, v)))
longitudinal_results[, FDR := p.adjust(p_value, "BH")]
fwrite(longitudinal_results, file.path(out, "within_subject_histology_associations.csv"))

coverage <- data.frame(signature = names(gene_sets),
                       genes_requested = sapply(gene_sets, length),
                       genes_observed = sapply(gene_sets, function(x) sum(toupper(x) %in% rownames(logcpm))))
coverage$coverage <- coverage$genes_observed / coverage$genes_requested
fwrite(coverage, file.path(out, "signature_coverage.csv"))
fwrite(md, file.path(out, "metadata_scores_long.csv"))
fwrite(outcome, file.path(out, "subject_outcomes.csv"))

png(file.path(out, "figures/baseline_outcome_effects.png"), width = 2200, height = 1500, res = 260)
par(mar = c(5, 11, 3, 2))
p <- group_results[order(group_results$mean_difference)]
y <- seq_len(nrow(p)); xr <- range(c(p$CI95_lower, p$CI95_upper), finite = TRUE)
plot(p$mean_difference, y, xlim = xr, ylim = c(.5, length(y)+.5), yaxt = "n", pch = 19,
     xlab = "Baseline mean difference: persistent minus regressive (95% CI)", ylab = "",
     main = "GSE79210 longitudinal PML field")
segments(p$CI95_lower, y, p$CI95_upper, y, lwd = 2)
axis(2, y, p$variable, las = 2); abline(v = 0, lty = 2, col = "grey40")
dev.off()

png(file.path(out, "figures/subject_imbalance_trajectories.png"), width = 1800, height = 1300, res = 240)
par(mar = c(5, 5, 3, 1))
plot(range(md$time_num), range(md$Repair_Cilia_Imbalance), type = "n",
     xlab = "Bronchoscopy time point", ylab = "Repair–cilia imbalance",
     main = "GSE79210 subject trajectories")
for (pid in unique(md$patient_id)) {
  z <- md[md$patient_id == pid, ]
  col <- if (unique(z$outcome) == "regressive") "#4C78A8" else if (unique(z$outcome) == "persistent") "#E45756" else "#999999"
  lines(z$time_num, z$Repair_Cilia_Imbalance, type = "b", col = paste0(col, "99"), pch = 16)
}
legend("topright", legend = c("persistent", "regressive", "stable normal"),
       col = c("#E45756", "#4C78A8", "#999999"), lty = 1, pch = 16, bty = "n")
dev.off()

key <- group_results[group_results$variable == "Repair_Cilia_Imbalance", ]
logkey <- logistic_results[logistic_results$variable == "Repair_Cilia_Imbalance" & logistic_results$model == "unadjusted_primary", ]
summary_lines <- c(
  "# GSE79210 longitudinal PML analysis", "",
  sprintf("Processed RNA-seq samples: %d; subjects: %d; repeated time points: 2 for 18 subjects and 3 for 5 subjects.", nrow(md), length(unique(md$patient_id))),
  "Expression is from bronchial brushing; public max histology summarizes the subject's contemporaneous lesion burden. This is a subject-level airway-field analysis, not same-lesion biopsy expression.", "",
  "## Outcome definition", "",
  sprintf("Among baseline-dysplastic subjects: persistent=%d, regressive=%d, progressive=%d.",
          sum(primary$outcome == "persistent"), sum(primary$outcome == "regressive"),
          sum(primary$outcome == "progressive")),
  "Persistent means dysplasia at baseline and last visit; regressive means dysplasia at baseline and no dysplasia at last visit. No baseline-dysplastic subject met a separate progressive category in the public binary dysplasia field.", "",
  "## Repair–cilia imbalance", "",
  sprintf("Baseline persistent minus regressive mean difference=%.3f (95%% CI %.3f to %.3f), Hedges g=%.3f, P=%.3g, FDR=%.3g.", key$mean_difference, key$CI95_lower, key$CI95_upper, key$hedges_g, key$p_value, key$FDR),
  sprintf("Unadjusted future persistent/regressive OR per 1-SD baseline imbalance=%.3f (95%% CI %.3f to %.3f), P=%.3g, FDR=%.3g.", logkey$OR_per_SD, logkey$CI95_lower, logkey$CI95_upper, logkey$p_value, logkey$FDR), "",
  "## Inference limits", "",
  "The primary endpoint has only 18 baseline-dysplastic subjects. Age/sex/pack-years adjusted logistic models are explicitly exploratory and may be unstable. COPD status was not added to that multivariable model because sparse/unknown levels would overfit. Associations do not establish that the airway-field signature causes lesion persistence."
)
writeLines(summary_lines, file.path(out, "analysis_summary.md"))
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
cat(sprintf("Completed GSE79210: %d samples, %d subjects; persistent=%d regressive=%d\n",
            nrow(md), length(unique(md$patient_id)), sum(primary$outcome == "persistent"),
            sum(primary$outcome == "regressive")))
