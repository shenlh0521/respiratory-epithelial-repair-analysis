#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

set.seed(20260811)
project <- Sys.getenv("PROJECT_ROOT", unset = normalizePath(".", mustWork = FALSE))
legacy <- file.path(project, "files-mentioned-by-the-user-users/outputs/lung_bulla_cilia_cancer")
base <- file.path(project, "01_PML_extended/GSE114489")
raw <- file.path(base, "raw")
out <- file.path(base, "results")
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)

parse_series <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  get_meta <- function(prefix) {
    hit <- lines[grepl(paste0("^", prefix, "\\t"), lines)]
    lapply(hit, function(x) gsub('^"|"$', '', strsplit(x, "\t", fixed = TRUE)[[1]][-1]))
  }
  gsm <- get_meta("!Sample_geo_accession")[[1]]
  md <- data.frame(gsm = gsm, title = get_meta("!Sample_title")[[1]], stringsAsFactors = FALSE)
  chars <- get_meta("!Sample_characteristics_ch1")
  for (z in chars) {
    key <- tolower(trimws(sub(":.*$", "", z[1])))
    key <- gsub("[^a-z0-9]+", "_", key)
    md[[key]] <- trimws(sub("^[^:]+:", "", z))
  }
  begin <- which(lines == "!series_matrix_table_begin") + 1L
  end <- which(lines == "!series_matrix_table_end") - 1L
  tab <- fread(text = paste(lines[begin:end], collapse = "\n"), check.names = FALSE,
               data.table = FALSE)
  names(tab)[1] <- "probe_id"
  names(tab)[-1] <- gsub('^"|"$', '', names(tab)[-1])
  list(metadata = md, expression = tab)
}

read_annotation <- function(path) {
  x <- fread(cmd = sprintf("gzip -cd %s", shQuote(path)), skip = "ID\tGene title",
             data.table = FALSE, check.names = FALSE)
  id_col <- intersect(c("ID", "ID_REF"), names(x))[1]
  sym_col <- intersect(c("Gene symbol", "Gene Symbol", "GENE_SYMBOL"), names(x))[1]
  if (is.na(id_col) || is.na(sym_col)) stop("GPL annotation lacks ID or gene symbol column")
  ans <- x[, c(id_col, sym_col), drop = FALSE]
  names(ans) <- c("probe_id", "symbol")
  ans$symbol <- toupper(trimws(sub(" ///.*$", "", ans$symbol)))
  ans <- ans[nzchar(ans$symbol) & ans$symbol != "---", , drop = FALSE]
  ans
}

series <- parse_series(file.path(raw, "GSE114489_series_matrix.txt.gz"))
md <- series$metadata
expr <- series$expression
ann <- read_annotation(file.path(raw, "GPL6244.annot.gz"))
expr <- merge(ann, expr, by = "probe_id")
sample_cols <- intersect(md$gsm, names(expr))
stopifnot(length(sample_cols) == nrow(md))

# Gene-level expression is the mean across transcript-cluster rows mapping to a symbol.
gene_dt <- as.data.table(expr[, c("symbol", sample_cols), drop = FALSE])
gene_dt <- gene_dt[, lapply(.SD, mean, na.rm = TRUE), by = symbol, .SDcols = sample_cols]
gene_mat <- as.matrix(gene_dt[, -1])
rownames(gene_mat) <- gene_dt$symbol
storage.mode(gene_mat) <- "numeric"
gene_z <- t(scale(t(gene_mat)))
gene_z[!is.finite(gene_z)] <- NA_real_

legacy_sets <- yaml.load_file(file.path(legacy, "config/gene_sets.yaml"))$gene_sets
gene_sets <- list(
  DNA_DAMAGE_RESPONSE = legacy_sets$DNA_DAMAGE_SENTINELS,
  ABNORMAL_REPAIR = legacy_sets$ABNORMAL_REPAIR,
  KRT14_KRT17_REPAIR = c("KRT14", "KRT17"),
  CILIA_CONSENSUS = legacy_sets$CILIA_CONSENSUS,
  MULTICILIOGENESIS = unique(c(legacy_sets$CILIA_REGULATORY, legacy_sets$DEUTEROSOMAL)),
  MATURE_CILIATED = c("FOXJ1", "PIFO", "TPPP3", "CAPS", "RSPH1"),
  KRT4_KRT13_TRANSITIONAL_FATE = c("KRT4", "KRT13")
)
scores <- sapply(gene_sets, function(gs) {
  hit <- intersect(toupper(gs), rownames(gene_z))
  if (length(hit) < 2) rep(NA_real_, ncol(gene_z)) else colMeans(gene_z[hit, , drop = FALSE], na.rm = TRUE)
})
scores <- as.data.frame(scores, check.names = FALSE)
scores$gsm <- rownames(scores)
md <- merge(md, scores, by = "gsm", all.x = TRUE, sort = FALSE)
for (g in c("FOXJ1", "KRT14", "KRT17", "MCIDAS", "GMNC", "KRT4", "KRT13")) {
  md[[g]] <- if (g %in% rownames(gene_z)) gene_z[g, md$gsm] else NA_real_
}

num_prefix <- function(x) suppressWarnings(as.numeric(sub("\\.1$", "", x)))
md$group <- as.integer(sub("Group ([1-4]).*$", "\\1", md$title))
md$lesion_site_id <- as.integer(gsub("[^0-9]", "", sub("^.*Biopsy", "", md$title)))
md$baseline_grade <- num_prefix(md$bl_frozen_dx)
md$baseline_ffpe_grade <- num_prefix(md$bl_ffpe_dx)
md$followup_grade <- num_prefix(md$f_u_bx_dx)
md$age_num <- suppressWarnings(as.numeric(md$age))
md$pack_years_num <- suppressWarnings(as.numeric(md$smoking_pack_yr))
md$outcome4 <- factor(md$group, levels = 1:4,
                      labels = c("persistent_dysplasia", "regressive_dysplasia",
                                 "progressive_nondysplasia", "stable_nondysplasia"))
md$outcome_binary <- as.integer(md$group %in% c(1, 3))
zscore <- function(x) as.numeric(scale(x))
rank_normal <- function(x) {
  ans <- rep(NA_real_, length(x)); ok <- is.finite(x); n <- sum(ok)
  ans[ok] <- qnorm((rank(x[ok], ties.method = "average") - 0.5) / n); ans
}
md$Repair_Cilia_Imbalance <- zscore(md$ABNORMAL_REPAIR) - zscore(md$CILIA_CONSENSUS)
md$Repair_Cilia_Imbalance_rank <- rank_normal(md$ABNORMAL_REPAIR) - rank_normal(md$CILIA_CONSENSUS)

analysis_vars <- c(names(gene_sets), "Repair_Cilia_Imbalance", "Repair_Cilia_Imbalance_rank",
                   "FOXJ1", "KRT14", "KRT17", "MCIDAS", "GMNC", "KRT4", "KRT13")

continuous_model <- function(d, variable, contrast, group_col = NULL, adjusted = TRUE) {
  vars <- unique(c(variable, group_col, "baseline_grade", "age_num", "pack_years_num"))
  z <- d[complete.cases(d[, vars, drop = FALSE]), , drop = FALSE]
  model_label <- if (adjusted) "adjusted_baseline_grade_age_packyears" else "unadjusted"
  if (nrow(z) < 8 || !is.finite(sd(z[[variable]])) || sd(z[[variable]]) == 0 ||
      (!is.null(group_col) && length(unique(z[[group_col]])) < 2)) {
    return(data.frame(variable = variable, contrast = contrast, model = model_label,
                      n_samples = nrow(z), estimate_SD = NA_real_, std_error = NA_real_,
                      CI95_lower = NA_real_, CI95_upper = NA_real_, statistic = NA_real_,
                      p_value = NA_real_, inference = "not estimable; insufficient observed data",
                      stringsAsFactors = FALSE))
  }
  z$outcome_z <- zscore(z[[variable]])
  if (!is.null(group_col)) z$case <- as.integer(z[[group_col]])
  form <- if (!is.null(group_col) && adjusted) outcome_z ~ case + baseline_grade + age_num + pack_years_num else
    if (!is.null(group_col)) outcome_z ~ case else
      if (adjusted) outcome_z ~ baseline_grade + age_num + pack_years_num else outcome_z ~ baseline_grade
  fit <- lm(form, data = z)
  term <- if (!is.null(group_col)) "case" else "baseline_grade"
  co <- summary(fit)$coefficients[term, ]
  ci <- confint(fit, term, level = 0.95)
  data.frame(variable = variable, contrast = contrast,
             model = model_label,
             n_samples = nobs(fit), estimate_SD = unname(co["Estimate"]),
             std_error = unname(co["Std. Error"]), CI95_lower = ci[1], CI95_upper = ci[2],
             statistic = unname(co["t value"]), p_value = unname(co["Pr(>|t|)"]),
             inference = "lesion-site level; patient identifier unavailable", stringsAsFactors = FALSE)
}

primary <- md[md$group %in% c(1, 2), , drop = FALSE]
primary$case_primary <- as.integer(primary$group == 1)
primary_effects <- rbindlist(lapply(analysis_vars, function(v) rbind(
  continuous_model(primary, v, "persistent_vs_regressive_dysplasia", "case_primary", FALSE),
  continuous_model(primary, v, "persistent_vs_regressive_dysplasia", "case_primary", TRUE)
)))
primary_effects[, FDR := p.adjust(p_value, "BH"), by = model]
fwrite(primary_effects, file.path(out, "persistent_vs_regressive_effects.csv"))

secondary <- md
secondary$case_combined <- md$outcome_binary
secondary_effects <- rbindlist(lapply(analysis_vars, function(v)
  continuous_model(secondary, v, "persistent_or_progressive_vs_regressive_or_stable", "case_combined", TRUE)))
secondary_effects[, FDR := p.adjust(p_value, "BH")]
fwrite(secondary_effects, file.path(out, "combined_outcome_sensitivity_effects.csv"))

severity_effects <- rbindlist(lapply(analysis_vars, function(v)
  continuous_model(md, v, "per_one_baseline_histology_grade", NULL, TRUE)))
severity_effects[, FDR := p.adjust(p_value, "BH")]
fwrite(severity_effects, file.path(out, "baseline_histology_trends.csv"))

logistic_one <- function(d, variable, adjusted) {
  vars <- c(variable, "case_primary", "baseline_grade", "age_num", "pack_years_num")
  z <- d[complete.cases(d[, vars, drop = FALSE]), , drop = FALSE]
  model_label <- if (adjusted) "adjusted_baseline_grade_age_packyears" else "unadjusted"
  if (nrow(z) < 8 || length(unique(z$case_primary)) < 2 ||
      !is.finite(sd(z[[variable]])) || sd(z[[variable]]) == 0) {
    return(data.frame(variable = variable, model = model_label, n_sites = nrow(z),
                      events_persistent = sum(z$case_primary), log_OR_per_SD = NA_real_,
                      OR_per_SD = NA_real_, CI95_lower = NA_real_, CI95_upper = NA_real_,
                      p_value = NA_real_, converged = FALSE,
                      inference = "not estimable; insufficient observed data", stringsAsFactors = FALSE))
  }
  z$score_std <- zscore(z[[variable]])
  form <- if (adjusted) case_primary ~ score_std + baseline_grade + age_num + pack_years_num else case_primary ~ score_std
  fit <- glm(form, family = binomial(), data = z)
  co <- summary(fit)$coefficients["score_std", ]
  est <- unname(co["Estimate"]); se <- unname(co["Std. Error"])
  data.frame(variable = variable, model = model_label,
             n_sites = nobs(fit), events_persistent = sum(z$case_primary), log_OR_per_SD = est,
             OR_per_SD = exp(est), CI95_lower = exp(est - 1.96 * se), CI95_upper = exp(est + 1.96 * se),
             p_value = unname(co["Pr(>|z|)"]), converged = fit$converged,
             inference = "lesion-site logistic model; patient identifier unavailable", stringsAsFactors = FALSE)
}
logistic <- rbindlist(lapply(analysis_vars, function(v) rbind(logistic_one(primary, v, FALSE), logistic_one(primary, v, TRUE))))
logistic[, FDR := p.adjust(p_value, "BH"), by = model]
fwrite(logistic, file.path(out, "baseline_future_outcome_logistic.csv"))

coverage <- rbindlist(lapply(names(gene_sets), function(nm) {
  req <- toupper(gene_sets[[nm]]); hit <- intersect(req, rownames(gene_mat))
  data.frame(signature = nm, requested_n = length(req), observed_n = length(hit),
             observed_genes = paste(hit, collapse = ";"), missing_genes = paste(setdiff(req, hit), collapse = ";"))
}))
fwrite(coverage, file.path(out, "signature_coverage.csv"))
fwrite(md, file.path(out, "metadata_scores.csv"))

plotdat <- primary_effects[model == "adjusted_baseline_grade_age_packyears" &
                           variable %in% c("ABNORMAL_REPAIR", "KRT14_KRT17_REPAIR", "CILIA_CONSENSUS",
                                           "MULTICILIOGENESIS", "MATURE_CILIATED", "Repair_Cilia_Imbalance",
                                           "KRT4_KRT13_TRANSITIONAL_FATE", "FOXJ1")]
png(file.path(out, "figures", "GSE114489_outcome_effects.png"), width = 2200, height = 1450, res = 260)
par(mar = c(5, 12, 3, 2))
plotdat <- plotdat[order(estimate_SD)]
y <- seq_len(nrow(plotdat)); xr <- range(c(plotdat$CI95_lower, plotdat$CI95_upper), finite = TRUE)
plot(plotdat$estimate_SD, y, xlim = xr, ylim = c(.5, length(y) + .5), yaxt = "n", pch = 19,
     xlab = "Adjusted standardized mean difference: persistent minus regressive (95% CI)", ylab = "",
     main = "GSE114489 baseline dysplasia and future lesion fate")
segments(plotdat$CI95_lower, y, plotdat$CI95_upper, y, lwd = 2)
axis(2, y, plotdat$variable, las = 2, cex.axis = .8); abline(v = 0, lty = 2, col = "grey40")
dev.off()

imb <- primary_effects[variable == "Repair_Cilia_Imbalance" & model == "adjusted_baseline_grade_age_packyears"]
orimb <- logistic[variable == "Repair_Cilia_Imbalance" & model == "adjusted_baseline_grade_age_packyears"]
writeLines(c(
  "# GSE114489 third-stage longitudinal PML analysis", "",
  sprintf("Baseline biopsy sites: %d; persistent dysplasia: %d; regressive dysplasia: %d; progressive nondysplasia: %d; stable nondysplasia: %d.",
          nrow(md), sum(md$group == 1), sum(md$group == 2), sum(md$group == 3), sum(md$group == 4)),
  "The public series does not expose a participant identifier. All inference is therefore lesion-site level and cannot adjust standard errors for multiple sites from the same person.", "",
  "## Primary lesion-fate test", "",
  sprintf("Adjusted persistent-minus-regressive Repair_Cilia_Imbalance = %.3f SD (95%% CI %.3f to %.3f), P=%.3g, FDR=%.3g.",
          imb$estimate_SD, imb$CI95_lower, imb$CI95_upper, imb$p_value, imb$FDR),
  sprintf("Adjusted OR for persistence per 1-SD baseline imbalance = %.3f (95%% CI %.3f to %.3f), P=%.3g, FDR=%.3g.",
          orimb$OR_per_SD, orimb$CI95_lower, orimb$CI95_upper, orimb$p_value, orimb$FDR), "",
  "Models adjust for baseline frozen histology grade, age and pack-years. Sex, smoking status, anatomical site and batch were not available as reusable covariates in GEO.",
  "The combined persistent/progressive versus regressive/stable comparison is secondary because it mixes baseline dysplastic and nondysplastic lesions.",
  "Associations support risk stratification but do not establish causality."
), file.path(out, "analysis_summary.md"))
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
cat(sprintf("Completed GSE114489: %d sites; primary %d persistent vs %d regressive\n",
            nrow(md), sum(primary$case_primary), sum(primary$case_primary == 0)))
