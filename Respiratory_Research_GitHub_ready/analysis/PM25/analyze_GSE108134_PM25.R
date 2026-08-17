#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

set.seed(20260811)
project <- Sys.getenv("PROJECT_ROOT", unset = normalizePath(".", mustWork = FALSE))
root <- file.path(project, "03_PM25/GSE108134")
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
target_genes <- unique(c(unlist(gene_sets), "FOXJ1", "MCIDAS", "GMNC", paste0("RFX", 1:7),
                         "KRT14", "KRT17"))

parse_series_matrix <- function(path) {
  x <- readLines(gzfile(path), warn = FALSE)
  get_fields <- function(prefix) {
    lines <- x[startsWith(x, prefix)]
    lapply(lines, function(z) {
      v <- strsplit(z, "\t", fixed = TRUE)[[1]][-1]
      gsub('^"|"$', "", v)
    })
  }
  titles <- get_fields("!Sample_title")[[1]]
  gsm <- get_fields("!Sample_geo_accession")[[1]]
  chars <- get_fields("!Sample_characteristics_ch1")
  extract_char <- function(label) {
    hit <- vapply(chars, function(v) startsWith(tolower(v[1]), tolower(label)), logical(1))
    if (!any(hit)) return(rep(NA_character_, length(titles)))
    sub(paste0("^", label, "\\s*"), "", chars[[which(hit)[1]]], ignore.case = TRUE)
  }
  titles <- sub(" \\(demographics_pollution\\)$", "", titles)
  data.frame(sample_id = titles, gsm = gsm,
             smoking_status = extract_char("smoking status:"),
             pm25_30day = as.numeric(extract_char("30-day mean pm2.5 exposure level:")),
             bronchoscopy_month = extract_char("time of serial bronchoscopy (month):"),
             stringsAsFactors = FALSE)
}

cluster_fit <- function(data, outcome, exposure = "pm25_30day", cluster = "subject_id") {
  z <- data[complete.cases(data[, c(outcome, exposure, cluster), drop = FALSE]), , drop = FALSE]
  X <- model.matrix(as.formula(paste("~", exposure)), z)
  y <- z[[outcome]]
  beta <- drop(solve(crossprod(X), crossprod(X, y)))
  names(beta) <- colnames(X)
  resid <- as.numeric(y - X %*% beta)
  bread <- solve(crossprod(X))
  groups <- split(seq_len(nrow(z)), factor(z[[cluster]]), drop = TRUE)
  meat <- matrix(0, ncol(X), ncol(X))
  for (idx in groups) {
    u <- as.numeric(crossprod(X[idx, , drop = FALSE], resid[idx]))
    meat <- meat + tcrossprod(u)
  }
  G <- length(groups); N <- nrow(z); P <- ncol(X)
  vc <- (G / (G - 1)) * ((N - 1) / (N - P)) * bread %*% meat %*% bread
  term <- exposure
  est <- unname(beta[term]); se <- sqrt(vc[term, term]); df <- G - 1
  stat <- est / se; p <- 2 * pt(abs(stat), df = df, lower.tail = FALSE)
  crit <- qt(0.975, df = df)
  data.frame(n_samples = N, n_subjects = G, estimate_per_ug_m3 = est,
             std_error = se, CI95_lower = est - crit * se, CI95_upper = est + crit * se,
             statistic = stat, p_value = p, stringsAsFactors = FALSE)
}

metadata <- parse_series_matrix(file.path(data_dir, "GSE108134_series_matrix.txt.gz"))
metadata$subject_id <- sub("_M.*$", "", metadata$sample_id)
pm_sd <- sd(metadata$pm25_30day, na.rm = TRUE)

cohorts <- list(
  smoker = file.path(data_dir, "GSE108134_307SAE-S_PM2.5_resid.csv.gz"),
  nonsmoker = file.path(data_dir, "GSE108134_98SAE-NS_PM2.5_resid.csv.gz")
)

all_effects <- list(); all_scores <- list(); coverage <- list()
for (cohort in names(cohorts)) {
  path <- cohorts[[cohort]]
  if (!file.exists(path)) stop("Missing processed matrix: ", path)
  mat <- fread(cmd = sprintf("gzip -cd %s", shQuote(path)), check.names = FALSE,
               data.table = FALSE)
  names(mat)[1] <- "gene"
  mat$gene <- toupper(mat$gene)
  mat <- mat[mat$gene %in% target_genes, , drop = FALSE]
  if (anyDuplicated(mat$gene)) {
    mat <- aggregate(. ~ gene, data = mat, FUN = mean)
  }
  rownames(mat) <- mat$gene
  expr <- as.matrix(mat[, -1, drop = FALSE]); storage.mode(expr) <- "double"
  zexpr <- t(scale(t(expr)))
  zexpr[!is.finite(zexpr)] <- NA_real_
  score <- function(genes) {
    use <- intersect(genes, rownames(zexpr))
    if (!length(use)) return(rep(NA_real_, ncol(zexpr)))
    colMeans(zexpr[use, , drop = FALSE], na.rm = TRUE)
  }
  score_df <- data.frame(sample_id = colnames(expr), stringsAsFactors = FALSE)
  for (nm in names(gene_sets)) score_df[[nm]] <- score(gene_sets[[nm]])
  for (gene in intersect(target_genes, rownames(zexpr))) score_df[[gene]] <- zexpr[gene, ]
  score_df$ABNORMAL_REPAIR_Z <- as.numeric(scale(score_df$ABNORMAL_REPAIR))
  score_df$CILIA_CONSENSUS_Z <- as.numeric(scale(score_df$CILIA_CONSENSUS))
  score_df$Repair_Cilia_Imbalance <- score_df$ABNORMAL_REPAIR_Z - score_df$CILIA_CONSENSUS_Z
  score_df <- merge(score_df, metadata, by = "sample_id", all.x = TRUE, sort = FALSE)
  score_df$cohort <- cohort
  all_scores[[cohort]] <- score_df

  for (nm in names(gene_sets)) {
    coverage[[length(coverage) + 1]] <- data.frame(
      cohort = cohort, signature = nm, n_defined = length(gene_sets[[nm]]),
      n_detected = length(intersect(gene_sets[[nm]], rownames(expr))),
      detected_genes = paste(intersect(gene_sets[[nm]], rownames(expr)), collapse = ";"),
      missing_genes = paste(setdiff(gene_sets[[nm]], rownames(expr)), collapse = ";"))
  }
  outcomes <- c(names(gene_sets), "Repair_Cilia_Imbalance", "FOXJ1", "MCIDAS", "GMNC",
                intersect(paste0("RFX", 1:7), names(score_df)), "KRT14", "KRT17")
  outcomes <- unique(outcomes[outcomes %in% names(score_df)])
  eff <- rbindlist(lapply(outcomes, function(nm) {
    ans <- cluster_fit(score_df, nm)
    ans$outcome <- nm
    ans
  }), fill = TRUE)
  eff$cohort <- cohort
  eff$estimate_per_PM25_SD <- eff$estimate_per_ug_m3 * pm_sd
  eff$CI95_lower_per_PM25_SD <- eff$CI95_lower * pm_sd
  eff$CI95_upper_per_PM25_SD <- eff$CI95_upper * pm_sd
  eff$FDR <- p.adjust(eff$p_value, "BH")
  eff$analysis_method <- "linear model on published covariate-residualized expression; subject-cluster robust SE"
  all_effects[[cohort]] <- eff
}

effects <- rbindlist(all_effects, fill = TRUE)
scores <- rbindlist(all_scores, fill = TRUE)
fwrite(effects, file.path(out, "pm25_signature_gene_effects.csv"))
fwrite(scores, file.path(out, "pm25_scores.csv"))
fwrite(rbindlist(coverage, fill = TRUE), file.path(out, "signature_coverage.csv"))

primary <- effects[outcome %in% c("DNA_DAMAGE_RESPONSE", "OXIDATIVE_STRESS", "ABNORMAL_REPAIR",
                                  "KRT14_KRT17_REPAIR", "CILIA_CONSENSUS", "MULTICILIOGENESIS",
                                  "MATURE_CILIATED", "Repair_Cilia_Imbalance")]
png(file.path(out, "figures/PM25_signature_effects.png"), width = 2200, height = 1300, res = 240)
par(mar = c(5, 15, 3, 2))
primary <- primary[order(primary$cohort, primary$outcome), ]
y <- seq_len(nrow(primary))
xr <- range(c(primary$CI95_lower_per_PM25_SD, primary$CI95_upper_per_PM25_SD), finite = TRUE)
plot(primary$estimate_per_PM25_SD, y, xlim = xr, ylim = c(0.5, length(y) + 0.5),
     pch = ifelse(primary$FDR < 0.05, 19, 1), yaxt = "n",
     xlab = "Effect per SD higher 30-day PM2.5 (95% CI)", ylab = "",
     main = "GSE108134: PM2.5 signature effects")
segments(primary$CI95_lower_per_PM25_SD, y, primary$CI95_upper_per_PM25_SD, y, lwd = 2)
axis(2, at = y, labels = paste(primary$cohort, primary$outcome, sep = " | "), las = 2, cex.axis = 0.75)
abline(v = 0, lty = 2, col = "grey40")
dev.off()

report_lines <- c(
  "# GSE108134 PM2.5 unified-signature analysis", "",
  sprintf("Public processed residual matrices: %d smoker samples and %d nonsmoker samples.",
          sum(scores$cohort == "smoker"), sum(scores$cohort == "nonsmoker")),
  sprintf("Unique subjects represented: %d smokers and %d nonsmokers.",
          unique(effects$n_subjects[effects$cohort == "smoker"])[1],
          unique(effects$n_subjects[effects$cohort == "nonsmoker"])[1]), "",
  "Expression values are the investigators' covariate-residualized matrices. Each prespecified signature is the mean of gene-wise z scores. Associations use 30-day mean PM2.5 and subject-cluster robust standard errors.", "",
  "## Primary effects", ""
)
for (i in seq_len(nrow(primary))) {
  z <- primary[i]
  report_lines <- c(report_lines, sprintf("- %s | %s: effect/SD %.3f (95%% CI %.3f to %.3f), P=%s, FDR=%s.",
    z$cohort, z$outcome, z$estimate_per_PM25_SD, z$CI95_lower_per_PM25_SD,
    z$CI95_upper_per_PM25_SD, format(z$p_value, digits = 3), format(z$FDR, digits = 3)))
}
report_lines <- c(report_lines, "", "Associations are observational and do not establish a causal PM2.5-to-repair-to-cancer pathway.")
writeLines(report_lines, file.path(out, "analysis_summary.md"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
