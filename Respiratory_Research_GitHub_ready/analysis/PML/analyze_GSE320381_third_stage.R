#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(yaml)
})

set.seed(20260811)
project <- Sys.getenv("PROJECT_ROOT", unset = normalizePath(".", mustWork = FALSE))
legacy <- file.path(project, "files-mentioned-by-the-user-users/outputs/lung_bulla_cilia_cancer")
base <- file.path(project, "01_PML_extended/GSE320381")
raw <- file.path(base, "raw")
out <- file.path(base, "results")
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)

parse_soft <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  starts <- which(grepl("^\\^SAMPLE = ", lines)); ends <- c(starts[-1] - 1L, length(lines))
  rbindlist(lapply(seq_along(starts), function(i) {
    z <- lines[starts[i]:ends[i]]
    one <- function(prefix) {
      hit <- z[grepl(prefix, z)]; if (!length(hit)) return(NA_character_)
      sub(prefix, "", hit[1])
    }
    chars <- sub("^!Sample_characteristics_ch1 = ", "", z[grepl("^!Sample_characteristics_ch1 = ", z)])
    keys <- tolower(gsub("[^a-z0-9]+", "_", sub(":.*$", "", chars)))
    vals <- setNames(sub("^[^:]+: ", "", chars), keys)
    data.frame(gsm = sub("^\\^SAMPLE = ", "", z[1]),
               sample_id = one("^!Sample_title = "), tissue = vals["tissue"],
               patient_id = vals["patient_id"], histologic_grade = vals["histologic_grade"],
               simplified_grade = vals["simplified_histologic_grade"],
               dysplasia = as.integer(vals["dysplasia"]), platform = one("^!Sample_platform_id = "),
               stringsAsFactors = FALSE)
  }), fill = TRUE)
}

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
map <- fread(file.path(legacy, "data/metadata/target_symbol_ensembl.csv"), data.table = FALSE)
map <- rbind(map, data.frame(symbol = c("KRT4", "KRT13"),
                             ensembl_id = c("ENSG00000170477", "ENSG00000171401")))
map <- unique(map)

make_scores <- function(count_path, md) {
  x <- fread(cmd = sprintf("gzip -cd %s", shQuote(count_path)), data.table = FALSE, check.names = FALSE)
  names(x)[1] <- "ensembl_id"
  x$ensembl_id <- sub("\\.[0-9]+$", "", x$ensembl_id)
  sample_cols <- intersect(md$sample_id, names(x))
  stopifnot(length(sample_cols) == nrow(md))
  mat <- as.matrix(x[, sample_cols, drop = FALSE]); storage.mode(mat) <- "numeric"; rownames(mat) <- x$ensembl_id
  keep <- rowSums(mat >= 5) >= max(3, floor(0.05 * ncol(mat)))
  dge <- calcNormFactors(DGEList(mat[keep, , drop = FALSE]))
  logcpm <- cpm(dge, log = TRUE, prior.count = 1)
  m <- map[map$ensembl_id %in% rownames(logcpm), , drop = FALSE]
  target <- logcpm[m$ensembl_id, , drop = FALSE]; rownames(target) <- toupper(m$symbol)
  if (anyDuplicated(rownames(target))) target <- target[!duplicated(rownames(target)), , drop = FALSE]
  gz <- t(scale(t(target))); gz[!is.finite(gz)] <- NA_real_
  sc <- sapply(gene_sets, function(gs) {
    hit <- intersect(toupper(gs), rownames(gz))
    if (length(hit) < 2) rep(NA_real_, ncol(gz)) else colMeans(gz[hit, , drop = FALSE], na.rm = TRUE)
  })
  sc <- as.data.frame(sc, check.names = FALSE); sc$sample_id <- rownames(sc)
  for (g in c("FOXJ1", "KRT14", "KRT17", "MCIDAS", "GMNC", "KRT4", "KRT13"))
    sc[[g]] <- if (g %in% rownames(gz)) gz[g, sc$sample_id] else NA_real_
  coverage <- rbindlist(lapply(names(gene_sets), function(nm) {
    req <- toupper(gene_sets[[nm]]); hit <- intersect(req, rownames(gz))
    data.frame(signature = nm, requested_n = length(req), observed_n = length(hit),
               observed_genes = paste(hit, collapse = ";"), missing_genes = paste(setdiff(req, hit), collapse = ";"))
  }))
  list(scores = sc, coverage = coverage)
}

md <- parse_soft(file.path(raw, "GSE320381_family.soft.gz"))
grade_map <- c("NAD" = 0, "Normal" = 0, "Inflammation" = 0, "Reactive Airway Disease" = 0,
               "Basal Cell Hyperplasia" = 1, "Metaplasia" = 2, "Early Squamous Metaplasia" = 2,
               "Squamous Metaplasia" = 2, "Mild Dysplasia" = 3, "Moderate Dysplasia" = 4,
               "Severe Dysplasia" = 5, "Carcinoma In Situ" = 6, "Invasive Carcinoma" = 7)
md$grade_ordinal <- unname(grade_map[md$histologic_grade])
bxmd <- md[md$tissue == "Bronchial Biopsy", , drop = FALSE]
brmd <- md[md$tissue == "Bronchial Brush", , drop = FALSE]
bxobj <- make_scores(file.path(raw, "GSE320381_BX_counts.tsv.gz"), bxmd)
brobj <- make_scores(file.path(raw, "GSE320381_BR_counts.tsv.gz"), brmd)
bx <- merge(bxmd, bxobj$scores, by = "sample_id", all.x = TRUE, sort = FALSE)
br <- merge(brmd, brobj$scores, by = "sample_id", all.x = TRUE, sort = FALSE)
bx <- as.data.frame(bx); br <- as.data.frame(br)

zscore <- function(x) as.numeric(scale(x))
rank_normal <- function(x) { ans <- rep(NA_real_, length(x)); ok <- is.finite(x); n <- sum(ok); ans[ok] <- qnorm((rank(x[ok]) - .5) / n); ans }
for (dname in c("bx", "br")) {
  d <- get(dname)
  d$Repair_Cilia_Imbalance <- zscore(d$ABNORMAL_REPAIR) - zscore(d$CILIA_CONSENSUS)
  d$Repair_Cilia_Imbalance_rank <- rank_normal(d$ABNORMAL_REPAIR) - rank_normal(d$CILIA_CONSENSUS)
  assign(dname, d)
}
analysis_vars <- c(names(gene_sets), "Repair_Cilia_Imbalance", "Repair_Cilia_Imbalance_rank",
                   "FOXJ1", "KRT14", "KRT17", "MCIDAS", "GMNC", "KRT4", "KRT13")

cluster_fit <- function(d, variable, rhs, term, contrast, cluster = "patient_id") {
  vars <- unique(c(variable, all.vars(as.formula(paste("~", rhs))), cluster))
  z <- d[complete.cases(d[, vars, drop = FALSE]), , drop = FALSE]
  if (nrow(z) < 12 || length(unique(z[[cluster]])) < 5 || !is.finite(sd(z[[variable]])) || sd(z[[variable]]) == 0)
    return(data.frame(variable = variable, contrast = contrast, n_samples = nrow(z), n_subjects = length(unique(z[[cluster]])),
                      estimate_SD = NA_real_, std_error = NA_real_, CI95_lower = NA_real_, CI95_upper = NA_real_,
                      statistic = NA_real_, p_value = NA_real_, inference = "not estimable", stringsAsFactors = FALSE))
  y <- zscore(z[[variable]]); X0 <- model.matrix(as.formula(paste("~", rhs)), z)
  q <- qr(X0); keep <- sort(q$pivot[seq_len(q$rank)]); X <- X0[, keep, drop = FALSE]
  if (!term %in% colnames(X)) return(data.frame(variable = variable, contrast = contrast, n_samples = nrow(z),
    n_subjects = length(unique(z[[cluster]])), estimate_SD = NA_real_, std_error = NA_real_, CI95_lower = NA_real_,
    CI95_upper = NA_real_, statistic = NA_real_, p_value = NA_real_, inference = "term aliased", stringsAsFactors = FALSE))
  beta <- solve(crossprod(X), crossprod(X, y)); e <- as.numeric(y - X %*% beta); bread <- solve(crossprod(X))
  groups <- split(seq_len(nrow(z)), factor(z[[cluster]]), drop = TRUE); meat <- matrix(0, ncol(X), ncol(X))
  for (idx in groups) { u <- as.numeric(crossprod(X[idx, , drop = FALSE], e[idx])); meat <- meat + tcrossprod(u) }
  G <- length(groups); N <- nrow(z); P <- ncol(X); vc <- (G/(G-1))*((N-1)/(N-P))*bread%*%meat%*%bread
  j <- match(term, colnames(X)); est <- beta[j]; se <- sqrt(vc[j,j]); crit <- qt(.975, G-1); stat <- est/se
  data.frame(variable = variable, contrast = contrast, n_samples = N, n_subjects = G, estimate_SD = est,
             std_error = se, CI95_lower = est-crit*se, CI95_upper = est+crit*se, statistic = stat,
             p_value = 2*pt(abs(stat), G-1, lower.tail = FALSE), inference = "patient-cluster robust", stringsAsFactors = FALSE)
}

severity <- rbindlist(lapply(analysis_vars, function(v) rbind(
  cluster_fit(bx, v, "grade_ordinal", "grade_ordinal", "per_one_histology_grade_primary"),
  cluster_fit(bx, v, "grade_ordinal + patient_id", "grade_ordinal", "within_patient_per_grade_sensitivity")
)))
severity[, FDR := p.adjust(p_value, "BH"), by = contrast]
fwrite(severity, file.path(out, "biopsy_histology_trends.csv"))

dysplasia_effects <- rbindlist(lapply(analysis_vars, function(v)
  cluster_fit(bx, v, "dysplasia", "dysplasia", "dysplasia_CIS_SCC_vs_normal_hyperplasia_metaplasia")))
dysplasia_effects[, FDR := p.adjust(p_value, "BH")]
fwrite(dysplasia_effects, file.path(out, "biopsy_dysplasia_effects.csv"))

maxgrade <- aggregate(grade_ordinal ~ patient_id, bx, max, na.rm = TRUE)
names(maxgrade)[2] <- "max_biopsy_grade"
brf <- merge(br, maxgrade, by = "patient_id", all.x = TRUE)
field <- rbindlist(lapply(analysis_vars, function(v)
  cluster_fit(brf, v, "max_biopsy_grade + platform", "max_biopsy_grade", "brushing_per_max_biopsy_grade")))
field[, FDR := p.adjust(p_value, "BH")]
fwrite(field, file.path(out, "brushing_field_injury_trends.csv"))

corrows <- list()
for (v in analysis_vars) {
  bxp <- aggregate(bx[[v]], list(patient_id = bx$patient_id), mean, na.rm = TRUE); names(bxp)[2] <- "biopsy_mean"
  brp <- aggregate(br[[v]], list(patient_id = br$patient_id), mean, na.rm = TRUE); names(brp)[2] <- "brush_mean"
  z <- merge(bxp, brp, by = "patient_id"); z <- z[complete.cases(z), ]
  for (method in c("pearson", "spearman")) {
    ct <- suppressWarnings(cor.test(z$biopsy_mean, z$brush_mean, method = method, exact = FALSE))
    ci <- if (method == "pearson") ct$conf.int else c(NA_real_, NA_real_)
    corrows[[length(corrows)+1]] <- data.frame(variable=v, method=method, n_subjects=nrow(z), correlation=unname(ct$estimate),
      CI95_lower=ci[1], CI95_upper=ci[2], p_value=ct$p.value, comparison="patient-mean biopsy vs patient-mean brushing")
  }
}
cors <- rbindlist(corrows); cors[, FDR := p.adjust(p_value, "BH")]
fwrite(cors, file.path(out, "biopsy_brushing_correlations.csv"))

bxobj$coverage[, tissue := "biopsy"]; brobj$coverage[, tissue := "brushing"]
fwrite(rbind(bxobj$coverage, brobj$coverage), file.path(out, "signature_coverage.csv"))
fwrite(bx, file.path(out, "biopsy_metadata_scores.csv")); fwrite(brf, file.path(out, "brushing_metadata_scores.csv"))

plotdat <- severity[contrast == "per_one_histology_grade_primary" & variable %in% c("ABNORMAL_REPAIR", "KRT14_KRT17_REPAIR",
  "CILIA_CONSENSUS", "MULTICILIOGENESIS", "MATURE_CILIATED", "Repair_Cilia_Imbalance", "KRT4_KRT13_TRANSITIONAL_FATE", "FOXJ1")]
png(file.path(out, "figures", "GSE320381_severity_effects.png"), width=2200, height=1450, res=260)
par(mar=c(5,12,3,2)); plotdat <- plotdat[order(estimate_SD)]; y <- seq_len(nrow(plotdat)); xr <- range(c(plotdat$CI95_lower,plotdat$CI95_upper),finite=TRUE)
plot(plotdat$estimate_SD,y,xlim=xr,ylim=c(.5,length(y)+.5),yaxt="n",pch=19,xlab="SD change per histology grade (95% CI)",ylab="",
     main="GSE320381 biopsy severity: patient-cluster robust")
segments(plotdat$CI95_lower,y,plotdat$CI95_upper,y,lwd=2); axis(2,y,plotdat$variable,las=2,cex.axis=.8); abline(v=0,lty=2,col="grey40"); dev.off()

imb <- severity[contrast == "per_one_histology_grade_primary" & variable == "Repair_Cilia_Imbalance"]
fld <- field[variable == "Repair_Cilia_Imbalance"]
writeLines(c("# GSE320381 third-stage analysis", "",
  sprintf("Biopsies: %d from %d subjects; brushings: %d from %d subjects.", nrow(bx), length(unique(bx$patient_id)), nrow(br), length(unique(br$patient_id))),
  sprintf("Biopsy Repair_Cilia_Imbalance per histology grade = %.3f SD (95%% CI %.3f to %.3f), P=%.3g, FDR=%.3g.", imb$estimate_SD,imb$CI95_lower,imb$CI95_upper,imb$p_value,imb$FDR),
  sprintf("Brushing imbalance per one-grade increase in the subject's maximum biopsy histology = %.3f SD (95%% CI %.3f to %.3f), P=%.3g, FDR=%.3g.", fld$estimate_SD,fld$CI95_lower,fld$CI95_upper,fld$p_value,fld$FDR), "",
  "Primary severity inference uses patient-cluster robust standard errors. All biopsies were sequenced on GPL30173, so no biopsy platform term was needed; a within-patient fixed-effect sensitivity model is reported separately.",
  "The brushing analysis tests a field-of-injury association; it is cross-sectional and does not imply that brushing expression causes biopsy severity."
), file.path(out, "analysis_summary.md"))
capture.output(sessionInfo(), file=file.path(out,"sessionInfo.txt"))
cat(sprintf("Completed GSE320381: %d biopsies, %d brushings, %d subjects total\n", nrow(bx),nrow(br),length(unique(md$patient_id))))
