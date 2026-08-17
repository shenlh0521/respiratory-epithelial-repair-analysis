#!/usr/bin/env Rscript
project <- Sys.getenv("PROJECT_ROOT", unset = normalizePath(".", mustWork = FALSE))
legacy_lib <- file.path(project, "files-mentioned-by-the-user-users/outputs/lung_bulla_cilia_cancer/work/Rlib")
.libPaths(c(legacy_lib, .libPaths()))
suppressPackageStartupMessages({library(data.table); library(AnnotationDbi); library(org.Hs.eg.db)})

root <- file.path(project, "04_smoke/GSE135188")
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
target_genes <- unique(c(unlist(gene_sets), "FOXJ1", "MCIDAS", "GMNC", paste0("RFX", 1:7), "KRT14", "KRT17"))

parse_series <- function(path) {
  x <- readLines(gzfile(path), warn = FALSE)
  fields <- function(prefix) {
    lapply(x[startsWith(x, prefix)], function(z) gsub('^"|"$', '', strsplit(z, "\t", fixed = TRUE)[[1]][-1]))
  }
  title <- fields("!Sample_title")[[1]]; gsm <- fields("!Sample_geo_accession")[[1]]
  ch <- fields("!Sample_characteristics_ch1")
  get_char <- function(label) {
    hit <- which(vapply(ch, function(v) startsWith(tolower(v[1]), tolower(label)), logical(1)))[1]
    if (is.na(hit)) return(rep(NA_character_, length(title)))
    sub(paste0("^", label, "\\s*"), "", ch[[hit]], ignore.case = TRUE)
  }
  data.table(sample_id = title, gsm = gsm, donor = get_char("donor:"),
             disease = get_char("disease state:"), treatment = get_char("treatment:"),
             age = as.numeric(get_char("age:")), sex = get_char("sex:"))
}

meta <- parse_series(file.path(root, "data/GSE135188_series_matrix.txt.gz"))
expr <- fread(cmd = sprintf("gzip -cd %s", shQuote(file.path(root, "data/GSE135188_uniqueRPKM.txt.gz"))),
              check.names = FALSE)
setnames(expr, 1, "ensembl_id")
mapping <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = expr$ensembl_id,
                                 columns = "SYMBOL", keytype = "ENSEMBL")
mapping <- mapping[!is.na(mapping$SYMBOL) & nzchar(mapping$SYMBOL), ]
mapping$SYMBOL <- toupper(mapping$SYMBOL)
expr[, SYMBOL := mapping$SYMBOL[match(ensembl_id, mapping$ENSEMBL)]]
expr <- expr[!is.na(SYMBOL) & SYMBOL %in% target_genes]
sample_cols <- intersect(meta$sample_id, names(expr))
expr[, mean_rpkm := rowMeans(.SD, na.rm = TRUE), .SDcols = sample_cols]
setorder(expr, SYMBOL, -mean_rpkm)
expr <- expr[!duplicated(SYMBOL)]
mat <- as.matrix(expr[, ..sample_cols]); storage.mode(mat) <- "double"; rownames(mat) <- expr$SYMBOL
mat <- log2(mat + 0.1)
z <- t(scale(t(mat))); z[!is.finite(z)] <- NA_real_
score <- function(genes) {
  use <- intersect(genes, rownames(z))
  if (!length(use)) return(rep(NA_real_, ncol(z)))
  colMeans(z[use, , drop = FALSE], na.rm = TRUE)
}
scores <- data.table(sample_id = colnames(z))
for (nm in names(gene_sets)) scores[[nm]] <- score(gene_sets[[nm]])
for (g in intersect(target_genes, rownames(z))) scores[[g]] <- z[g, ]
scores[, ABNORMAL_REPAIR_Z := as.numeric(scale(ABNORMAL_REPAIR))]
scores[, CILIA_CONSENSUS_Z := as.numeric(scale(CILIA_CONSENSUS))]
scores[, Repair_Cilia_Imbalance := ABNORMAL_REPAIR_Z - CILIA_CONSENSUS_Z]
scores <- merge(scores, meta, by = "sample_id", all.x = TRUE)
fwrite(scores, file.path(out, "smoke_well_scores.csv"))

coverage <- rbindlist(lapply(names(gene_sets), function(nm) data.frame(
  signature = nm, n_defined = length(gene_sets[[nm]]),
  n_detected = length(intersect(gene_sets[[nm]], rownames(z))),
  detected_genes = paste(intersect(gene_sets[[nm]], rownames(z)), collapse = ";"),
  missing_genes = paste(setdiff(gene_sets[[nm]], rownames(z)), collapse = ";"))))
fwrite(coverage, file.path(out, "signature_coverage.csv"))

outcomes <- unique(c(names(gene_sets), "Repair_Cilia_Imbalance", "FOXJ1", "MCIDAS", "GMNC",
                     intersect(paste0("RFX", 1:7), names(scores)), "KRT14", "KRT17"))
donor_scores <- scores[, lapply(.SD, mean, na.rm = TRUE), by = .(donor, disease, treatment), .SDcols = outcomes]
fwrite(donor_scores, file.path(out, "donor_treatment_scores.csv"))

paired_effect <- function(data, outcome, stratum) {
  w <- dcast(data, donor + disease ~ treatment, value.var = outcome)
  w[, delta := CS - air]
  tt <- t.test(w$delta)
  data.frame(stratum = stratum, outcome = outcome, n_donors = nrow(w),
             estimate_CS_minus_air = mean(w$delta), CI95_lower = tt$conf.int[1],
             CI95_upper = tt$conf.int[2], paired_effect_dz = mean(w$delta) / sd(w$delta),
             statistic = unname(tt$statistic), p_value = tt$p.value,
             analysis_method = "donor-paired t-test after averaging three technical/ALI wells",
             stringsAsFactors = FALSE)
}

effect_list <- lapply(outcomes, function(nm) paired_effect(donor_scores, nm, "all_donors"))
for (dis in unique(donor_scores$disease)) {
  effect_list <- c(effect_list,
                   lapply(outcomes, function(nm) paired_effect(donor_scores[disease == dis], nm, dis)))
}
effects <- rbindlist(effect_list, fill = TRUE)
effects[, FDR := p.adjust(p_value, "BH"), by = stratum]
fwrite(effects, file.path(out, "smoke_signature_gene_effects.csv"))

primary <- effects[stratum == "all_donors" & outcome %in% c("DNA_DAMAGE_RESPONSE", "OXIDATIVE_STRESS",
  "ABNORMAL_REPAIR", "KRT14_KRT17_REPAIR", "CILIA_CONSENSUS", "MULTICILIOGENESIS",
  "MATURE_CILIATED", "Repair_Cilia_Imbalance")][order(outcome)]
png(file.path(out, "figures/smoke_signature_effects.png"), width = 1900, height = 1050, res = 230)
par(mar = c(5, 14, 3, 2)); y <- seq_len(nrow(primary))
plot(primary$estimate_CS_minus_air, y, xlim = range(c(primary$CI95_lower, primary$CI95_upper)),
     yaxt = "n", pch = ifelse(primary$FDR < .05, 19, 1),
     xlab = "Whole cigarette smoke minus air (95% CI)", ylab = "",
     main = "GSE135188 donor-paired SAEC ALI smoke effects")
segments(primary$CI95_lower, y, primary$CI95_upper, y, lwd = 2)
axis(2, y, primary$outcome, las = 2); abline(v = 0, lty = 2, col = "grey40"); dev.off()

txt <- c("# GSE135188 whole-cigarette-smoke validation", "",
         "Primary human small-airway epithelial ALI cultures from six donors (three healthy, three COPD) were exposed to whole cigarette smoke or air. Three wells per donor-treatment were averaged before donor-paired inference.",
         "Canonical project signatures are unchanged from the PM2.5 analysis.", "", "## Primary effects", "")
for (i in seq_len(nrow(primary))) {
  q <- primary[i]
  txt <- c(txt, sprintf("- %s: smoke-air %.3f (95%% CI %.3f to %.3f), dz=%.3f, P=%s, FDR=%s.",
    q$outcome, q$estimate_CS_minus_air, q$CI95_lower, q$CI95_upper, q$paired_effect_dz,
    format(q$p_value, digits = 3), format(q$FDR, digits = 3)))
}
txt <- c(txt, "", "Disease-stratified three-donor results are sensitivity analyses only. The experiment supports exposure-response interpretation in vitro but does not by itself establish a lung-cancer mechanism.")
writeLines(txt, file.path(out, "analysis_summary.md"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
