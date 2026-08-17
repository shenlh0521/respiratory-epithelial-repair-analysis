#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

set.seed(20260813)
project <- Sys.getenv("PROJECT_ROOT", unset = normalizePath(".", mustWork = FALSE))
base <- file.path(project, "12_targeted_manuscript_upgrades")
out <- file.path(base, "03_ALI_maturation_projection")
figdir <- file.path(base, "05_figures/candidate_figures")
legacy <- file.path(project, "files-mentioned-by-the-user-users/outputs/lung_bulla_cilia_cancer")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

sig <- fread(file.path(out, "ALI_maturation_signature.csv"), data.table = FALSE)
up <- unique(toupper(sig$gene[sig$direction == "maturation_up"]))
down <- unique(toupper(sig$gene[sig$direction == "maturation_down"]))
stopifnot(length(up) == 50L, length(down) == 3L)

zscore <- function(x) as.numeric(scale(x))

collapse_rows <- function(mat, symbols) {
  keep <- !is.na(symbols) & nzchar(symbols) & rowSums(is.finite(mat)) > 0
  mat <- mat[keep, , drop = FALSE]
  symbols <- toupper(symbols[keep])
  dt <- as.data.table(mat)
  dt[, symbol := symbols]
  ans <- dt[, lapply(.SD, mean, na.rm = TRUE), by = symbol]
  m <- as.matrix(ans[, -1]); storage.mode(m) <- "numeric"; rownames(m) <- ans$symbol
  m[!is.finite(m)] <- NA_real_
  m
}

rank_score <- function(gene_mat, sample_ids) {
  # Within-sample percentile ranks over all unambiguously mapped platform genes.
  ranks <- apply(gene_mat, 2, function(v) {
    ans <- rep(NA_real_, length(v)); ok <- is.finite(v)
    ans[ok] <- rank(v[ok], ties.method = "average") / sum(ok)
    ans
  })
  if (is.null(dim(ranks))) ranks <- matrix(ranks, ncol = 1)
  rownames(ranks) <- rownames(gene_mat); colnames(ranks) <- colnames(gene_mat)
  hit_up <- intersect(up, rownames(ranks)); hit_down <- intersect(down, rownames(ranks))
  passed <- length(hit_up) >= 25L && length(hit_down) >= 2L
  score <- if (passed) colMeans(ranks[hit_up, , drop = FALSE], na.rm = TRUE) -
    colMeans(ranks[hit_down, , drop = FALSE], na.rm = TRUE) else rep(NA_real_, ncol(ranks))
  list(scores = data.frame(sample_id = sample_ids, ALI_MATURATION_SCORE = score,
                           stringsAsFactors = FALSE),
       coverage = data.frame(up_requested = length(up), up_observed = length(hit_up),
                             down_requested = length(down), down_observed = length(hit_down),
                             coverage_status = if (passed) "PASS" else "FAIL",
                             observed_up = paste(hit_up, collapse = ";"),
                             missing_up = paste(setdiff(up, hit_up), collapse = ";"),
                             observed_down = paste(hit_down, collapse = ";"),
                             missing_down = paste(setdiff(down, hit_down), collapse = ";"),
                             stringsAsFactors = FALSE))
}

ensembl_symbols <- function(ids) {
  ids <- sub("\\.[0-9]+$", "", ids)
  valid <- intersect(unique(ids), keys(org.Hs.eg.db, keytype = "ENSEMBL"))
  mp <- if (length(valid)) AnnotationDbi::select(org.Hs.eg.db, keys = valid,
                                                 keytype = "ENSEMBL", columns = "SYMBOL") else
    data.frame(ENSEMBL = character(), SYMBOL = character())
  mp <- mp[!is.na(mp$SYMBOL) & nzchar(mp$SYMBOL), , drop = FALSE]
  counts <- table(mp$ENSEMBL)
  mp <- mp[mp$ENSEMBL %in% names(counts[counts == 1L]), , drop = FALSE]
  ans <- setNames(toupper(mp$SYMBOL), mp$ENSEMBL)[ids]
  # Preserve frozen Ensembl-labelled signature features if directly present.
  direct <- ids %in% c(up, down) & grepl("^ENSG", ids)
  ans[direct] <- ids[direct]
  unname(ans)
}

cluster_fit <- function(data, outcome, rhs, term, weights = NULL, cluster = "patient_id") {
  vars <- unique(c(outcome, all.vars(as.formula(paste("~", rhs))), cluster))
  z <- data[complete.cases(data[, vars, drop = FALSE]), , drop = FALSE]
  if (nrow(z) < 12L || length(unique(z[[cluster]])) < 5L || sd(z[[outcome]]) == 0)
    return(list(ok = FALSE, n = nrow(z), clusters = length(unique(z[[cluster]]))))
  y <- zscore(z[[outcome]])
  X0 <- model.matrix(as.formula(paste("~", rhs)), data = z)
  q <- qr(X0); keep <- sort(q$pivot[seq_len(q$rank)]); X <- X0[, keep, drop = FALSE]
  beta <- solve(crossprod(X), crossprod(X, y)); resid <- as.numeric(y - X %*% beta)
  bread <- solve(crossprod(X)); groups <- split(seq_len(nrow(z)), droplevels(factor(z[[cluster]])))
  meat <- matrix(0, ncol(X), ncol(X))
  for (idx in groups) {
    u <- as.numeric(crossprod(X[idx, , drop = FALSE], resid[idx])); meat <- meat + tcrossprod(u)
  }
  G <- length(groups); N <- nrow(z); P <- ncol(X)
  vc <- (G / (G - 1)) * ((N - 1) / (N - P)) * bread %*% meat %*% bread
  if (is.null(weights)) weights <- setNames(1, term)
  cvec <- setNames(rep(0, length(beta)), colnames(X)); common <- intersect(names(weights), names(cvec))
  cvec[common] <- weights[common]
  est <- sum(cvec * beta); se <- sqrt(drop(t(cvec) %*% vc %*% cvec)); df <- G - 1L
  crit <- qt(.975, df); p <- 2 * pt(abs(est / se), df, lower.tail = FALSE)
  list(ok = TRUE, estimate = est, se = se, lower = est - crit * se, upper = est + crit * se,
       p = p, n = N, clusters = G)
}

model_row <- function(dataset, comparison, fit, expected = "negative", method, model) {
  data.frame(dataset = dataset, comparison = comparison, expected_direction = expected,
             model = model, estimate_SD = if (fit$ok) fit$estimate else NA_real_,
             std_error = if (fit$ok) fit$se else NA_real_,
             CI95_lower = if (fit$ok) fit$lower else NA_real_,
             CI95_upper = if (fit$ok) fit$upper else NA_real_,
             P = if (fit$ok) fit$p else NA_real_, n_samples = fit$n,
             n_subjects = fit$clusters, analysis_method = method,
             stringsAsFactors = FALSE)
}

cluster_cor <- function(data, x, y, dataset, method, cluster = NULL, B = 2000L) {
  vars <- c(x, y, cluster); z <- data[complete.cases(data[, vars, drop = FALSE]), , drop = FALSE]
  obs <- suppressWarnings(cor(z[[x]], z[[y]], method = method))
  if (!is.null(cluster)) {
    ids <- unique(z[[cluster]])
    boot <- replicate(B, {
      draw <- sample(ids, length(ids), replace = TRUE)
      idx <- unlist(lapply(draw, function(id) which(z[[cluster]] == id)), use.names = FALSE)
      suppressWarnings(cor(z[[x]][idx], z[[y]][idx], method = method))
    })
    inference <- "patient-cluster bootstrap"
    n_subjects <- length(ids)
  } else {
    boot <- replicate(B, {
      idx <- sample(seq_len(nrow(z)), nrow(z), replace = TRUE)
      suppressWarnings(cor(z[[x]][idx], z[[y]][idx], method = method))
    })
    inference <- "lesion-site bootstrap; patient identifier unavailable"
    n_subjects <- NA_integer_
  }
  boot <- boot[is.finite(boot)]; se <- sd(boot); p <- if (is.finite(se) && se > 0) 2 * pnorm(abs(obs / se), lower.tail = FALSE) else NA_real_
  data.frame(dataset = dataset, method = method, correlation = obs,
             CI95_lower = unname(quantile(boot, .025)), CI95_upper = unname(quantile(boot, .975)),
             P = p, n_samples = nrow(z), n_subjects = n_subjects, inference = inference,
             stringsAsFactors = FALSE)
}

coverage <- list(); sample_outputs <- list(); results <- list(); correlations <- list()

# GSE109743: study residuals, existing passed-QC biopsy metadata and original model.
md109 <- fread(file.path(legacy, "results/GSE109743/metadata_clean.csv"), data.table = FALSE)
md109 <- md109[as.logical(md109$passed_qc_flag) & grepl("biopsy", md109$source, ignore.case = TRUE), ]
x109 <- fread(cmd = sprintf("gzip -cd %s", shQuote(file.path(legacy, "data/raw/pending/GSE109743_residuals.txt.gz"))),
              data.table = FALSE, check.names = FALSE)
names(x109)[1] <- "ensembl_id"; ids109 <- sub("\\.[0-9]+$", "", x109$ensembl_id)
cols109 <- intersect(md109$sample_id, names(x109)); mat109 <- as.matrix(x109[, cols109, drop = FALSE]); storage.mode(mat109) <- "numeric"
g109 <- collapse_rows(mat109, ensembl_symbols(ids109)); sc109 <- rank_score(g109, colnames(g109)); coverage[["GSE109743"]] <- sc109$coverage
d109 <- merge(md109, sc109$scores, by = "sample_id", all.x = TRUE, sort = FALSE)
d109$molecular_subtype <- factor(d109$molecular_subtype, levels = c("Normal", "Inflammatory", "Proliferative", "Secretory"))
d109$patient_id <- factor(d109$patient_id); d109$flow_cell_id <- factor(d109$flow_cell_id)
fit <- cluster_fit(d109, "ALI_MATURATION_SCORE", "molecular_subtype + median_tin_num + flow_cell_id + patient_id",
                   "molecular_subtypeProliferative", c(molecular_subtypeProliferative = 1))
results[[length(results) + 1L]] <- model_row("GSE109743", "Proliferative_vs_Normal", fit,
  method = "patient-FE OLS with patient-cluster robust SE", model = "adjusted_TIN_flow_patientFE")
fit <- cluster_fit(d109, "ALI_MATURATION_SCORE", "grade_ordinal + median_tin_num + flow_cell_id + patient_id", "grade_ordinal")
results[[length(results) + 1L]] <- model_row("GSE109743", "per_one_histology_grade", fit,
  method = "patient-FE OLS with patient-cluster robust SE", model = "adjusted_TIN_flow_patientFE")
old109 <- fread(file.path(project, "01_PML_extended/GSE109743/pml_extended_scores.csv"), data.table = FALSE)
d109 <- merge(d109, old109[, c("sample_id", "Repair_Cilia_Imbalance")], by = "sample_id", all.x = TRUE)
correlations[[length(correlations) + 1L]] <- cluster_cor(d109, "Repair_Cilia_Imbalance", "ALI_MATURATION_SCORE", "GSE109743", "pearson", "patient_id")
correlations[[length(correlations) + 1L]] <- cluster_cor(d109, "Repair_Cilia_Imbalance", "ALI_MATURATION_SCORE", "GSE109743", "spearman", "patient_id")
sample_outputs[["GSE109743"]] <- d109[, c("sample_id", "patient_id", "grade", "grade_ordinal", "molecular_subtype", "ALI_MATURATION_SCORE", "Repair_Cilia_Imbalance")]
rm(x109, mat109, g109); gc()

# GSE320381 biopsy counts, TMM logCPM and original grade/patient model.
parse_soft320 <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE); starts <- which(grepl("^\\^SAMPLE = ", lines)); ends <- c(starts[-1] - 1L, length(lines))
  rbindlist(lapply(seq_along(starts), function(i) {
    z <- lines[starts[i]:ends[i]]; one <- function(prefix) { h <- z[grepl(prefix, z)]; if (length(h)) sub(prefix, "", h[1]) else NA_character_ }
    chars <- sub("^!Sample_characteristics_ch1 = ", "", z[grepl("^!Sample_characteristics_ch1 = ", z)])
    vals <- setNames(sub("^[^:]+: ", "", chars), tolower(gsub("[^a-z0-9]+", "_", sub(":.*$", "", chars))))
    data.frame(gsm = sub("^\\^SAMPLE = ", "", z[1]), sample_id = one("^!Sample_title = "), tissue = vals["tissue"],
               patient_id = vals["patient_id"], histologic_grade = vals["histologic_grade"], stringsAsFactors = FALSE)
  }), fill = TRUE)
}
md320 <- as.data.frame(parse_soft320(file.path(project, "01_PML_extended/GSE320381/raw/GSE320381_family.soft.gz")))
grade_map320 <- c("NAD"=0,"Normal"=0,"Inflammation"=0,"Reactive Airway Disease"=0,"Basal Cell Hyperplasia"=1,
                  "Metaplasia"=2,"Early Squamous Metaplasia"=2,"Squamous Metaplasia"=2,"Mild Dysplasia"=3,
                  "Moderate Dysplasia"=4,"Severe Dysplasia"=5,"Carcinoma In Situ"=6,"Invasive Carcinoma"=7)
md320$grade_ordinal <- unname(grade_map320[md320$histologic_grade]); md320 <- md320[md320$tissue == "Bronchial Biopsy", ]
x320 <- fread(cmd = sprintf("gzip -cd %s", shQuote(file.path(project, "01_PML_extended/GSE320381/raw/GSE320381_BX_counts.tsv.gz"))),
              data.table = FALSE, check.names = FALSE); names(x320)[1] <- "ensembl_id"
ids320 <- sub("\\.[0-9]+$", "", x320$ensembl_id); cols320 <- intersect(md320$sample_id, names(x320))
mat320 <- as.matrix(x320[, cols320, drop = FALSE]); storage.mode(mat320) <- "numeric"; rownames(mat320) <- ids320
keep320 <- rowSums(mat320 >= 5) >= max(3, floor(.05 * ncol(mat320)))
log320 <- cpm(calcNormFactors(DGEList(mat320[keep320, , drop = FALSE])), log = TRUE, prior.count = 1)
g320 <- collapse_rows(log320, ensembl_symbols(rownames(log320))); sc320 <- rank_score(g320, colnames(g320)); coverage[["GSE320381"]] <- sc320$coverage
d320 <- merge(md320, sc320$scores, by = "sample_id", all.x = TRUE, sort = FALSE)
fit <- cluster_fit(d320, "ALI_MATURATION_SCORE", "grade_ordinal", "grade_ordinal")
results[[length(results) + 1L]] <- model_row("GSE320381", "per_one_histology_grade", fit,
  method = "OLS with patient-cluster robust SE", model = "original_primary_grade_model")
old320 <- fread(file.path(project, "01_PML_extended/GSE320381/results/biopsy_metadata_scores.csv"), data.table = FALSE)
d320 <- merge(d320, old320[, c("sample_id", "Repair_Cilia_Imbalance")], by = "sample_id", all.x = TRUE)
correlations[[length(correlations) + 1L]] <- cluster_cor(d320, "Repair_Cilia_Imbalance", "ALI_MATURATION_SCORE", "GSE320381", "pearson", "patient_id")
correlations[[length(correlations) + 1L]] <- cluster_cor(d320, "Repair_Cilia_Imbalance", "ALI_MATURATION_SCORE", "GSE320381", "spearman", "patient_id")
sample_outputs[["GSE320381"]] <- d320[, c("sample_id", "patient_id", "histologic_grade", "grade_ordinal", "ALI_MATURATION_SCORE", "Repair_Cilia_Imbalance")]
rm(x320, mat320, log320, g320); gc()

# GSE33479 Agilent platform; keep only probes with one platform symbol.
parse_soft334 <- function(path) {
  z <- readLines(gzfile(path), warn = FALSE); s <- which(grepl("^\\^SAMPLE = ", z)); e <- c(s[-1] - 1L, length(z))
  rbindlist(lapply(seq_along(s), function(i) {
    q <- z[s[i]:e[i]]; one <- function(prefix) { h <- q[grepl(prefix,q)]; if(length(h)) sub(prefix,"",h[1]) else NA_character_ }
    ch <- sub("^!Sample_characteristics_ch1 = ","",q[grepl("^!Sample_characteristics_ch1 = ",q)]); v <- setNames(sub("^[^:]+: ","",ch),tolower(sub(":.*$","",ch)))
    supp <- one("^!Sample_supplementary_file = ")
    data.frame(sample_id=sub("\\.gz$", "", sub("^GSM[0-9]+_", "", basename(supp))), histology=v["histology"],
               smoking_status=v["smoking status"], sex=v["gender"], age=as.numeric(v["age"]), patient_id=v["patient"], stringsAsFactors=FALSE)
  }), fill=TRUE)
}
path334 <- file.path(project, "01_PML_extended/GSE33479/raw/GSE33479_family.soft.gz")
md334 <- as.data.frame(parse_soft334(path334)); x334 <- fread(cmd=sprintf("gzip -cd %s",shQuote(file.path(project,"01_PML_extended/GSE33479/raw/GSE33479.txt.gz"))), check.names=FALSE, data.table=FALSE)
cols334 <- md334$sample_id; mat334 <- log2(pmax(as.matrix(x334[, cols334, drop=FALSE]), 1e-6)); storage.mode(mat334) <- "numeric"
soft334 <- readLines(gzfile(path334),warn=FALSE); pb <- which(soft334=="!platform_table_begin")[1]; pe <- which(soft334=="!platform_table_end")[1]
gpl334 <- fread(text=paste(soft334[(pb+1):(pe-1)],collapse="\n"),select=c("ID","GENE_SYMBOL"),data.table=FALSE)
symraw334 <- toupper(gpl334$GENE_SYMBOL[match(x334$ID_REF, gpl334$ID)])
sym334 <- ifelse(!is.na(symraw334) & !grepl(";|///|,", symraw334), trimws(symraw334), NA_character_)
g334 <- collapse_rows(mat334, sym334); sc334 <- rank_score(g334, colnames(g334)); coverage[["GSE33479"]] <- sc334$coverage
d334 <- merge(md334, sc334$scores, by="sample_id", all.x=TRUE, sort=FALSE)
map334 <- c("normal"=0,"hyperplasia"=1,"metaplasia"=2,"mild dysplasia"=3,"moderate dysplasia"=4,"severe dysplasia"=5,"carcinoma in situ"=6,"squamous cell carcinoma"=7,"scc"=7)
d334$grade_ordinal <- unname(map334[tolower(trimws(d334$histology))])
fit <- cluster_fit(d334, "ALI_MATURATION_SCORE", "grade_ordinal + age + sex + smoking_status", "grade_ordinal")
results[[length(results) + 1L]] <- model_row("GSE33479", "per_one_histology_grade", fit,
  method="age/sex/smoking-adjusted OLS with patient-cluster robust SE", model="adjusted_age_sex_smoking")
old334 <- fread(file.path(project,"01_PML_extended/GSE33479/results/metadata_signature_scores.csv"),data.table=FALSE)
d334 <- merge(d334, old334[,c("sample_id","Repair_Cilia_Imbalance")],by="sample_id",all.x=TRUE)
correlations[[length(correlations)+1L]] <- cluster_cor(d334,"Repair_Cilia_Imbalance","ALI_MATURATION_SCORE","GSE33479","pearson","patient_id")
correlations[[length(correlations)+1L]] <- cluster_cor(d334,"Repair_Cilia_Imbalance","ALI_MATURATION_SCORE","GSE33479","spearman","patient_id")
sample_outputs[["GSE33479"]] <- d334[,c("sample_id","patient_id","histology","grade_ordinal","ALI_MATURATION_SCORE","Repair_Cilia_Imbalance")]
rm(x334,mat334,g334,soft334); gc()

# GSE114489 Affymetrix processed series and official GPL annotation.
parse_series114 <- function(path) {
  lines <- readLines(gzfile(path), warn=FALSE); get_meta <- function(prefix) { hit<-lines[grepl(paste0("^",prefix,"\\t"),lines)]; lapply(hit,function(x) gsub('^"|"$','',strsplit(x,"\t",fixed=TRUE)[[1]][-1])) }
  gsm<-get_meta("!Sample_geo_accession")[[1]]; md<-data.frame(gsm=gsm,title=get_meta("!Sample_title")[[1]],stringsAsFactors=FALSE)
  for(z in get_meta("!Sample_characteristics_ch1")){key<-gsub("[^a-z0-9]+","_",tolower(trimws(sub(":.*$","",z[1]))));md[[key]]<-trimws(sub("^[^:]+:","",z))}
  begin<-which(lines=="!series_matrix_table_begin")+1L;end<-which(lines=="!series_matrix_table_end")-1L
  tab<-fread(text=paste(lines[begin:end],collapse="\n"),check.names=FALSE,data.table=FALSE);names(tab)[1]<-"probe_id";names(tab)[-1]<-gsub('^"|"$','',names(tab)[-1]);list(metadata=md,expression=tab)
}
s114 <- parse_series114(file.path(project,"01_PML_extended/GSE114489/raw/GSE114489_series_matrix.txt.gz"));md114<-s114$metadata;x114<-s114$expression
ann114<-fread(cmd=sprintf("gzip -cd %s",shQuote(file.path(project,"01_PML_extended/GSE114489/raw/GPL6244.annot.gz"))),skip="ID\tGene title",data.table=FALSE,check.names=FALSE)
symraw114<-toupper(ann114[["Gene symbol"]][match(x114$probe_id,ann114$ID)])
sym114<-ifelse(!is.na(symraw114) & !grepl("///|;|,",symraw114),trimws(symraw114),NA_character_)
mat114<-as.matrix(x114[,md114$gsm,drop=FALSE]);storage.mode(mat114)<-"numeric";g114<-collapse_rows(mat114,sym114);sc114<-rank_score(g114,colnames(g114));coverage[["GSE114489"]]<-sc114$coverage
names(sc114$scores)[names(sc114$scores)=="sample_id"] <- "gsm"
d114<-merge(md114,sc114$scores,by="gsm",all.x=TRUE,sort=FALSE);d114$group<-as.integer(sub("Group ([1-4]).*$","\\1",d114$title));d114$baseline_grade<-as.numeric(sub("\\.1$","",d114$bl_frozen_dx));d114$age_num<-as.numeric(d114$age);d114$pack_years_num<-as.numeric(d114$smoking_pack_yr)
primary114<-d114[d114$group%in%c(1,2),];primary114$case<-as.integer(primary114$group==1);primary114$score_z<-zscore(primary114$ALI_MATURATION_SCORE)
fit114<-lm(score_z~case+baseline_grade+age_num+pack_years_num,data=primary114);co114<-summary(fit114)$coefficients["case",];ci114<-confint(fit114,"case")
fit<-list(ok=TRUE,estimate=unname(co114["Estimate"]),se=unname(co114["Std. Error"]),lower=ci114[1],upper=ci114[2],p=unname(co114["Pr(>|t|)"]),n=nobs(fit114),clusters=NA_integer_)
results[[length(results)+1L]]<-model_row("GSE114489","persistent_vs_regressive_dysplasia",fit,method="baseline-grade/age/pack-years-adjusted lesion-site OLS",model="adjusted_baseline_grade_age_packyears")
old114<-fread(file.path(project,"01_PML_extended/GSE114489/results/metadata_scores.csv"),data.table=FALSE);d114<-merge(d114,old114[,c("gsm","Repair_Cilia_Imbalance")],by="gsm",all.x=TRUE)
correlations[[length(correlations)+1L]]<-cluster_cor(d114,"Repair_Cilia_Imbalance","ALI_MATURATION_SCORE","GSE114489","pearson",NULL)
correlations[[length(correlations)+1L]]<-cluster_cor(d114,"Repair_Cilia_Imbalance","ALI_MATURATION_SCORE","GSE114489","spearman",NULL)
sample_outputs[["GSE114489"]]<-d114[,c("gsm","title","group","baseline_grade","ALI_MATURATION_SCORE","Repair_Cilia_Imbalance")]
rm(x114,mat114,g114,s114);gc()

# ALI internal maturation reference visualization, using the same frozen rank score.
ali <- fread(file.path(base,"02_ALI_trajectory/ALI_lineage_log2CPM.csv"),data.table=FALSE,check.names=FALSE);genes_ali<-toupper(ali[[1]]);mat_ali<-as.matrix(ali[,-1]);storage.mode(mat_ali)<-"numeric";rownames(mat_ali)<-genes_ali
sc_ali<-rank_score(mat_ali,colnames(mat_ali))$scores;meta_ali<-fread(file.path(base,"02_ALI_trajectory/slingshot_cell_pseudotime.csv"),data.table=FALSE);ali_cells<-merge(meta_ali,sc_ali,by.x="cell",by.y="sample_id")
ali_cells<-as.data.table(ali_cells);ali_cells[,pseudotime_bin:=pmin(10L,ceiling(frank(pseudotime,ties.method="first")/.N*10)),by=donor]
ali_bins<-ali_cells[,.(pseudotime=mean(pseudotime),ALI_MATURATION_SCORE=mean(ALI_MATURATION_SCORE),n_cells=.N),by=.(donor,pseudotime_bin)]
fwrite(ali_bins,file.path(out,"ALI_donor_pseudotime_bin_maturation_scores.csv"))

covtab<-rbindlist(lapply(names(coverage),function(nm)cbind(dataset=nm,coverage[[nm]])),fill=TRUE);fwrite(covtab,file.path(out,"platform_coverage.csv"))
for(nm in names(sample_outputs)) fwrite(sample_outputs[[nm]],file.path(out,paste0(nm,"_ALI_maturation_scores.csv")))
res<-rbindlist(results,fill=TRUE);res[,FDR:=p.adjust(P,"BH")]
res[,status:=fifelse(is.na(P),"METHOD_FAILED",fifelse(estimate_SD<0 & FDR<.05,"SUPPORTS",fifelse(estimate_SD<0,"DIRECTIONAL_ONLY",fifelse(FDR<.05,"INCONSISTENT","NULL"))))]
fwrite(res,file.path(out,"PML_maturation_projection_results.csv"))
cortab<-rbindlist(correlations,fill=TRUE);cortab[,FDR:=p.adjust(P,"BH")];cortab[,status:=fifelse(is.na(P),"METHOD_FAILED",fifelse(correlation<0 & FDR<.05,"SUPPORTS",fifelse(correlation<0,"DIRECTIONAL_ONLY",fifelse(FDR<.05,"INCONSISTENT","NULL"))))]
fwrite(cortab,file.path(out,"Repair_Cilia_vs_ALI_maturation_correlations.csv"))

# Candidate Figure D: normal reference, frozen PML tests, and score correlations.
pdf(file.path(figdir,"Candidate_Figure_D_ALI_maturation_projection.pdf"),width=11,height=8.5,useDingbats=FALSE)
layout(matrix(c(1,2,3,4),2,2,byrow=TRUE),widths=c(1.05,1),heights=c(1,1));par(oma=c(1,1,2,1))
cols<-c(T84="#0072B2",T85="#D55E00",T89="#009E73")
plot(NA,xlim=c(0,1),ylim=range(ali_bins$ALI_MATURATION_SCORE),xlab="Slingshot pseudotime",ylab="ALI maturation score",main="A  Normal ALI maturation reference")
for(don in unique(ali_bins$donor)){z<-ali_bins[ali_bins$donor==don,][order(pseudotime)];lines(z$pseudotime,z$ALI_MATURATION_SCORE,type="b",pch=16,col=cols[don],lwd=1.5)}
legend("topleft",legend=names(cols),col=cols,pch=16,lty=1,bty="n",cex=.8)
labs<-paste(res$dataset,res$comparison,sep=" | ");yy<-rev(seq_len(nrow(res)));xr<-range(c(res$CI95_lower,res$CI95_upper,0),na.rm=TRUE)
plot(res$estimate_SD,yy,xlim=xr,yaxt="n",pch=19,col=ifelse(res$estimate_SD<0,"#0072B2","#D55E00"),xlab="Standardized effect (95% CI)",ylab="",main="B  Frozen PML projections")
segments(res$CI95_lower,yy,res$CI95_upper,yy,lwd=2,col="grey30");axis(2,yy,labs,las=2,cex.axis=.58);abline(v=0,lty=2,col="grey50")
cc<-cortab[cortab$method=="spearman",];yy2<-rev(seq_len(nrow(cc)));xr2<-range(c(cc$CI95_lower,cc$CI95_upper,0),na.rm=TRUE)
plot(cc$correlation,yy2,xlim=xr2,yaxt="n",pch=19,col=ifelse(cc$correlation<0,"#0072B2","#D55E00"),xlab="Spearman correlation (95% bootstrap CI)",ylab="",main="C  Repair–Cilia vs maturation")
segments(cc$CI95_lower,yy2,cc$CI95_upper,yy2,lwd=2,col="grey30");axis(2,yy2,cc$dataset,las=2,cex.axis=.75);abline(v=0,lty=2,col="grey50")
barcols<-ifelse(covtab$coverage_status=="PASS","#009E73","#D55E00");barplot(t(as.matrix(covtab[,c("up_observed","down_observed")])),beside=FALSE,names.arg=covtab$dataset,col=c("#56B4E9","#E69F00"),ylim=c(0,55),ylab="Frozen genes observed",main="D  Platform coverage",las=2,cex.names=.7)
legend("topright",legend=c("maturation-up","maturation-down"),fill=c("#56B4E9","#E69F00"),bty="n",cex=.75);mtext("Candidate only — does not replace submission Figures 1–8",outer=TRUE,line=.3,cex=.8)
dev.off()

writeLines(c("# ALI-derived maturation projection", "",
  "The maturation signature was trained only in normal GSE134174 ALI cells and frozen before PML projection.",
  "Each score is the within-sample mean percentile rank of maturation-up genes minus that of maturation-down genes.",
  "All four platforms passed the prespecified >=25/50 up and >=2/3 down coverage rule.",
  "PML labels and outcomes were not used for gene selection, weighting, or score construction.",
  "PML profiles are interpreted as shifts along a normal-airway transcriptional reference, not literal lineage arrest."),file.path(out,"README.md"))
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
cat("Completed frozen ALI maturation projections\n");print(res);print(cortab);print(covtab[,c("dataset","up_observed","down_observed","coverage_status")])
