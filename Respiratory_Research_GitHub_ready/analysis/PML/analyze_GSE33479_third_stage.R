#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(yaml)})
set.seed(20260811)
project <- Sys.getenv("PROJECT_ROOT", unset = normalizePath(".", mustWork = FALSE))
base <- file.path(project, "01_PML_extended/GSE33479")
out <- file.path(base, "results"); dir.create(file.path(out,"figures"), recursive=TRUE, showWarnings=FALSE)
legacy <- file.path(project, "files-mentioned-by-the-user-users/outputs/lung_bulla_cilia_cancer")

parse_soft <- function(path) {
  z <- readLines(gzfile(path), warn=FALSE)
  s <- which(grepl("^\\^SAMPLE = ", z)); e <- c(s[-1]-1L, length(z))
  rbindlist(lapply(seq_along(s), function(i) {
    q <- z[s[i]:e[i]]
    one <- function(prefix) {h <- q[grepl(prefix,q)]; if(length(h)) sub(prefix,"",h[1]) else NA_character_}
    ch <- sub("^!Sample_characteristics_ch1 = ","",q[grepl("^!Sample_characteristics_ch1 = ",q)])
    v <- setNames(sub("^[^:]+: ","",ch),tolower(sub(":.*$","",ch)))
    supp <- one("^!Sample_supplementary_file = ")
    data.frame(gsm=one("^!Sample_geo_accession = "), title=one("^!Sample_title = "),
      sample_id=sub("\\.gz$", "", sub("^GSM[0-9]+_", "", basename(supp))), histology=v["histology"],
      phenotypical_stage=v["phenotypical stage"], molecular_group=v["molecular group"],
      smoking_status=v["smoking status"], sex=v["gender"], age=as.numeric(v["age"]),
      patient_id=v["patient"], stringsAsFactors=FALSE)
  }), fill=TRUE)
}

md <- parse_soft(file.path(base,"raw/GSE33479_family.soft.gz"))
x <- fread(cmd=sprintf("gzip -cd %s",shQuote(file.path(base,"raw/GSE33479.txt.gz"))), check.names=FALSE, data.table=FALSE)
stopifnot(nrow(md)==122, all(md$sample_id %in% names(x)))
sample_cols <- md$sample_id
mat <- as.matrix(x[,sample_cols,drop=FALSE]); storage.mode(mat) <- "numeric"
# Values are positive normalized ratios; log2 stabilizes their right-skew before gene-wise z scoring.
mat <- log2(pmax(mat, 1e-6))
soft_lines <- readLines(gzfile(file.path(base,"raw/GSE33479_family.soft.gz")),warn=FALSE)
pb <- which(soft_lines=="!platform_table_begin")[1]; pe <- which(soft_lines=="!platform_table_end")[1]
gpl <- fread(text=paste(soft_lines[(pb+1):(pe-1)],collapse="\n"),select=c("ID","GENE_SYMBOL"),data.table=FALSE)
probe_symbol <- setNames(toupper(gpl$GENE_SYMBOL),gpl$ID)
aliases <- strsplit(ifelse(is.na(probe_symbol[x[["ID_REF"]]]),"",probe_symbol[x[["ID_REF"]]]),";[ ]*")

legacy_sets <- yaml.load_file(file.path(legacy,"config/gene_sets.yaml"))$gene_sets
gene_sets <- list(
  DNA_DAMAGE_RESPONSE=legacy_sets$DNA_DAMAGE_SENTINELS,
  ABNORMAL_REPAIR=legacy_sets$ABNORMAL_REPAIR,
  KRT14_KRT17_REPAIR=c("KRT14","KRT17"),
  CILIA_CONSENSUS=legacy_sets$CILIA_CONSENSUS,
  MULTICILIOGENESIS=unique(c(legacy_sets$CILIA_REGULATORY,legacy_sets$DEUTEROSOMAL)),
  MATURE_CILIATED=c("FOXJ1","PIFO","TPPP3","CAPS","RSPH1"),
  KRT4_KRT13_TRANSITIONAL_FATE=c("KRT4","KRT13")
)
targets <- unique(toupper(c(unlist(gene_sets),"FOXJ1","KRT14","KRT17","MCIDAS","GMNC","KRT4","KRT13")))
gene_expr <- lapply(targets, function(g) {
  idx <- which(vapply(aliases, function(a) g %in% trimws(a), logical(1)))
  if(!length(idx)) return(NULL)
  if(length(idx)==1) mat[idx,] else apply(mat[idx,,drop=FALSE],2,median,na.rm=TRUE)
})
names(gene_expr) <- targets; gene_expr <- gene_expr[!vapply(gene_expr,is.null,logical(1))]
gmat <- do.call(rbind,gene_expr); rownames(gmat) <- names(gene_expr)
gz <- t(scale(t(gmat))); gz[!is.finite(gz)] <- NA_real_
scores <- sapply(gene_sets, function(gs) {
  hit <- intersect(toupper(gs),rownames(gz)); if(length(hit)<2) rep(NA_real_,ncol(gz)) else colMeans(gz[hit,,drop=FALSE],na.rm=TRUE)
})
scores <- as.data.frame(scores,check.names=FALSE); scores$sample_id <- rownames(scores)
for(g in c("FOXJ1","KRT14","KRT17","MCIDAS","GMNC","KRT4","KRT13")) scores[[g]] <- if(g %in% rownames(gz)) gz[g,] else NA_real_
d <- merge(md,scores,by="sample_id",all.x=TRUE,sort=FALSE)
d <- as.data.frame(d)

norm_hist <- tolower(trimws(d$histology))
grade_map <- c("normal"=0,"hyperplasia"=1,"metaplasia"=2,"mild dysplasia"=3,
               "moderate dysplasia"=4,"severe dysplasia"=5,"carcinoma in situ"=6,
               "squamous cell carcinoma"=7,"scc"=7)
d$grade_ordinal <- unname(grade_map[norm_hist])
if(anyNA(d$grade_ordinal)) stop("Unmapped histology: ",paste(unique(d$histology[is.na(d$grade_ordinal)]),collapse=", "))
zscore <- function(v) as.numeric(scale(v))
rank_normal <- function(v) {o<-rep(NA_real_,length(v)); ok<-is.finite(v); o[ok]<-qnorm((rank(v[ok])-.5)/sum(ok));o}
d$Repair_Cilia_Imbalance <- zscore(d$ABNORMAL_REPAIR)-zscore(d$CILIA_CONSENSUS)
d$Repair_Cilia_Imbalance_rank <- rank_normal(d$ABNORMAL_REPAIR)-rank_normal(d$CILIA_CONSENSUS)
vars <- c(names(gene_sets),"Repair_Cilia_Imbalance","Repair_Cilia_Imbalance_rank","FOXJ1","KRT14","KRT17","MCIDAS","GMNC","KRT4","KRT13")

cluster_fit <- function(variable, adjusted=TRUE) {
  rhs <- if(adjusted) "grade_ordinal + age + sex + smoking_status" else "grade_ordinal"
  keepvars <- unique(c(variable,all.vars(as.formula(paste("~",rhs))),"patient_id"))
  z <- d[complete.cases(d[,keepvars,drop=FALSE]),,drop=FALSE]
  if(nrow(z)<12 || length(unique(z$patient_id))<5 || !is.finite(sd(z[[variable]])) || sd(z[[variable]])==0)
    return(data.frame(variable=variable,model=if(adjusted)"adjusted_age_sex_smoking" else "unadjusted",
      comparison="per_one_histology_grade",n_samples=nrow(z),n_subjects=length(unique(z$patient_id)),
      effect_size_SD=NA_real_,std_error=NA_real_,CI95_lower=NA_real_,CI95_upper=NA_real_,P=NA_real_,
      analysis_method="not estimable",stringsAsFactors=FALSE))
  z$sex <- droplevels(factor(z$sex)); z$smoking_status <- droplevels(factor(z$smoking_status))
  y <- zscore(z[[variable]]); X0 <- model.matrix(as.formula(paste("~",rhs)),z)
  q <- qr(X0); X <- X0[,sort(q$pivot[seq_len(q$rank)]),drop=FALSE]
  j <- match("grade_ordinal",colnames(X)); beta <- solve(crossprod(X),crossprod(X,y)); e <- as.numeric(y-X%*%beta)
  bread <- solve(crossprod(X)); groups <- split(seq_len(nrow(z)),factor(z$patient_id)); meat <- matrix(0,ncol(X),ncol(X))
  for(idx in groups){u<-as.numeric(crossprod(X[idx,,drop=FALSE],e[idx]));meat<-meat+tcrossprod(u)}
  G<-length(groups);N<-nrow(z);P<-ncol(X);vc<-(G/(G-1))*((N-1)/(N-P))*bread%*%meat%*%bread
  est<-beta[j];se<-sqrt(vc[j,j]);crit<-qt(.975,G-1);stat<-est/se
  data.frame(variable=variable,model=if(adjusted)"adjusted_age_sex_smoking" else "unadjusted",
    comparison="per_one_histology_grade",n_samples=N,n_subjects=G,effect_size_SD=est,std_error=se,
    CI95_lower=est-crit*se,CI95_upper=est+crit*se,P=2*pt(abs(stat),G-1,lower.tail=FALSE),
    analysis_method="OLS with patient-cluster robust SE",stringsAsFactors=FALSE)
}
res <- rbindlist(lapply(vars,function(v) rbind(cluster_fit(v,TRUE),cluster_fit(v,FALSE))))
res[,FDR:=p.adjust(P,"BH"),by=model]
fwrite(res,file.path(out,"histologic_continuum_trends.csv"))
fwrite(d,file.path(out,"metadata_signature_scores.csv"))

coverage <- rbindlist(lapply(names(gene_sets),function(nm){req<-toupper(gene_sets[[nm]]);hit<-intersect(req,rownames(gz));data.frame(signature=nm,requested_n=length(req),observed_n=length(hit),observed_genes=paste(hit,collapse=";"),missing_genes=paste(setdiff(req,hit),collapse=";"))}))
fwrite(coverage,file.path(out,"signature_coverage.csv"))
counts <- as.data.table(d)[,.(n_samples=.N,n_subjects=uniqueN(patient_id)),by=.(grade_ordinal,histology)][order(grade_ordinal)]
fwrite(counts,file.path(out,"histology_sample_counts.csv"))

p <- res[model=="adjusted_age_sex_smoking" & variable %in% c("ABNORMAL_REPAIR","KRT14_KRT17_REPAIR","CILIA_CONSENSUS","MULTICILIOGENESIS","MATURE_CILIATED","KRT4_KRT13_TRANSITIONAL_FATE","Repair_Cilia_Imbalance","FOXJ1")]
p <- p[order(effect_size_SD)]; png(file.path(out,"figures/GSE33479_histologic_continuum_effects.png"),2200,1450,res=260)
par(mar=c(5,12,3,2)); y<-seq_len(nrow(p)); xr<-range(c(p$CI95_lower,p$CI95_upper),finite=TRUE)
plot(p$effect_size_SD,y,xlim=xr,ylim=c(.5,length(y)+.5),yaxt="n",pch=19,xlab="SD change per histology grade (95% CI)",ylab="",main="GSE33479 full histologic continuum")
segments(p$CI95_lower,y,p$CI95_upper,y,lwd=2);axis(2,y,p$variable,las=2,cex.axis=.8);abline(v=0,lty=2,col="grey40");dev.off()

imb <- res[model=="adjusted_age_sex_smoking" & variable=="Repair_Cilia_Imbalance"]
writeLines(c("# GSE33479 third-stage analysis","",sprintf("Analyzed %d biopsies from %d participants across %d ordered histologic grades.",nrow(d),uniqueN(d$patient_id),uniqueN(d$grade_ordinal)),
 sprintf("Adjusted Repair_Cilia_Imbalance trend = %.3f SD per grade (95%% CI %.3f to %.3f), P=%.3g, FDR=%.3g.",imb$effect_size_SD,imb$CI95_lower,imb$CI95_upper,imb$P,imb$FDR),"",
 "The primary model adjusts age, sex, and smoking status and uses patient-cluster robust standard errors.",
 "This cross-sectional continuum supports association with histology severity but cannot establish temporal progression or causality."),file.path(out,"analysis_summary.md"))
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
print(counts);print(res[model=="adjusted_age_sex_smoking" & variable %in% c("ABNORMAL_REPAIR","CILIA_CONSENSUS","Repair_Cilia_Imbalance","KRT4_KRT13_TRANSITIONAL_FATE","FOXJ1")])
