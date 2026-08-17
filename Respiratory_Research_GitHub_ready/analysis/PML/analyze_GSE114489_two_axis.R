#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table)})
set.seed(20260811)
project <- Sys.getenv("PROJECT_ROOT", unset = normalizePath(".", mustWork = FALSE))
base <- file.path(project,"01_PML_extended/GSE114489")
raw <- file.path(base,"raw"); out <- file.path(base,"two_axis_results")
dir.create(file.path(out,"figures"),recursive=TRUE,showWarnings=FALSE)

parse_series_expr <- function(path) {
  lines <- readLines(gzfile(path),warn=FALSE); b<-which(lines=="!series_matrix_table_begin")+1L;e<-which(lines=="!series_matrix_table_end")-1L
  x<-fread(text=paste(lines[b:e],collapse="\n"),check.names=FALSE,data.table=FALSE);names(x)[1]<-"probe_id";names(x)[-1]<-gsub('^"|"$','',names(x)[-1]);x
}
read_annotation <- function(path){x<-fread(cmd=sprintf("gzip -cd %s",shQuote(path)),skip="ID\tGene title",data.table=FALSE,check.names=FALSE);a<-x[,c("ID","Gene symbol")];names(a)<-c("probe_id","symbol");a$symbol<-toupper(trimws(sub(" ///.*$","",a$symbol)));a[nzchar(a$symbol)&a$symbol!="---",]}
sig <- fread(file.path(project,"01_PML_extended/signatures/IMMUNE_MODULE9_ENRICHMENT_CORE.csv"),data.table=FALSE)$gene
expr <- merge(read_annotation(file.path(raw,"GPL6244.annot.gz")),parse_series_expr(file.path(raw,"GSE114489_series_matrix.txt.gz")),by="probe_id")
samples <- fread(file.path(base,"results/metadata_scores.csv"),data.table=FALSE)$gsm
dt <- as.data.table(expr[expr$symbol %in% sig,c("symbol",samples),drop=FALSE]);dt<-dt[,lapply(.SD,mean,na.rm=TRUE),by=symbol,.SDcols=samples]
mat<-as.matrix(dt[,-1]);rownames(mat)<-dt$symbol;storage.mode(mat)<-"numeric";z<-t(scale(t(mat)));z[!is.finite(z)]<-NA_real_
module9 <- colMeans(z,na.rm=TRUE)
d <- fread(file.path(base,"results/metadata_scores.csv"),data.table=FALSE)
d$IMMUNE_MODULE9_ENRICHMENT_CORE <- module9[d$gsm]
d$immune_z <- as.numeric(scale(d$IMMUNE_MODULE9_ENRICHMENT_CORE))
d$imbalance_z <- as.numeric(scale(d$Repair_Cilia_Imbalance))
d$case_primary <- as.integer(d$group==1)
primary <- d[d$group %in% c(1,2),]

continuous <- function(v) {
  q<-primary[complete.cases(primary[,c(v,"case_primary","baseline_grade","age_num","pack_years_num")]),]
  q$y<-as.numeric(scale(q[[v]]));fit<-lm(y~case_primary+baseline_grade+age_num+pack_years_num,q);co<-coef(summary(fit))["case_primary",];ci<-confint(fit,"case_primary")
  data.frame(axis=v,contrast="persistent_minus_regressive",n_sites=nrow(q),effect_size_SD=co["Estimate"],CI95_lower=ci[1],CI95_upper=ci[2],P=co["Pr(>|t|)"],analysis_method="adjusted lesion-site OLS; patient ID unavailable")
}
effects <- rbind(continuous("imbalance_z"),continuous("immune_z"));effects$FDR<-p.adjust(effects$P,"BH");fwrite(effects,file.path(out,"two_axis_outcome_effects.csv"))

fit_logistic <- function(label,form) {
  fit<-glm(form,family=binomial(),data=primary);terms<-setdiff(names(coef(fit)),c("(Intercept)","baseline_grade","age_num","pack_years_num"))
  rbindlist(lapply(terms,function(term){co<-coef(summary(fit))[term,];est<-co["Estimate"];se<-co["Std. Error"];data.frame(model=label,term=term,n_sites=nobs(fit),events=sum(model.response(model.frame(fit))),OR_per_SD=exp(est),CI95_lower=exp(est-1.96*se),CI95_upper=exp(est+1.96*se),P=co["Pr(>|z|)"],AIC=AIC(fit),converged=fit$converged,analysis_method="adjusted lesion-site logistic; patient ID unavailable")}))
}
mods <- rbind(
 fit_logistic("epithelial_axis_only",case_primary~imbalance_z+baseline_grade+age_num+pack_years_num),
 fit_logistic("immune_axis_only",case_primary~immune_z+baseline_grade+age_num+pack_years_num),
 fit_logistic("two_axis",case_primary~imbalance_z+immune_z+baseline_grade+age_num+pack_years_num))
mods$FDR<-p.adjust(mods$P,"BH");fwrite(mods,file.path(out,"two_axis_logistic_models.csv"))

ct1<-cor.test(d$imbalance_z,d$immune_z,method="pearson");ct2<-cor.test(d$imbalance_z,d$immune_z,method="spearman",exact=FALSE)
cors<-data.frame(method=c("Pearson","Spearman"),n_sites=nrow(d),correlation=c(ct1$estimate,ct2$estimate),P=c(ct1$p.value,ct2$p.value));cors$FDR<-p.adjust(cors$P,"BH");fwrite(cors,file.path(out,"epithelial_immune_axis_correlations.csv"))
coverage<-data.frame(signature="IMMUNE_MODULE9_ENRICHMENT_CORE",requested_n=length(sig),observed_n=nrow(mat),observed_genes=paste(sort(rownames(mat)),collapse=";"),missing_genes=paste(sort(setdiff(sig,rownames(mat))),collapse=";"));fwrite(coverage,file.path(out,"module9_coverage.csv"));fwrite(d,file.path(out,"metadata_two_axis_scores.csv"))

cols<-ifelse(primary$group==1,"#B23A48","#2878B5");png(file.path(out,"figures/GSE114489_two_axis_outcome.png"),1800,1550,res=250)
plot(primary$imbalance_z,primary$immune_z,pch=19,col=adjustcolor(cols,.75),xlab="Repair–Cilia Imbalance (z)",ylab="Immune Module 9 core (z)",main="GSE114489 baseline dysplasia: two-axis lesion fate")
abline(v=0,h=0,lty=2,col="grey60");legend("topright",legend=c("Persistent","Regressive"),col=c("#B23A48","#2878B5"),pch=19,bty="n");dev.off()

im<-effects[effects$axis=="imbalance_z",];mu<-effects[effects$axis=="immune_z",];two<-mods[mods$model=="two_axis",]
writeLines(c("# GSE114489 epithelial/immune two-axis analysis","",sprintf("The external Module 9 enrichment core covered %d/%d genes on GPL6244.",nrow(mat),length(sig)),
 sprintf("Persistent-minus-regressive epithelial imbalance = %.3f SD (95%% CI %.3f to %.3f), P=%.3g, FDR=%.3g.",im$effect_size_SD,im$CI95_lower,im$CI95_upper,im$P,im$FDR),
 sprintf("Persistent-minus-regressive immune Module 9 core = %.3f SD (95%% CI %.3f to %.3f), P=%.3g, FDR=%.3g.",mu$effect_size_SD,mu$CI95_lower,mu$CI95_upper,mu$P,mu$FDR),"",
 "The joint logistic model is exploratory and lesion-site based because public patient identifiers are unavailable.",
 "The 95-gene immune score is an externally reconstructed enrichment-annotated core of the reported 110-gene Module 9, not a claim of exact full membership."),file.path(out,"analysis_summary.md"))
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"));print(effects);print(mods)
