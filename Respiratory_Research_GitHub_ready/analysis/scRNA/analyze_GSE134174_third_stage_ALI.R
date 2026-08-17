#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(mgcv) })
set.seed(20260811)
project <- Sys.getenv("PROJECT_ROOT", unset = normalizePath(".", mustWork = FALSE))
root <- file.path(project, "02_scRNA/GSE134174")
out <- file.path(root, "third_stage_results")
dir.create(file.path(out, "figures"), recursive=TRUE, showWarnings=FALSE)

gene_sets <- list(
  DNA_DAMAGE_RESPONSE=c("H2AX","TP53BP1","ATM","ATR","CHEK1","CHEK2"),
  ABNORMAL_REPAIR=c("KRT8","KRT17","KRT19","CLDN4","SFN","LGALS3","KRT14"),
  KRT14_KRT17_REPAIR=c("KRT14","KRT17"),
  CILIA_CONSENSUS=c("FOXJ1","MCIDAS","GMNC","MYB","TP73","RFX2","RFX3","CCNO","CDC20B","DEUP1","TUBB4B","DNAH5","DNAH9","DNAH11","DNAI1","DNAI2","CCDC39","CCDC40","SPEF2","HYDIN","PIFO","CFAP43","CFAP44","RSPH1","RSPH4A","RSPH9"),
  MULTICILIOGENESIS=c("FOXJ1","MCIDAS","GMNC","MYB","TP73","RFX2","RFX3","CCNO","CDC20B","DEUP1","CEP78","CETN2","PLK4","STIL"),
  MATURE_CILIATED=c("FOXJ1","PIFO","TPPP3","CAPS","RSPH1"),
  KRT4_KRT13_TRANSITIONAL_FATE=c("KRT4","KRT13")
)
pb <- fread(file.path(root,"derived/invitro_target_pseudobulk_third_stage.csv"))
lib <- unique(pb[,.(donor,day,cluster_ident,library_size,n_cells)])[,.(library_size=sum(library_size),n_cells=sum(n_cells)),by=.(donor,day)]
cnt <- pb[,.(count=sum(count)),by=.(donor,day,gene)]; cnt <- merge(cnt,lib,by=c("donor","day"))
cnt[,log2_cpm:=log2((count+.5)/(library_size+1)*1e6)]
w <- dcast(cnt,donor+day+n_cells+library_size~gene,value.var="log2_cpm")
w[,day_numeric:=fifelse(day=="seed_day",-3,fifelse(day=="day_minus2",-2,as.numeric(sub("day_","",day))))]
for(g in intersect(unique(unlist(gene_sets)),names(w))) w[[paste0(g,"_Z")]] <- as.numeric(scale(w[[g]]))
for(nm in names(gene_sets)) {
  use <- paste0(intersect(gene_sets[[nm]],names(w)),"_Z")
  w[[nm]] <- if(length(use)>=2) rowMeans(w[,..use],na.rm=TRUE) else NA_real_
}
w[,Repair_Cilia_Imbalance:=as.numeric(scale(ABNORMAL_REPAIR))-as.numeric(scale(CILIA_CONSENSUS))]
for(g in c("FOXJ1","KRT14","KRT17","MCIDAS","GMNC","KRT4","KRT13")) if(g %in% names(w)) w[[paste0(g,"_score")]] <- as.numeric(scale(w[[g]]))
score_vars <- c(names(gene_sets),"Repair_Cilia_Imbalance",paste0(intersect(c("FOXJ1","KRT14","KRT17","MCIDAS","GMNC","KRT4","KRT13"),names(w)),"_score"))
fwrite(w[,c("donor","day","day_numeric","n_cells",score_vars),with=FALSE],file.path(out,"ALI_full_timecourse_scores.csv"))

models <- rbindlist(lapply(score_vars,function(v){
  z <- w[is.finite(get(v))]
  fit <- gam(as.formula(paste(v,"~ s(day_numeric,k=7) + donor")),data=z,method="REML")
  st <- summary(fit)$s.table[1,]
  data.frame(variable=v,n_donors=uniqueN(z$donor),n_observations=nrow(z),edf=unname(st["edf"]),F=unname(st["F"]),
             p_smooth=unname(st["p-value"]),adjusted_R2=summary(fit)$r.sq,method="GAM smooth plus donor fixed effect")
}))
models[,FDR:=p.adjust(p_smooth,"BH")]; fwrite(models,file.path(out,"ALI_full_timecourse_models.csv"))

mean_curve <- w[,lapply(.SD,mean,na.rm=TRUE),by=day_numeric,.SDcols=score_vars][order(day_numeric)]
peaks <- rbindlist(lapply(score_vars,function(v){x<-mean_curve[[v]]; data.frame(variable=v,max_day=mean_curve$day_numeric[which.max(x)],max_score=max(x),min_day=mean_curve$day_numeric[which.min(x)],min_score=min(x))}))
fwrite(peaks,file.path(out,"ALI_key_peak_timing.csv"))

w[,stage:=cut(day_numeric,breaks=c(-Inf,2,14,Inf),labels=c("early","intermediate","late"))]
stage_mean <- w[,lapply(.SD,mean,na.rm=TRUE),by=.(donor,stage),.SDcols=score_vars]
contrasts <- rbindlist(lapply(score_vars,function(v){
  wide<-dcast(stage_mean,donor~stage,value.var=v)
  diffs<-list(intermediate_vs_early=wide$intermediate-wide$early,late_vs_intermediate=wide$late-wide$intermediate,late_vs_early=wide$late-wide$early)
  rbindlist(lapply(names(diffs),function(k){d<-diffs[[k]];tt<-tryCatch(t.test(d),error=function(e)NULL);data.frame(variable=v,contrast=k,n_donors=sum(is.finite(d)),estimate=mean(d,na.rm=TRUE),CI95_lower=if(is.null(tt))NA_real_ else tt$conf.int[1],CI95_upper=if(is.null(tt))NA_real_ else tt$conf.int[2],p_value=if(is.null(tt))NA_real_ else tt$p.value)}))
}))
contrasts[,FDR:=p.adjust(p_value,"BH")]; fwrite(contrasts,file.path(out,"ALI_stage_contrasts_third_stage.csv"))

cor_targets <- c("ABNORMAL_REPAIR","KRT14_KRT17_REPAIR","CILIA_CONSENSUS","MULTICILIOGENESIS","MATURE_CILIATED","Repair_Cilia_Imbalance")
corrows <- rbindlist(lapply(cor_targets,function(v){
  donor_r <- w[,.(r=suppressWarnings(cor(KRT4_KRT13_TRANSITIONAL_FATE,get(v),method="spearman"))),by=donor]
  rz <- atanh(pmax(pmin(donor_r$r,.999),-.999)); tt<-tryCatch(t.test(rz),error=function(e)NULL)
  data.frame(reference="KRT4_KRT13_TRANSITIONAL_FATE",variable=v,n_donors=nrow(donor_r),mean_spearman=tanh(mean(rz)),
             CI95_lower=if(is.null(tt))NA_real_ else tanh(tt$conf.int[1]),CI95_upper=if(is.null(tt))NA_real_ else tanh(tt$conf.int[2]),
             p_value=if(is.null(tt))NA_real_ else tt$p.value,method="donor-specific Spearman; Fisher-z mean")
}))
corrows[,FDR:=p.adjust(p_value,"BH")]; fwrite(corrows,file.path(out,"ALI_transitional_cilia_correlations.csv"))

coverage <- rbindlist(lapply(names(gene_sets),function(nm){hit<-intersect(gene_sets[[nm]],names(w));data.frame(signature=nm,n_defined=length(gene_sets[[nm]]),n_detected=length(hit),detected=paste(hit,collapse=";"),missing=paste(setdiff(gene_sets[[nm]],hit),collapse=";"))}))
fwrite(coverage,file.path(out,"signature_coverage.csv"))

focus <- c("KRT14_KRT17_REPAIR","KRT4_KRT13_TRANSITIONAL_FATE","MULTICILIOGENESIS","MATURE_CILIATED","Repair_Cilia_Imbalance")
cols <- c("#E45756","#F28E2B","#59A14F","#4C78A8","#B279A2")
png(file.path(out,"figures/ALI_third_stage_full_trajectory.png"),width=2100,height=1400,res=250)
plot(range(mean_curve$day_numeric),range(as.matrix(mean_curve[,..focus])),type="n",xlab="ALI day",ylab="Mean gene-wise z score",main="GSE134174: full 20-timepoint repair-to-ciliogenesis trajectory")
for(i in seq_along(focus)) lines(mean_curve$day_numeric,mean_curve[[focus[i]]],type="b",lwd=2,pch=16,col=cols[i])
legend("topleft",legend=focus,col=cols,lwd=2,pch=16,cex=.78,bty="n");dev.off()

tr <- peaks[variable=="KRT4_KRT13_TRANSITIONAL_FATE"]; mc<-peaks[variable=="MULTICILIOGENESIS"]; mat<-peaks[variable=="MATURE_CILIATED"]
writeLines(c("# GSE134174 full ALI time-course — third stage","",
  sprintf("The processed ALI experiment contains %d donor-day observations from %d donors across %d distinct time points.",nrow(w),uniqueN(w$donor),uniqueN(w$day_numeric)),
  sprintf("KRT4/KRT13 transitional fate peaks at day %s; multiciliogenesis peaks at day %s; mature-ciliated score peaks at day %s.",tr$max_day,mc$max_day,mat$max_day),"",
  "GAMs use the full time axis with a nonlinear smooth and donor fixed effects. With only three donors, trajectory P values and confidence intervals are supportive rather than definitive.",
  "This is an unexposed differentiation series: it validates temporal ordering but cannot demonstrate exposure-induced trajectory blockade.",
  "The KRT4/KRT13 score is a prespecified two-gene suprabasal/transitional marker set; it was not trained on lesion outcome."
),file.path(out,"analysis_summary.md"))
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
cat(sprintf("Completed ALI third-stage trajectory: %d donor-days, %d time points\n",nrow(w),uniqueN(w$day_numeric)))
