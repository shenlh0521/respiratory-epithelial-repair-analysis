#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(yaml)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

set.seed(20260813)
project <- Sys.getenv("PROJECT_ROOT", unset = normalizePath(".", mustWork = FALSE))
upgrade <- file.path(project, "12_targeted_manuscript_upgrades")
base <- file.path(upgrade, "04_DECAMP_validation")
figdir <- file.path(upgrade, "05_figures/candidate_figures")
legacy <- file.path(project, "files-mentioned-by-the-user-users/outputs/lung_bulla_cilia_cancer")

parse_soft <- function(path) {
  lines <- readLines(gzfile(path), warn=FALSE); s <- which(grepl("^\\^SAMPLE = ",lines)); e <- c(s[-1]-1L,length(lines))
  rbindlist(lapply(seq_along(s),function(i){
    z<-lines[s[i]:e[i]];one<-function(p){h<-z[grepl(p,z)];if(length(h))sub(p,"",h[1])else NA_character_}
    ch<-sub("^!Sample_characteristics_ch1 = ","",z[grepl("^!Sample_characteristics_ch1 = ",z)])
    v<-setNames(sub("^[^:]+: ","",ch),tolower(sub(":.*$","",ch)))
    data.frame(gsm=one("^!Sample_geo_accession = "),sample_id=one("^!Sample_title = "),histology=v["genotype"],
               binary=v["treatment"],batch=v["batch"],platform=one("^!Sample_platform_id = "),stringsAsFactors=FALSE)
  }),fill=TRUE)
}

cluster_fit <- function(d,outcome) {
  z<-d[complete.cases(d[,c(outcome,"dysplasia","platform","patient_id")]),]
  y<-as.numeric(scale(z[[outcome]]));X0<-model.matrix(~dysplasia+platform,z);q<-qr(X0);X<-X0[,sort(q$pivot[seq_len(q$rank)]),drop=FALSE]
  j<-match("dysplasia",colnames(X));beta<-solve(crossprod(X),crossprod(X,y));resid<-as.numeric(y-X%*%beta);bread<-solve(crossprod(X))
  groups<-split(seq_len(nrow(z)),factor(z$patient_id));meat<-matrix(0,ncol(X),ncol(X));for(idx in groups){u<-as.numeric(crossprod(X[idx,,drop=FALSE],resid[idx]));meat<-meat+tcrossprod(u)}
  G<-length(groups);N<-nrow(z);P<-ncol(X);vc<-(G/(G-1))*((N-1)/(N-P))*bread%*%meat%*%bread;est<-beta[j];se<-sqrt(vc[j,j]);crit<-qt(.975,G-1)
  data.frame(endpoint=outcome,comparison="dysplasia_or_worse_vs_non_dysplasia",estimate_SD=est,std_error=se,CI95_lower=est-crit*se,CI95_upper=est+crit*se,P=2*pt(abs(est/se),G-1,lower.tail=FALSE),n_samples=N,n_subjects=G,model="standardized score ~ dysplasia + platform",inference="patient-cluster robust SE",stringsAsFactors=FALSE)
}

md<-as.data.frame(parse_soft(file.path(base,"GSE300258_family.soft.gz")));md$patient_id<-sub("-.*$","",md$sample_id);md$dysplasia<-as.integer(md$binary=="dysplasia")
x<-fread(cmd=sprintf("gzip -cd %s",shQuote(file.path(base,"GSE300258_DECAMP_bronchial_biopsies_WSI_GE_filtered_counts.tsv.gz"))),data.table=FALSE,check.names=FALSE);names(x)[1]<-"ensembl_id"
ids<-sub("\\.[0-9]+$","",x$ensembl_id);cols<-intersect(md$sample_id,names(x));counts<-as.matrix(x[,cols,drop=FALSE]);storage.mode(counts)<-"numeric";rownames(counts)<-ids
dge<-calcNormFactors(DGEList(counts));logcpm<-cpm(dge,log=TRUE,prior.count=1)
valid<-intersect(unique(rownames(logcpm)),keys(org.Hs.eg.db,keytype="ENSEMBL"));mp<-AnnotationDbi::select(org.Hs.eg.db,keys=valid,keytype="ENSEMBL",columns="SYMBOL");mp<-mp[!is.na(mp$SYMBOL)&nzchar(mp$SYMBOL),];cnt<-table(mp$ENSEMBL);mp<-mp[mp$ENSEMBL%in%names(cnt[cnt==1]),]
sym<-toupper(setNames(mp$SYMBOL,mp$ENSEMBL)[rownames(logcpm)]);keep<-!is.na(sym)&nzchar(sym);dt<-as.data.table(logcpm[keep,,drop=FALSE]);dt[,symbol:=sym[keep]];gdt<-dt[,lapply(.SD,mean,na.rm=TRUE),by=symbol];gmat<-as.matrix(gdt[,-1]);storage.mode(gmat)<-"numeric";rownames(gmat)<-gdt$symbol

# Existing frozen Repair-Cilia definition.
sets<-yaml.load_file(file.path(legacy,"config/gene_sets.yaml"))$gene_sets
gene_z<-t(scale(t(gmat)));gene_z[!is.finite(gene_z)]<-NA_real_
repair_hit<-intersect(toupper(sets$ABNORMAL_REPAIR),rownames(gene_z));cilia_hit<-intersect(toupper(sets$CILIA_CONSENSUS),rownames(gene_z))
repair<-colMeans(gene_z[repair_hit,,drop=FALSE],na.rm=TRUE);cilia<-colMeans(gene_z[cilia_hit,,drop=FALSE],na.rm=TRUE)
imb<-as.numeric(scale(repair))-as.numeric(scale(cilia));names(imb)<-colnames(gmat)

# Frozen ALI maturation percentile-rank score.
sig<-fread(file.path(upgrade,"03_ALI_maturation_projection/ALI_maturation_signature.csv"),data.table=FALSE);up<-toupper(sig$gene[sig$direction=="maturation_up"]);down<-toupper(sig$gene[sig$direction=="maturation_down"])
ranks<-apply(gmat,2,function(v)rank(v,ties.method="average")/length(v));rownames(ranks)<-rownames(gmat);colnames(ranks)<-colnames(gmat);up_hit<-intersect(up,rownames(ranks));down_hit<-intersect(down,rownames(ranks))
if(length(up_hit)<25||length(down_hit)<2)stop("Frozen ALI signature failed coverage")
maturation<-colMeans(ranks[up_hit,,drop=FALSE])-colMeans(ranks[down_hit,,drop=FALSE]);names(maturation)<-colnames(gmat)

md$Repair_Cilia_Imbalance<-imb[md$sample_id];md$ALI_MATURATION_SCORE<-maturation[md$sample_id]
res<-rbind(cluster_fit(md,"Repair_Cilia_Imbalance"),cluster_fit(md,"ALI_MATURATION_SCORE"));res$FDR<-p.adjust(res$P,"BH")
expected<-c(Repair_Cilia_Imbalance="positive",ALI_MATURATION_SCORE="negative")
res$expected_direction<-expected[res$endpoint]
res$status<-ifelse(res$endpoint=="Repair_Cilia_Imbalance",ifelse(res$estimate_SD>0 & res$FDR<.05,"SUPPORTS",ifelse(res$estimate_SD>0,"DIRECTIONAL_ONLY",ifelse(res$FDR<.05,"INCONSISTENT","NULL"))),ifelse(res$estimate_SD<0 & res$FDR<.05,"SUPPORTS",ifelse(res$estimate_SD<0,"DIRECTIONAL_ONLY",ifelse(res$FDR<.05,"INCONSISTENT","NULL"))))
res$critical_limitation<-"GEO batch perfectly collinear with dysplasia; biology and batch are not identifiable"

fwrite(md,file.path(base,"GSE300258_frozen_scores.csv"));fwrite(res,file.path(base,"GSE300258_validation_results.csv"))
fwrite(data.frame(signature=c("ABNORMAL_REPAIR","CILIA_CONSENSUS","ALI_MATURATION_UP","ALI_MATURATION_DOWN"),requested_n=c(length(sets$ABNORMAL_REPAIR),length(sets$CILIA_CONSENSUS),length(up),length(down)),observed_n=c(length(repair_hit),length(cilia_hit),length(up_hit),length(down_hit)),observed_genes=c(paste(repair_hit,collapse=";"),paste(cilia_hit,collapse=";"),paste(up_hit,collapse=";"),paste(down_hit,collapse=";"))),file.path(base,"signature_coverage.csv"))

pdf(file.path(figdir,"Candidate_Figure_E_GSE300258_external_validation.pdf"),width=11,height=6.5,useDingbats=FALSE)
par(mfrow=c(1,2),mar=c(5,5,5,1),oma=c(0,0,2.5,0))
for(v in c("Repair_Cilia_Imbalance","ALI_MATURATION_SCORE")){
  z<-md[,c(v,"dysplasia")];boxplot(z[[v]]~factor(z$dysplasia,levels=c(0,1),labels=c("Non-dysplasia","Dysplasia")),col=c("#BDBDBD","#D55E00"),ylab=v,xlab="",main=v,cex.main=.95);stripchart(z[[v]]~factor(z$dysplasia,levels=c(0,1)),vertical=TRUE,method="jitter",add=TRUE,pch=16,col="#00000066")
  rr<-res[res$endpoint==v,];mtext(sprintf("Effect %.2f SD (95%% CI %.2f to %.2f); FDR %.3g",rr$estimate_SD,rr$CI95_lower,rr$CI95_upper,rr$FDR),side=3,line=.2,cex=.75)
}
mtext("Candidate only: platform-adjusted patient-cluster inference; GEO batch is perfectly collinear with dysplasia",outer=TRUE,side=3,line=.7,cex=.85,col="#B2182B")
dev.off()
capture.output(sessionInfo(),file=file.path(base,"sessionInfo.txt"))
cat("Completed one-time frozen GSE300258 validation\n");print(res);print(table(md$histology));print(table(md$binary,md$batch))
