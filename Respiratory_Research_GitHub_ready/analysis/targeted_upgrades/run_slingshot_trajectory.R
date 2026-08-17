#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(SingleCellExperiment);library(slingshot)})
set.seed(20260813)
project <- Sys.getenv("PROJECT_ROOT", unset = normalizePath(".", mustWork = FALSE));out<-file.path(project,"12_targeted_manuscript_upgrades/02_ALI_trajectory")
md<-fread(file.path(out,"ALI_lineage_metadata.csv"));x<-fread(file.path(out,"ALI_lineage_log2CPM.csv"));genes<-x$gene;mat<-as.matrix(x[,-1]);rownames(mat)<-genes
hvg<-fread(file.path(out,"ALI_HVG_1000.csv"))$gene;z<-t(scale(t(mat[hvg,,drop=FALSE])));z[!is.finite(z)]<-0
pca<-prcomp(t(z),rank.=20,center=FALSE,scale.=FALSE)$x
sce<-SingleCellExperiment(assays=list(logcounts=mat));reducedDim(sce,"PCA")<-pca
sce<-slingshot(sce,clusterLabels=md$cluster,reducedDim="PCA",start.clus="basal.colonies",end.clus="ciliated",allow.breaks=FALSE)
pt<-slingPseudotime(sce);w<-slingCurveWeights(sce);lineage<-which.max(colSums(w,na.rm=TRUE));pseudo<-pt[,lineage];weight<-w[,lineage]
if(cor(pseudo,md$day_numeric,use="complete.obs")<0)pseudo<-max(pseudo,na.rm=TRUE)-pseudo
pseudo<-(pseudo-min(pseudo,na.rm=TRUE))/(max(pseudo,na.rm=TRUE)-min(pseudo,na.rm=TRUE));md[,pseudotime:=pseudo];md[,lineage_weight:=weight]

sets<-list(ABNORMAL_REPAIR=c("KRT8","KRT17","KRT19","CLDN4","SFN","LGALS3","KRT14"),KRT14_KRT17_REPAIR=c("KRT14","KRT17"),CILIA_CONSENSUS=c("FOXJ1","MCIDAS","GMNC","MYB","TP73","RFX2","RFX3","CCNO","CDC20B","DEUP1","TUBB4B","DNAH5","DNAH9","DNAH11","DNAI1","DNAI2","CCDC39","CCDC40","SPEF2","HYDIN","PIFO","CFAP43","CFAP44","RSPH1","RSPH4A","RSPH9"),MULTICILIOGENESIS=c("FOXJ1","MCIDAS","GMNC","MYB","TP73","RFX2","RFX3","CCNO","CDC20B","DEUP1","CEP78","CETN2","PLK4","STIL"),MATURE_CILIATED=c("FOXJ1","PIFO","TPPP3","CAPS","RSPH1"),KRT4_KRT13_TRANSITIONAL_FATE=c("KRT4","KRT13"))
gene_z<-t(scale(t(mat)));gene_z[!is.finite(gene_z)]<-NA_real_
scores<-data.table(cell=md$cell)
for(nm in names(sets)){hit<-intersect(sets[[nm]],rownames(gene_z));scores[[nm]]<-colMeans(gene_z[hit,,drop=FALSE],na.rm=TRUE)}
scores[,ABNORMAL_REPAIR_Z:=as.numeric(scale(ABNORMAL_REPAIR))];scores[,CILIA_CONSENSUS_Z:=as.numeric(scale(CILIA_CONSENSUS))];scores[,Repair_Cilia_Imbalance:=ABNORMAL_REPAIR_Z-CILIA_CONSENSUS_Z]
for(g in c("FOXJ1","MCIDAS","GMNC"))scores[[g]]<-as.numeric(scale(mat[g,]))
scores<-cbind(md[,.(cell,donor,day,day_numeric,cluster,pseudotime,lineage_weight)],scores[,-1])
scores<-scores[is.finite(pseudotime)&lineage_weight>.20]
scores[,pseudotime_bin:=pmin(10L,floor(pseudotime*10)+1L)]
programs<-c(names(sets),"Repair_Cilia_Imbalance","FOXJ1","MCIDAS","GMNC")
bins<-scores[,c(list(n_cells=.N,pseudotime=mean(pseudotime),day_numeric=mean(day_numeric)),lapply(.SD,mean,na.rm=TRUE)),by=.(donor,pseudotime_bin),.SDcols=programs]
tests<-rbindlist(lapply(programs,function(v){ds<-bins[,.(slope=coef(lm(get(v)~pseudotime))[2],spearman=cor(get(v),pseudotime,method="spearman")),by=donor];tt<-tryCatch(t.test(ds$slope),error=function(e)NULL);data.table(program=v,n_donors=nrow(ds),slope_per_unit_pseudotime=mean(ds$slope),CI95_lower=if(is.null(tt))NA_real_ else tt$conf.int[1],CI95_upper=if(is.null(tt))NA_real_ else tt$conf.int[2],P=if(is.null(tt))NA_real_ else tt$p.value,mean_donor_spearman=mean(ds$spearman),donor_direction_consistency=ifelse(all(ds$slope>0),"all_positive",ifelse(all(ds$slope<0),"all_negative","mixed")))}))
tests[,FDR:=p.adjust(P,"BH")]
fwrite(md,file.path(out,"slingshot_cell_pseudotime.csv"));fwrite(scores,file.path(out,"slingshot_cell_program_scores.csv"));fwrite(bins,file.path(out,"donor_pseudotime_bin_scores.csv"));fwrite(tests,file.path(out,"trajectory_program_tests.csv"))

# ALI-only maturation candidate statistics; PML data are not read here.
excluded<-unique(c(unlist(sets[c("ABNORMAL_REPAIR","CILIA_CONSENSUS")]),"KRT14","KRT17","KRT4","KRT13","FOXJ1"));candidates<-setdiff(intersect(hvg,rownames(mat)),excluded)
gene_rows<-rbindlist(lapply(candidates,function(g){cell_values<-mat[g,match(scores$cell,md$cell)];tmp<-data.table(donor=scores$donor,bin=scores$pseudotime_bin,pseudotime=scores$pseudotime,value=cell_values)[,.(pseudotime=mean(pseudotime),value=mean(value)),by=.(donor,bin)];dr<-tmp[,.(spearman=cor(value,pseudotime,method="spearman"),slope=coef(lm(value~pseudotime))[2]),by=donor];fit<-lm(value~pseudotime+donor,data=tmp);co<-summary(fit)$coef["pseudotime",];data.table(gene=g,slope=co["Estimate"],P=co["Pr(>|t|)"],mean_spearman=mean(dr$spearman),min_abs_donor_spearman=min(abs(dr$spearman)),direction_consistent=(all(dr$spearman>0)|all(dr$spearman<0)),donor_spearman=paste(sprintf("%s:%.4f",dr$donor,dr$spearman),collapse=";"))}))
gene_rows[,FDR:=p.adjust(P,"BH")];gene_rows[,direction:=ifelse(slope>0,"maturation_up","maturation_down")]
selected<-gene_rows[direction_consistent==TRUE&abs(mean_spearman)>=.60&FDR<.05][order(direction,-abs(slope),-min_abs_donor_spearman)][,head(.SD,50),by=direction]
fwrite(gene_rows,file.path(out,"ALI_maturation_gene_statistics.csv"));fwrite(selected,file.path(project,"12_targeted_manuscript_upgrades/03_ALI_maturation_projection/ALI_maturation_signature.csv"))

pdf(file.path(project,"12_targeted_manuscript_upgrades/05_figures/candidate_figures/Candidate_Figure_B_ALI_trajectory.pdf"),width=11,height=8.5,useDingbats=FALSE)
layout(matrix(1:4,2,2));cols<-c("#A85A44","#D08C60","#E9C46A","#62A7A1","#326273","#235789","#1D3557","#0B132B");plot(pca[,1],pca[,2],col=cols[match(md$cluster,unique(md$cluster))],pch=16,cex=.35,xlab="PC1",ylab="PC2",main="A  Slingshot ALI ciliated lineage");legend("topright",legend=unique(md$cluster),col=cols[seq_along(unique(md$cluster))],pch=16,cex=.55,bty="n");plot(pca[,1],pca[,2],col=hcl.colors(100,"viridis")[pmax(1,round(md$pseudotime*99)+1)],pch=16,cex=.35,xlab="PC1",ylab="PC2",main="B  Slingshot pseudotime");focus<-c("Repair_Cilia_Imbalance","KRT14_KRT17_REPAIR","KRT4_KRT13_TRANSITIONAL_FATE","MULTICILIOGENESIS","MATURE_CILIATED");plot(0,0,type="n",xlim=c(0,1),ylim=range(as.matrix(bins[,..focus])),xlab="Pseudotime",ylab="Donor-bin mean score",main="C  Repair-to-cilia programs");for(i in seq_along(focus)){curve<-bins[,.(x=mean(pseudotime),y=mean(get(focus[i]))),by=pseudotime_bin][order(x)];lines(curve$x,curve$y,col=c("#8C2F39","#C8553D","#E9C46A","#2A9D8F","#264653")[i],lwd=2)};legend("topright",legend=focus,col=c("#8C2F39","#C8553D","#E9C46A","#2A9D8F","#264653"),lwd=2,cex=.55,bty="n");ri<-bins[,.(x=mean(pseudotime),y=mean(Repair_Cilia_Imbalance),lo=min(Repair_Cilia_Imbalance),hi=max(Repair_Cilia_Imbalance)),by=pseudotime_bin][order(x)];plot(ri$x,ri$y,type="b",pch=16,lwd=2,col="#8C2F39",xlab="Pseudotime",ylab="Repair-Cilia Imbalance",main="D  Donor-aware trajectory");segments(ri$x,ri$lo,ri$x,ri$hi,col="#8C2F39");dev.off()
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"));cat("Slingshot lineage cells:",nrow(scores),"; maturation genes:",nrow(selected),"\n")
