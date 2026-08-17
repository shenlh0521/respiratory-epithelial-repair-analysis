#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(DWLS); library(quadprog)})
set.seed(20260813)
project <- Sys.getenv("PROJECT_ROOT", unset = normalizePath(".", mustWork = FALSE))
out <- file.path(project,"12_targeted_manuscript_upgrades/01_deconvolution_robustness")
states <- c("Basal_repair","Ciliated","Secretory","Transitional")

nnls_qp <- function(S,b){D<-crossprod(S)+diag(1e-8,ncol(S));d<-crossprod(S,b);z<-solve.QP(D,d,diag(ncol(S)),rep(0,ncol(S)))$solution;names(z)<-colnames(S);z/sum(z)}
row_minmax <- function(x){lo<-apply(x,1,min);hi<-apply(x,1,max);den<-hi-lo;den[den==0]<-NA;z<-(x-lo)/den;z[!is.finite(z)]<-0;z}
quiet_dwls <- function(S,b,seed){set.seed(seed);capture.output(ans<-tryCatch(solveDampenedWLS(S,b),error=function(e)e));if(inherits(ans,"error"))return(rep(NA_real_,ncol(S)));ans<-as.numeric(ans);names(ans)<-colnames(S);ans}

sig<-fread(file.path(out,"lodo_reference_signatures.csv"));mix<-fread(file.path(out,"lodo_mixture_expression.csv"));truth<-fread(file.path(out,"lodo_mixture_truth.csv"));pred<-list();counter<-0L
for(donor in unique(sig$held_out_donor)){
 s<-sig[held_out_donor==donor]
 sl<-dcast(s,gene~state,value.var="log2_cpm");sl_genes<-sl$gene;Sl<-as.matrix(sl[,..states]);rownames(Sl)<-sl_genes
 sc<-dcast(s,gene~state,value.var="cpm");sc_genes<-sc$gene;Sc<-as.matrix(sc[,..states]);rownames(Sc)<-sc_genes
 m<-mix[held_out_donor==donor]
 ml<-dcast(m,gene~mixture_id,value.var="log2_cpm");mix_columns<-setdiff(names(ml),"gene");Ml<-as.matrix(ml[match(rownames(Sl),gene),..mix_columns]);rownames(Ml)<-rownames(Sl)
 mc<-dcast(m,gene~mixture_id,value.var="cpm");mix_columns_cpm<-setdiff(names(mc),"gene");Mc<-as.matrix(mc[match(rownames(Sc),gene),..mix_columns_cpm]);rownames(Mc)<-rownames(Sc)
 A<-row_minmax(Sl);Y<-row_minmax(Ml)
 for(j in seq_len(ncol(Ml))){counter<-counter+1L;id<-colnames(Ml)[j];zn<-tryCatch(nnls_qp(A,Y[,j]),error=function(e)rep(NA_real_,4));names(zn)<-states;zd<-quiet_dwls(Sc,Mc[,j],20260813+counter);names(zd)<-states;pred[[length(pred)+1]]<-data.table(mixture_id=id,held_out_donor=donor,method="NNLS_LODO",state=states,predicted=zn);pred[[length(pred)+1]]<-data.table(mixture_id=id,held_out_donor=donor,method="DWLS",state=states,predicted=zd)}
}
pred<-rbindlist(pred);bench<-merge(truth,pred,by=c("mixture_id","held_out_donor","state"),all.x=TRUE)
metric<-function(z,label){ok<-is.finite(z$observed)&is.finite(z$predicted);if(sum(ok)<3)return(data.table(scope=label,n=nrow(z),n_complete=sum(ok),pearson=NA_real_,spearman=NA_real_,MAE=NA_real_,RMSE=NA_real_));data.table(scope=label,n=nrow(z),n_complete=sum(ok),pearson=cor(z$observed[ok],z$predicted[ok]),spearman=cor(z$observed[ok],z$predicted[ok],method="spearman"),MAE=mean(abs(z$predicted[ok]-z$observed[ok])),RMSE=sqrt(mean((z$predicted[ok]-z$observed[ok])^2)))}
metrics<-rbindlist(lapply(unique(bench$method),function(method_name){z<-bench[method==method_name];rbind(cbind(method=method_name,metric(z,"overall")),rbindlist(lapply(states,function(state_name)cbind(method=method_name,metric(z[state==state_name],state_name)))))}))
status<-merge(metrics[scope=="overall",.(method,overall_pearson=pearson,overall_spearman=spearman,overall_MAE=MAE)],metrics[scope!="overall",.(min_state_pearson=min(pearson,na.rm=TRUE),max_state_MAE=max(MAE,na.rm=TRUE)),by=method],by="method")
status[,benchmark_status:=ifelse(!is.finite(overall_pearson)|!is.finite(overall_spearman)|!is.finite(overall_MAE)|!is.finite(min_state_pearson)|!is.finite(max_state_MAE),"METHOD_FAILED",ifelse(overall_pearson>=.70&overall_spearman>=.70&overall_MAE<=.15&min_state_pearson>=.30&max_state_MAE<=.25,"PASS","FAIL"))]
status<-rbind(status,data.table(method="MuSiC",overall_pearson=NA_real_,overall_spearman=NA_real_,overall_MAE=NA_real_,min_state_pearson=NA_real_,max_state_MAE=NA_real_,benchmark_status="METHOD_FAILED"),fill=TRUE)
fwrite(pred,file.path(out,"pseudobulk_predictions.csv"));fwrite(bench,file.path(out,"pseudobulk_benchmark_full_predictions.csv"));fwrite(metrics,file.path(out,"pseudobulk_benchmark.csv"));fwrite(status,file.path(out,"benchmark_method_status.csv"))

# GSE320381 projection only for benchmark-passing DWLS; original NNLS remains unchanged.
orig<-fread(file.path(project,"01_PML_extended/GSE320381/composition_adjusted_results/composition_adjusted_severity_models.csv"));meta<-fread(file.path(project,"01_PML_extended/GSE320381/results/biopsy_metadata_scores.csv"));gexpr<-fread(file.path(out,"GSE320381_marker_expression.csv"));ref<-fread(file.path(project,"02_scRNA/GSE134174/third_stage_results/deconvolution_reference/GSE134174_never_smoker_reference_logCPM.csv"))
genes<-intersect(ref$gene,unique(gexpr$gene));S<-2^as.matrix(ref[match(genes,gene),..states])-1;rownames(S)<-genes
fractions<-list()
if(status[method=="DWLS",benchmark_status]=="PASS")for(i in seq_along(unique(gexpr$sample_id))){sid<-unique(gexpr$sample_id)[i];z<-gexpr[sample_id==sid&gene%chin%genes];b<-z$cpm[match(genes,z$gene)];names(b)<-genes;est<-quiet_dwls(S,b,30300000+i);fractions[[i]]<-data.table(sample_id=sid,method="DWLS",state=states,contribution=est)}
frac_long<-if(length(fractions))rbindlist(fractions)else data.table(sample_id=character(),method=character(),state=character(),contribution=numeric());fwrite(frac_long,file.path(out,"GSE320381_relative_contributions_new_methods.csv"))

cluster_fit<-function(d,covs){vars<-c("imbalance_z","grade_ordinal","patient_id",covs);z<-d[complete.cases(d[,..vars])];X<-cbind(intercept=1,as.matrix(z[,c("grade_ordinal",covs),with=FALSE]));y<-z$imbalance_z;bread<-solve(crossprod(X));beta<-bread%*%crossprod(X,y);resid<-as.numeric(y-X%*%beta);groups<-as.character(z$patient_id);meat<-matrix(0,ncol(X),ncol(X));for(g in unique(groups)){idx<-which(groups==g);u<-crossprod(X[idx,,drop=FALSE],resid[idx]);meat<-meat+tcrossprod(u)};G<-uniqueN(groups);N<-nrow(z);P<-ncol(X);vc<-(G/(G-1))*((N-1)/(N-P))*bread%*%meat%*%bread;est<-beta[2];se<-sqrt(vc[2,2]);crit<-qt(.975,G-1);data.table(beta=as.numeric(est),lower=as.numeric(est-crit*se),upper=as.numeric(est+crit*se),P=as.numeric(2*pt(-abs(est/se),G-1)),n_samples=N,n_subjects=G)}
u<-orig[model=="unadjusted"];a<-orig[model=="composition_adjusted"];unadj<-u$effect_size_SD_per_grade[1]
comparison<-rbind(data.table(method="NNLS",benchmark_status="ORIGINAL_METHOD",severity_beta_unadjusted=unadj,severity_CI_unadjusted=sprintf("%.6f to %.6f",u$CI95_lower,u$CI95_upper),severity_P_unadjusted=u$P,severity_beta_adjusted=a$effect_size_SD_per_grade,severity_CI_adjusted=sprintf("%.6f to %.6f",a$CI95_lower,a$CI95_upper),severity_P_adjusted=a$P,attenuation_percent=100*(unadj-a$effect_size_SD_per_grade)/unadj,notes="Original 24-marker NNLS relative projection"),data.table(method="MuSiC",benchmark_status="METHOD_FAILED",severity_beta_unadjusted=unadj,severity_CI_unadjusted=sprintf("%.6f to %.6f",u$CI95_lower,u$CI95_upper),severity_P_unadjusted=u$P,severity_beta_adjusted=NA_real_,severity_CI_adjusted=NA_character_,severity_P_adjusted=NA_real_,attenuation_percent=NA_real_,notes="Official package source unavailable after bounded retries"))
if(nrow(frac_long)){fw<-dcast(frac_long,sample_id~state,value.var="contribution");d<-merge(meta,fw,by="sample_id");d[,imbalance_z:=as.numeric(scale(Repair_Cilia_Imbalance))];fit<-cluster_fit(d,c("Basal_repair","Ciliated","Secretory"));dw<-data.table(method="DWLS",benchmark_status=status[method=="DWLS",benchmark_status],severity_beta_unadjusted=unadj,severity_CI_unadjusted=sprintf("%.6f to %.6f",u$CI95_lower,u$CI95_upper),severity_P_unadjusted=u$P,severity_beta_adjusted=fit$beta,severity_CI_adjusted=sprintf("%.6f to %.6f",fit$lower,fit$upper),severity_P_adjusted=fit$P,attenuation_percent=100*(unadj-fit$beta)/unadj,notes="DWLS relative estimated contributions; Transitional omitted")}else dw<-data.table(method="DWLS",benchmark_status=status[method=="DWLS",benchmark_status],severity_beta_unadjusted=unadj,severity_CI_unadjusted=sprintf("%.6f to %.6f",u$CI95_lower,u$CI95_upper),severity_P_unadjusted=u$P,severity_beta_adjusted=NA_real_,severity_CI_adjusted=NA_character_,severity_P_adjusted=NA_real_,attenuation_percent=NA_real_,notes="Excluded because benchmark did not pass")
comparison<-rbind(comparison,dw,fill=TRUE);comparison[,FDR:=p.adjust(severity_P_adjusted,"BH")];fwrite(comparison,file.path(out,"deconvolution_method_comparison.csv"));print(status);print(comparison)
