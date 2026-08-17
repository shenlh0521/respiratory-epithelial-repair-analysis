#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(openxlsx)})
args <- commandArgs(trailingOnly=FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
code_dir <- dirname(normalizePath(script_arg))
base <- dirname(code_dir)
repo_root <- normalizePath(file.path(base, "..", ".."), mustWork=FALSE)
root <- Sys.getenv("PROJECT_ROOT", unset = repo_root)
out <- file.path(base, "outputs"); dir.create(out, recursive=TRUE, showWarnings=FALSE)
C<-c(repair="#D55E00",imbalance="#A44200",transition="#CC79A7",cilia="#0072B2",ciliating="#56B4E9",maturation="#005F73",neutral="#666666",null="#C9CDD0",stress="#7A5195",ink="#1F2529")
theme_set(theme_classic(base_family="Arial",base_size=9)+theme(plot.title=element_text(size=11,face="bold"),plot.tag=element_text(size=13,face="bold"),plot.tag.position=c(.01,.99),axis.title=element_text(size=9),axis.text=element_text(size=8),legend.text=element_text(size=8),legend.title=element_text(size=8),panel.grid.major.y=element_line(color="#E6E8EA",linewidth=.25)))
rd<-function(rel)read.csv(file.path(root,rel),check.names=FALSE)
saveS<-function(p,n,w=10,h=6){if(inherits(p,"ggplot")){p<-p+theme(plot.margin=margin(16,16,16,20))};base<-file.path(out,paste0("Figure_S",n));cairo_pdf(paste0(base,".pdf"),width=w,height=h,family="Arial");print(p);dev.off();svg(paste0(base,".svg"),width=w,height=h,family="Arial",onefile=TRUE);print(p);dev.off();png(paste0(base,".png"),width=w,height=h,units="in",res=600,type="cairo",family="Arial");print(p);dev.off()}
forest<-function(d,label,est="estimate",lo="lower",hi="upper",fdr="FDR",title="",xlabel="Effect (95% CI)",color=C["imbalance"]){d$lab<-d[[label]];d$e<-as.numeric(d[[est]]);d$l<-as.numeric(d[[lo]]);d$h<-as.numeric(d[[hi]]);d$lab<-factor(d$lab,levels=rev(unique(d$lab)));ggplot(d,aes(e,lab))+geom_vline(xintercept=0,linetype=2,color="#777",linewidth=.35)+geom_errorbarh(aes(xmin=l,xmax=h),height=.12,color=color,linewidth=.55)+geom_point(color=color,size=2)+labs(title=title,x=xlabel,y=NULL)+theme(panel.grid.major.y=element_blank())}
boxpts<-function(d,g,v,title,ylabel){pal<-c("Non-dysplasia"=unname(C["neutral"]),Dysplasia=unname(C["repair"]),Regressive=unname(C["cilia"]),Persistent=unname(C["repair"]));ggplot(d,aes(.data[[g]],as.numeric(.data[[v]]),color=.data[[g]],fill=.data[[g]]))+geom_boxplot(outlier.shape=NA,alpha=.28,width=.55,linewidth=.45)+geom_jitter(width=.08,size=.8,alpha=.35)+scale_color_manual(values=pal,guide="none")+scale_fill_manual(values=pal,guide="none")+labs(title=title,x=NULL,y=ylabel)}

# Frozen sources.
prov<-rd("11_BMC_Pulmonary_Medicine_submission/10_source_archive/provenance/final_result_provenance.csv");s<-rd("01_PML_extended/GSE109743/pml_extended_scores.csv")
eff<-rd("01_PML_extended/GSE109743/pml_subtype_effects.csv");corr<-rd("01_PML_extended/GSE109743/repair_cilia_correlations.csv")

# S1 PML overview/QC
datasets<-data.frame(dataset=c("GSE109743","GSE320381","GSE33479","GSE114489","GSE79210"),samples=c(295,108,122,38,51),participants=c(49,33,77,NA,23),role=c("Discovery","Severity replication","Severity replication","Lesion outcome","Field-brushing sensitivity"))
p1a<-ggplot(datasets,aes(samples,reorder(dataset,samples),fill=role))+geom_col(width=.65)+geom_text(aes(label=samples),hjust=-.2,size=3)+scale_fill_manual(values=c(Discovery=C["imbalance"],"Severity replication"=C["maturation"],"Lesion outcome"=C["transition"],"Field-brushing sensitivity"=C["neutral"]))+coord_cartesian(xlim=c(0,330))+labs(title="PML datasets",x="Samples",y=NULL,fill=NULL)
cnt<-as.data.frame(table(s$molecular_subtype));names(cnt)<-c("subtype","n");p1b<-ggplot(cnt,aes(n,reorder(subtype,n)))+geom_col(fill=C["imbalance"],width=.65)+geom_text(aes(label=n),hjust=-.2,size=3)+coord_cartesian(xlim=c(0,max(cnt$n)*1.2))+labs(title="GSE109743 molecular subtypes",x="Biopsies",y=NULL)
saveS(p1a+p1b+plot_annotation(tag_levels="A"),1)

# S2 score construction and rank sensitivity
flow<-data.frame(x=1:3,label=c("Z(ABNORMAL_REPAIR)","minus","Z(CILIA_CONSENSUS)"));p2a<-ggplot(flow)+geom_rect(aes(xmin=x-.38,xmax=x+.38,ymin=.42,ymax=.70),fill="white",color=c(C["repair"],C["neutral"],C["cilia"]),linewidth=.6)+geom_text(aes(x=x,y=.56,label=label),size=3)+coord_cartesian(xlim=c(.45,3.55),ylim=c(.2,.85))+labs(title="Repair-Cilia Imbalance")+theme_void(base_family="Arial")+theme(plot.title=element_text(size=11,face="bold"))
q<-eff[eff$outcome=="Repair_Cilia_Imbalance_rank"&grepl("_vs_Normal",eff$contrast),];q$feature<-gsub("_vs_Normal","",q$contrast);q$lower<-q$CI95_lower;q$upper<-q$CI95_upper
p2b<-forest(q,"feature",title="Rank-based subtype sensitivity",xlabel="Rank-normalized effect")
saveS(p2a+p2b+plot_annotation(tag_levels="A"),2)

# S3 subtype markers
q<-eff[eff$contrast=="Proliferative_vs_Normal"&eff$outcome%in%c("CILIA_CONSENSUS","ABNORMAL_REPAIR","FOXJ1","KRT14","KRT17","Repair_Cilia_Imbalance"),];q$feature<-q$outcome;q$lower<-q$CI95_lower;q$upper<-q$CI95_upper
saveS(forest(q,"feature",title="GSE109743: proliferative vs normal",xlabel="Adjusted effect (95% CI)"),3)

# S4 correlations and severity sensitivity
cc<-corr[corr$method=="pearson",];cc$feature<-cc$stratum;cc$estimate<-cc$correlation;cc$lower<-cc$CI95_lower;cc$upper<-cc$CI95_upper
gr<-rd("01_PML_extended/GSE109743/pml_grade_trends.csv");gr<-gr[gr$outcome%in%c("CILIA_CONSENSUS","ABNORMAL_REPAIR","FOXJ1","KRT14","KRT17","Repair_Cilia_Imbalance"),];gr$feature<-gr$outcome;gr$lower<-gr$CI95_lower;gr$upper<-gr$CI95_upper
saveS(forest(cc,"feature",title="Repair versus cilia correlation",xlabel="Pearson correlation")+forest(gr,"feature",title="Per-grade effects",xlabel="Effect per grade")+plot_annotation(tag_levels="A"),4,11,6)

# S5/S6 independent cohort full effects
for(z in list(c(5,"GSE320381"),c(6,"GSE33479"))){n<-as.integer(z[1]);ds<-z[2];q<-prov[prov$dataset==ds&grepl("histologic severity",prov$analysis),];q$feature<-q$feature;q$lower<-q$lower_CI;q$upper<-q$upper_CI;saveS(forest(q,"feature",title=paste(ds,"full severity results"),xlabel="Effect per histologic grade"),n,10,7)}

# S7 field null
q<-rd("01_PML_extended/GSE79210/results/baseline_persistent_vs_regressive.csv");q$feature<-q$variable;q$estimate<-q$mean_difference;q$lower<-q$CI95_lower;q$upper<-q$CI95_upper
saveS(forest(q,"feature",title="GSE79210 field-brushing outcome",xlabel="Persistent − regressive"),7,10,7)

# S8 persistence models
q<-rd("01_PML_extended/GSE114489/results/persistent_vs_regressive_effects.csv");q<-q[q$model=="adjusted_baseline_grade_age_packyears",];q$feature<-q$variable;q$estimate<-q$estimate_SD;q$lower<-q$CI95_lower;q$upper<-q$CI95_upper
saveS(forest(q,"feature",title="GSE114489 adjusted lesion-fate effects",xlabel="Standardized difference"),8,10,7)

# S9 scRNA state QC
st<-rd("02_scRNA/GSE134174/results/cell_state_signature_means.csv");p9a<-ggplot(st,aes(n_donors,reorder(subcluster_ident,n_donors)))+geom_col(fill=C["neutral"],width=.65)+labs(title="Donor representation by state",x="Donors",y=NULL)
comp<-rd("02_scRNA/GSE134174/results/smoking_cell_state_composition_effects.csv");comp<-comp[comp$state!="All_epithelial",];comp$feature<-comp$state;comp$estimate<-comp$estimate_percentage_points;comp$lower<-comp$CI95_lower_percentage_points;comp$upper<-comp$CI95_upper_percentage_points
p9b<-forest(comp,"feature",title="Heavy vs never-smoker composition",xlabel="Percentage-point difference",color=C["neutral"])
saveS(p9a+p9b+plot_annotation(tag_levels="A"),9,11,6)

# S10 scRNA state signatures
features<-c("ABNORMAL_REPAIR","KRT14_KRT17_REPAIR","CILIA_CONSENSUS","MULTICILIOGENESIS","MATURE_CILIATED","Repair_Cilia_Imbalance");mm<-do.call(rbind,lapply(features,function(f)data.frame(state=st$subcluster_ident,feature=f,value=st[[f]])));p10<-ggplot(mm,aes(feature,state,fill=value))+geom_tile(color="white")+scale_fill_gradient2(low=C["cilia"],mid="white",high=C["repair"],midpoint=0)+scale_x_discrete(labels=function(x)gsub("_"," ",x))+labs(title="GSE134174 state-level programs",x=NULL,y=NULL,fill="Mean score")+theme(axis.text.x=element_text(angle=30,hjust=1),panel.grid=element_blank())
saveS(p10,10,10,7)

# S11 supportive GSE233145
pr<-rd("02_scRNA/SECOND_VALIDATION/GSE233145/results/cell_state_signature_profiles.csv");pr<-pr[pr$time_point=="day 28",];mm<-do.call(rbind,lapply(features[1:5],function(f)data.frame(state=pr$cell_type,feature=f,value=pr[[f]])));p11<-ggplot(mm,aes(feature,state,fill=value))+geom_tile(color="white")+scale_fill_gradient2(low=C["cilia"],mid="white",high=C["repair"],midpoint=0)+labs(title="GSE233145 supportive state profiles",x=NULL,y=NULL,fill="Score")+theme(axis.text.x=element_text(angle=30,hjust=1),panel.grid=element_blank())
saveS(p11,11,10,6)

# S12 original ALI timepoint summary
ali<-rd("02_scRNA/GSE134174/third_stage_results/ALI_full_timecourse_scores.csv");sel<-c("KRT14_KRT17_REPAIR","KRT4_KRT13_TRANSITIONAL_FATE","MULTICILIOGENESIS","MATURE_CILIATED");al<-do.call(rbind,lapply(sel,function(f)data.frame(donor=ali$donor,day=ali$day_numeric,program=f,value=ali[[f]])));p12<-ggplot(al,aes(day,value,group=donor,color=program))+geom_line(alpha=.35)+stat_summary(aes(group=program),fun=mean,geom="line",linewidth=1.1)+scale_color_manual(values=c(KRT14_KRT17_REPAIR=C["repair"],KRT4_KRT13_TRANSITIONAL_FATE=C["transition"],MULTICILIOGENESIS=C["ciliating"],MATURE_CILIATED=C["cilia"]))+labs(title="Original ALI timepoint summaries",x="ALI day",y="Score",color=NULL)+theme(legend.position="bottom")
saveS(p12,12)

# S13 trajectory QC
pt<-rd("12_targeted_manuscript_upgrades/02_ALI_trajectory/slingshot_cell_pseudotime.csv");p13a<-ggplot(pt,aes(pseudotime,fill=donor))+geom_histogram(bins=25,position="identity",alpha=.35)+labs(title="Pseudotime coverage by donor",x="Slingshot pseudotime",y="Cells",fill="Donor")
qc<-read.csv(file.path(repo_root,"figures","main","source_data","Figure6_equal_cell_bin_QC.csv"),check.names=FALSE);p13b<-ggplot(qc,aes(equal_cell_bin,n_cells))+geom_col(fill=C["neutral"],width=.7)+geom_hline(yintercept=150,linetype=2,color=C["repair"])+geom_text(aes(label=paste0(n_donors," donors")),vjust=-.4,size=3)+coord_cartesian(ylim=c(0,175))+labs(title="Equal-cell display bins",subtitle="Prespecified display rule: at least 150 cells and at least 2 donors",x="Pseudotime decile",y="Cells")
saveS(p13a+p13b+plot_annotation(tag_levels="A"),13,11,5.8)

# S14 donor trajectories
bn<-read.csv(file.path(repo_root,"figures","main","source_data","staged_csv","Fig6_ALI_trajectory.csv"),check.names=FALSE);bn<-bn[bn$record_type=="donor_equal_cell_bin",];p14<-ggplot(bn,aes(pseudotime,Repair_Cilia_Imbalance,color=donor))+geom_line(linewidth=.75)+geom_point(size=1.2)+labs(title="Donor-specific Repair-Cilia trajectories",subtitle="Equal-cell display bins; frozen pseudotime",x="Slingshot pseudotime",y="Repair-Cilia score",color="Donor")
saveS(p14,14)

# S15 regulators and KRT4/KRT13 nuance
rg<-do.call(rbind,lapply(c("FOXJ1","MCIDAS","GMNC","KRT4_KRT13_TRANSITIONAL_FATE"),function(f)data.frame(donor=bn$donor,pseudotime=bn$pseudotime,program=f,value=bn[[f]])));p15<-ggplot(rg,aes(pseudotime,value,color=program,group=interaction(program,donor)))+geom_line(alpha=.30)+stat_summary(aes(group=program),fun=mean,geom="line",linewidth=1.1)+scale_color_manual(values=c(FOXJ1=C["cilia"],MCIDAS=C["ciliating"],GMNC=C["neutral"],KRT4_KRT13_TRANSITIONAL_FATE=C["transition"]))+labs(title="Regulator-specific and transitional trajectories",x="Pseudotime",y="Donor-bin score",color=NULL)+theme(legend.position="bottom")
saveS(p15,15)

# S16 maturation derivation
genes<-rd("12_targeted_manuscript_upgrades/03_ALI_maturation_projection/ALI_maturation_signature.csv");cov<-rd("12_targeted_manuscript_upgrades/03_ALI_maturation_projection/platform_coverage.csv");p16a<-ggplot(as.data.frame(table(genes$direction)),aes(Freq,Var1,fill=Var1))+geom_col()+geom_text(aes(label=Freq),hjust=-.2)+scale_fill_manual(values=c(maturation_down=C["repair"],maturation_up=C["maturation"]),guide="none")+coord_cartesian(xlim=c(0,max(table(genes$direction))*1.2))+labs(title="Frozen gene directions",x="Genes",y=NULL)
flow<-data.frame(x=1:4,label=c("Normal ALI\nonly","Donor-aware\npseudotime","53-gene frozen\nscore","Projection to\nPML"));p16b<-ggplot(flow)+geom_rect(aes(xmin=x-.38,xmax=x+.38,ymin=.42,ymax=.70),fill="white",color=C["maturation"],linewidth=.6)+geom_segment(data=data.frame(x=1:3),aes(x=x+.4,xend=x+.58,y=.56,yend=.56),arrow=arrow(length=unit(1.5,"mm")),color=C["neutral"])+geom_text(aes(x=x,y=.56,label=label),size=3)+annotate("text",x=2.5,y=.24,label="No PML labels used in score derivation",fontface="bold",color=C["maturation"],size=3)+coord_cartesian(xlim=c(.45,4.55),ylim=c(.12,.82))+labs(title="Normal-reference derivation")+theme_void(base_family="Arial")+theme(plot.title=element_text(size=11,face="bold"))
saveS(p16a+p16b+plot_annotation(tag_levels="A"),16)

# S17 cross-platform coverage
value_col<-grep("observed|detected",names(cov),value=TRUE)[1];p17<-ggplot(cov,aes(.data[[value_col]],reorder(dataset,.data[[value_col]])))+geom_col(fill=C["maturation"])+geom_text(aes(label=.data[[value_col]]),hjust=-.2)+labs(title="Cross-platform maturation-gene coverage",x="Observed genes",y=NULL)
saveS(p17,17)

# S17 maturation projection sensitivities
proj<-rd("12_targeted_manuscript_upgrades/03_ALI_maturation_projection/PML_maturation_projection_results.csv");proj$feature<-paste(proj$dataset,proj$comparison,sep=" | ");proj$estimate<-proj$estimate_SD;proj$lower<-proj$CI95_lower;proj$upper<-proj$CI95_upper
saveS(forest(proj,"feature",title="Frozen ALI maturation projections",xlabel="Standardized effect",color=C["maturation"]),18,11,7)

# S18 maturation persistence directional
g114<-rd("12_targeted_manuscript_upgrades/03_ALI_maturation_projection/GSE114489_ALI_maturation_scores.csv");g114$outcome<-ifelse(g114$group==1,"Persistent","Regressive")
p18<-boxpts(g114,"outcome","ALI_MATURATION_SCORE","GSE114489 maturation score","ALI maturation score")+annotate("text",x=1.5,y=Inf,label="Adjusted effect −0.582 SD; 95% CI −1.234 to 0.069; FDR 0.0782\nDirectional only; no significance symbol",vjust=1.5,size=3)
saveS(p18,19)

# S19 original NNLS sensitivity
dc<-rd("12_targeted_manuscript_upgrades/01_deconvolution_robustness/deconvolution_method_comparison.csv");nn<-dc[dc$method=="NNLS",];p19<-ggplot(data.frame(model=c("Unadjusted","NNLS-adjusted"),beta=c(nn$severity_beta_unadjusted,nn$severity_beta_adjusted)),aes(model,beta,fill=model))+geom_col(width=.6)+geom_text(aes(label=sprintf("%.3f",beta)),vjust=-.3)+scale_fill_manual(values=c(Unadjusted=C["imbalance"],"NNLS-adjusted"=C["neutral"]),guide="none")+labs(title="Single-method NNLS sensitivity",subtitle=sprintf("Estimated attenuation %.1f%%",nn$attenuation_percent),x=NULL,y="Severity beta per grade")
saveS(p19,20)

# S20 deconvolution benchmark
bm<-rd("12_targeted_manuscript_upgrades/01_deconvolution_robustness/benchmark_method_status.csv");p20a<-ggplot(bm[!is.na(bm$overall_pearson),],aes(overall_pearson,reorder(method,overall_pearson),color=benchmark_status))+geom_point(size=3)+geom_vline(xintercept=.70,linetype=2,color="#B2182B")+scale_color_manual(values=c(FAIL=C["repair"]))+coord_cartesian(xlim=c(0,1))+labs(title="LODO benchmark",x="Overall Pearson correlation",y=NULL,color=NULL)
status<-data.frame(method=bm$method,label=c(sprintf("NNLS_LODO: FAIL (r=%.3f, rho=%.3f, MAE=%.3f)",bm$overall_pearson[bm$method=="NNLS_LODO"],bm$overall_spearman[bm$method=="NNLS_LODO"],bm$overall_MAE[bm$method=="NNLS_LODO"]),sprintf("DWLS: FAIL (r=%.3f, rho=%.3f, MAE=%.3f)",bm$overall_pearson[bm$method=="DWLS"],bm$overall_spearman[bm$method=="DWLS"],bm$overall_MAE[bm$method=="DWLS"]),"MuSiC: METHOD_FAILED"));p20b<-ggplot(status)+geom_text(aes(x=0,y=rev(seq_len(nrow(status))),label=label),hjust=0,size=3.1)+annotate("text",x=0,y=.35,label="Multi-method robustness was not established.",hjust=0,fontface="bold",color=C["repair"],size=3.2)+coord_cartesian(xlim=c(0,1),ylim=c(.1,3.8))+labs(title="Prespecified disposition")+theme_void(base_family="Arial")+theme(plot.title=element_text(size=11,face="bold"))
saveS(p20a+p20b+plot_annotation(tag_levels="A"),21,11,5.5)

# S21 DECAMP batch confounding
de<-rd("12_targeted_manuscript_upgrades/04_DECAMP_validation/GSE300258_frozen_scores.csv");de$group<-factor(ifelse(de$dysplasia==1,"Dysplasia","Non-dysplasia"),levels=c("Non-dysplasia","Dysplasia"));p21a<-boxpts(de,"group","Repair_Cilia_Imbalance","Repair-Cilia Imbalance","Score");p21b<-boxpts(de,"group","ALI_MATURATION_SCORE","ALI maturation score","Score")
banner<-grid::textGrob("Histology was completely confounded with GEO batch in the available metadata.",gp=grid::gpar(col="#B2182B",fontface="bold",fontsize=12))
p21<-wrap_elements(full=banner)/(p21a+p21b+plot_annotation(tag_levels="A"))+plot_layout(heights=c(.08,.92))
saveS(p21,22,11,6)

# S22–S24 exposures
for(z in list(c(23,"PMD","PM2.5 controlled dose response"),c(24,"SMK","Whole smoke donor-paired effects"))){n<-as.integer(z[1]);pref<-z[2];title<-z[3];q<-prov[grepl(paste0("^",pref),prov$result_id),];q$feature<-q$feature;q$lower<-q$lower_CI;q$upper<-q$upper_CI;saveS(forest(q,"feature",title=title,xlabel="Dataset-specific effect"),n,10,7)}

# S25 stretch, upstream, and immune negative results
sq<-prov[grepl("^STR",prov$result_id),];sq$lower<-sq$lower_CI;sq$upper<-sq$upper_CI;p25a<-forest(sq,"feature",title="Mechanical stretch 24-h effects",xlabel="Dataset-specific effect")
stef<-rd("05_stretch/GSE59128_unified/paired_signature_effects.csv");stef<-stef[stef$time_h==24&stef$signature%in%c("YAP_TAZ_MECHANOTRANSDUCTION","PIEZO_MECHANICAL","NOTCH"),];stef$feature<-gsub("_"," ",stef$signature);stef$estimate<-stef$mean_delta_stretch_minus_sham;stef$lower<-stef$CI95_lower;stef$upper<-stef$CI95_upper
p25b<-forest(stef,"feature",title="YAP / PIEZO / NOTCH at 24 h",xlabel="Stretch - sham")
im<-prov[prov$result_id=="IMMUNE01",];im$feature<-"Immune Module 9 enrichment core";im$lower<-im$lower_CI;im$upper<-im$upper_CI;p25c<-forest(im,"feature",title="Exploratory immune lesion-fate analysis",xlabel="Persistent - regressive",color=C["null"])+annotate("text",x=Inf,y=-Inf,label="Non-replicated exploratory result",hjust=1.05,vjust=-.6,size=3)
saveS(p25a/(p25b+p25c)+plot_annotation(tag_levels="A"),25,11,8)
cat("Built supplementary Figures S1-S25.\n")
