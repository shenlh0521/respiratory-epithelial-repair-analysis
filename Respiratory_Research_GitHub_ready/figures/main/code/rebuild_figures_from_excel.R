#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(openxlsx)
  library(ggplot2)
  library(patchwork)
  library(grid)
})

args <- commandArgs(trailingOnly=FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
code_dir <- dirname(normalizePath(script_arg))
base <- dirname(code_dir)
workbook <- file.path(base, "source_data", "Figure_1-8_Source_Data.xlsx")
extract_dir <- file.path(base, "qa", "extracted_from_workbook_for_plot")
out <- normalizePath(file.path(base, "outputs"), mustWork=FALSE)
dir.create(out, recursive=TRUE, showWarnings=FALSE)

C <- c(repair="#D55E00", repair_light="#E69F73", imbalance="#A44200", transition="#CC79A7",
       transition_dark="#8E4B87", cilia="#0072B2", ciliating="#56B4E9", maturation="#005F73",
       secretory="#5F8F7B", persistent="#B54A2B", regressive="#4C78A8", neutral="#666666",
       null="#C9CDD0", stress="#7A5195", oxidative="#E69F00", ink="#1F2529", grid="#E6E8EA")

theme_pub <- function() theme_classic(base_family="Arial", base_size=8) +
  theme(axis.title=element_text(size=9, color=C["ink"]), axis.text=element_text(size=8, color=C["ink"]),
        plot.title=element_text(size=10, face="bold", hjust=0, margin=margin(b=5)),
        legend.text=element_text(size=8), legend.title=element_text(size=8),
        plot.tag=element_text(size=13, face="bold"), plot.tag.position=c(-0.035,1.02),
        panel.grid.major.y=element_line(color=C["grid"], linewidth=.25), panel.grid.minor=element_blank(),
        plot.margin=margin(5,6,5,6))
theme_set(theme_pub())

rd <- function(sheet) read.csv(file.path(extract_dir, paste0(sheet, ".csv")), check.names=FALSE, na.strings=c("", "NA", "NaN"))
num <- function(x) suppressWarnings(as.numeric(x))
label_fdr <- function(x) ifelse(is.na(x), "", ifelse(x == 0, "FDR <1e-300", sprintf("FDR %.2g", x)))

save_plot <- function(p, n, width, height) {
  figure_only <- Sys.getenv("FIGURE_ONLY", "")
  if (nzchar(figure_only) && as.integer(figure_only) != n) return(invisible(NULL))
  base <- file.path(out, paste0("Figure_", n))
  cairo_pdf(paste0(base, ".pdf"), width=width, height=height, family="Arial", onefile=TRUE); print(p); dev.off()
  svg(paste0(base, ".svg"), width=width, height=height, family="Arial", onefile=TRUE); print(p); dev.off()
  png(paste0(base, ".png"), width=width, height=height, units="in", res=600, type="cairo", family="Arial"); print(p); dev.off()
}

forest_plot <- function(d, label, estimate="estimate", lower="lower_CI", upper="upper_CI", color=C["imbalance"], title="", xlabel="Effect (95% CI)", fdr="FDR") {
  d$lab <- d[[label]]; d$est <- num(d[[estimate]]); d$lo <- num(d[[lower]]); d$hi <- num(d[[upper]])
  d$fdrlab <- if (fdr %in% names(d)) label_fdr(num(d[[fdr]])) else ""
  xr <- range(c(d$lo,d$hi),na.rm=TRUE); d$label_x <- d$hi + 0.03*diff(xr)
  d$lab <- factor(d$lab, levels=rev(unique(d$lab)))
  ggplot(d, aes(est, lab)) + geom_vline(xintercept=0, linetype=2, color="#777777", linewidth=.35) +
    geom_errorbarh(aes(xmin=lo, xmax=hi), height=.12, color=color, linewidth=.55) +
    geom_point(color=color, size=2) + geom_text(aes(x=label_x, label=fdrlab), hjust=0, vjust=-.15, size=2.35, family="Arial") +
    labs(x=xlabel, y=NULL, title=title) + scale_x_continuous(expand=expansion(mult=c(.05,.18))) + theme_pub() + theme(panel.grid.major.y=element_blank())
}

box_points <- function(d, group, value, order, palette, title, ylab) {
  d <- d[d[[group]] %in% order & !is.na(d[[value]]),]
  d$grp <- factor(d[[group]], levels=order); d$val <- num(d[[value]])
  ggplot(d, aes(grp, val, color=grp, fill=grp)) +
    geom_boxplot(width=.55, outlier.shape=NA, alpha=.28, linewidth=.45) +
    geom_jitter(width=.09, height=0, size=.65, alpha=.32, stroke=0) +
    scale_color_manual(values=setNames(unname(palette), names(palette)), guide="none") + scale_fill_manual(values=setNames(unname(palette), names(palette)), guide="none") +
    labs(x=NULL, y=ylab, title=title) + theme_pub() + theme(axis.text.x=element_text(angle=20,hjust=1))
}

schematic_panel <- function(lines, edge, title) {
  n <- nrow(lines); lines$y <- rev(seq_len(n));
  ggplot(lines) + geom_rect(aes(xmin=0,xmax=1,ymin=y-.34,ymax=y+.34),fill="#F7F8F9",color="#D1D5D8",linewidth=.35) +
    geom_text(aes(x=.04,y=y,label=dataset),hjust=0,fontface="bold",family="Arial",size=2.45,color=edge) +
    geom_text(aes(x=.36,y=y,label=paste(role,sample_summary,sep="\n")),hjust=0,family="Arial",size=2.25,lineheight=.92) +
    coord_cartesian(xlim=c(0,1),ylim=c(.4,n+.6),clip="off") + labs(title=title) + theme_void(base_family="Arial") +
    theme(plot.title=element_text(size=10,face="bold",hjust=0,margin=margin(b=5)),plot.margin=margin(5,6,5,6))
}

# Figure 1
d <- rd("Fig1_design")
p1a <- schematic_panel(d[d$record_type=="dataset" & d$domain=="PML",], C["imbalance"], "Human PML cohorts")
p1b <- schematic_panel(d[d$record_type=="dataset" & d$domain=="Single-cell",], C["maturation"], "Single-cell and normal differentiation")
p1c <- schematic_panel(d[d$record_type=="dataset" & d$domain=="Perturbation",], C["stress"], "Controlled airway perturbation")
fw <- d[d$record_type=="framework",][1:5,]; fw$x <- 1:5; fw$short<-c("Airway\nperturbation","Repair-\ndominant","Transitional","Multi-\nciliogenesis","Mature\nciliated")
p1d <- ggplot(fw) + geom_rect(aes(xmin=x-.46,xmax=x+.46,ymin=.44,ymax=.76),fill="white",color=c(C["stress"],C["repair"],C["transition"],C["ciliating"],C["cilia"]),linewidth=.65) +
  geom_segment(data=data.frame(x=1:4),aes(x=x+.47,xend=x+.53,y=.60,yend=.60),arrow=arrow(length=unit(1.6,"mm")),linewidth=.38,color="#777777") +
  geom_text(aes(x=x,y=.60,label=short),size=2.48,family="Arial",lineheight=.88) +
  annotate("text",x=3,y=.30,label="PML: association with a repair-high, less-mature shift",family="Arial",fontface="bold",size=2.82,color=C["maturation"]) +
  coord_cartesian(xlim=c(.48,5.52),ylim=c(.14,.86),clip="off") + labs(title="Evidence framework") + theme_void(base_family="Arial") + theme(plot.title=element_text(size=10,face="bold"))
p <- (p1a+p1b)/(p1c+p1d) + plot_annotation(tag_levels="A") & theme(plot.tag=element_text(size=13,face="bold"),plot.tag.position=c(-.035,1.02))
save_plot(p,1,7.09,6.0)

# Figure 2
d <- rd("Fig2_discovery"); s <- d[d$record_type=="sample_score",]; e <- d[d$record_type=="effect",]
order <- c("Normal","Secretory","Inflammatory","Proliferative"); pal <- c(Normal=unname(C["neutral"]),Secretory=unname(C["secretory"]),Inflammatory=unname(C["transition"]),Proliferative=unname(C["repair"]))
mk <- function(feat,title) box_points(s[s$feature==feat,],"molecular_subtype","value",order,pal,title,"Score")
p2a<-mk("CILIA_CONSENSUS","Cilia consensus");p2b<-mk("ABNORMAL_REPAIR","Abnormal repair");p2c<-mk("Repair_Cilia_Imbalance","Repair-Cilia imbalance")+theme(axis.text.x=element_text(size=7.4,angle=20,hjust=1))
e$color <- ifelse(e$feature=="FOXJ1",C["cilia"],ifelse(e$feature=="Repair_Cilia_Imbalance",C["imbalance"],C["repair"]))
e$est<-num(e$estimate);e$lo<-num(e$lower_CI);e$hi<-num(e$upper_CI);e$feature<-factor(e$feature,levels=rev(e$feature));e$lab<-label_fdr(num(e$FDR))
p2d<-ggplot(e,aes(est,feature,color=color))+geom_vline(xintercept=0,linetype=2,color="#777",linewidth=.35)+geom_errorbarh(aes(xmin=lo,xmax=hi),height=.12,linewidth=.55)+geom_point(size=2)+geom_text(aes(x=hi,label=lab),hjust=1,vjust=-.75,size=2.55,family="Arial")+scale_color_identity()+labs(title="Proliferative vs normal",x="Adjusted effect (95% CI)",y=NULL)+theme_pub()+theme(panel.grid.major.y=element_blank())
p <- (p2a+p2b)/(p2c+p2d)+plot_annotation(tag_levels="A")
save_plot(p,2,7.09,5.5)

# Figure 3
d<-rd("Fig3_severity"); eff<-d[d$record_type=="effect",]; pts<-d[d$record_type=="sample_score",]
p3a<-forest_plot(eff,"dataset",color=C["imbalance"],title="Severity effect across cohorts",xlabel="Repair-Cilia change per grade (95% CI)")+labs(tag="A")+theme(plot.tag.position=c(.01,.99),plot.margin=margin(8,6,5,10))
trend <- function(ds,title) {
  x<-pts[pts$dataset==ds & !is.na(pts$grade_ordinal),]; x$grade_ordinal<-num(x$grade_ordinal);x$value<-num(x$value)
  agg<-aggregate(value~patient+grade_ordinal,x,mean)
  sm_mean<-aggregate(value~grade_ordinal,agg,mean); sm_sd<-aggregate(value~grade_ordinal,agg,sd); sm_n<-aggregate(value~grade_ordinal,agg,length)
  sm<-data.frame(grade_ordinal=sm_mean$grade_ordinal,mean=sm_mean$value,lo=sm_mean$value-1.96*sm_sd$value/sqrt(sm_n$value),hi=sm_mean$value+1.96*sm_sd$value/sqrt(sm_n$value))
  ggplot(agg,aes(grade_ordinal,value))+geom_jitter(width=.06,height=0,size=.9,color=C["imbalance"],alpha=.38)+geom_line(data=sm,aes(grade_ordinal,mean),inherit.aes=FALSE,color=C["imbalance"],linewidth=.7)+geom_point(data=sm,aes(grade_ordinal,mean),inherit.aes=FALSE,color=C["imbalance"],size=1.6)+geom_errorbar(data=sm,aes(x=grade_ordinal,ymin=lo,ymax=hi),inherit.aes=FALSE,width=.08,color=C["imbalance"],linewidth=.45)+labs(title=title,x="Histologic grade",y="Repair-Cilia score")+theme_pub()
}
p3b<-trend("GSE109743","GSE109743")+labs(tag="B");p3c<-trend("GSE320381","GSE320381")+labs(tag="C");p3d<-trend("GSE33479","GSE33479")+labs(tag="D")
p <- p3a/(p3b+p3c+p3d)+plot_layout(heights=c(.72,1.28))
save_plot(p,3,7.09,5.7)

# Figure 4
d<-rd("Fig4_persistence"); pdat<-d[d$record_type=="sample_score",]
p4a<-ggplot()+annotate("rect",xmin=.02,xmax=.39,ymin=.43,ymax=.76,fill="white",color=C["neutral"],linewidth=.55)+annotate("rect",xmin=.61,xmax=.98,ymin=.43,ymax=.76,fill="white",color=C["neutral"],linewidth=.55)+annotate("segment",x=.40,xend=.60,y=.595,yend=.595,arrow=arrow(length=unit(1.6,"mm")),color="#777",linewidth=.4)+annotate("text",x=.205,y=.595,label="Baseline\nlesion-site biopsy",family="Arial",size=3.0,lineheight=.9)+annotate("text",x=.795,y=.595,label="Follow-up\npathology",family="Arial",size=3.0,lineheight=.9)+annotate("text",x=.50,y=.27,label="Persistent or regressive dysplasia",family="Arial",fontface="bold",size=2.85)+coord_cartesian(xlim=c(0,1),ylim=c(.1,.85))+labs(title="Baseline-to-follow-up design")+theme_void(base_family="Arial")+theme(plot.title=element_text(size=10,face="bold"))
p4b<-box_points(pdat,"outcome","value",c("Regressive","Persistent"),c(Regressive=C["regressive"],Persistent=C["persistent"]),"Baseline Repair-Cilia score","Repair-Cilia score")
di<-d[d$record_type=="adjusted_difference",][1,]; od<-d[d$record_type=="adjusted_or",][1,]
ce<-data.frame(label=c("Difference (SD)","Log odds ratio"),est=c(num(di$estimate),log(num(od$estimate))),lo=c(num(di$CI95_lower),log(num(od$CI95_lower))),hi=c(num(di$CI95_upper),log(num(od$CI95_upper))),fdr=c(num(di$FDR),num(od$FDR)))
p4c<-forest_plot(ce,"label",estimate="est",lower="lo",upper="hi",color=C["persistent"],title="Adjusted lesion-outcome estimates",xlabel="Estimate (95% CI; null=0)",fdr="fdr")+scale_x_continuous(expand=expansion(mult=c(.05,.32)))+theme(plot.margin=margin(5,16,5,6))
fld<-d[d$record_type=="field_null",][1,]; fd<-data.frame(dataset="GSE79210",estimate=num(fld$estimate),lower=num(fld$CI95_lower),upper=num(fld$CI95_upper),FDR=num(fld$FDR))
p4d<-forest_plot(fd,"dataset",lower="lower",upper="upper",color=C["null"],title="Airway field-brushing sensitivity",xlabel="Persistent − regressive (95% CI)")+scale_x_continuous(expand=expansion(mult=c(.05,.32)))+labs(subtitle="Different sampling compartment; null result")+theme(plot.subtitle=element_text(size=7.5,color=C["neutral"],margin=margin(b=4)),plot.margin=margin(5,16,5,6))
p<-(p4a+p4b)/(p4c+p4d)+plot_annotation(tag_levels="A")
save_plot(p,4,7.09,5.2)

# Figure 5
d<-rd("Fig5_scRNA"); st<-d[d$record_type=="state_score",]; comp<-d[d$record_type=="composition_effect",]
state_group<-function(x) ifelse(grepl("ciliated",x,ignore.case=TRUE),"Mature ciliated",ifelse(grepl("ciliating",x,ignore.case=TRUE),"Ciliating",ifelse(grepl("mucus|secretory",x,ignore.case=TRUE),"Secretory",ifelse(grepl("Differentiating|KRT8",x),"Transitional","Basal/repair"))))
st$group<-state_group(st$subcluster_ident); spal<-c("Basal/repair"=unname(C["repair"]),Transitional=unname(C["transition"]),Secretory=unname(C["secretory"]),Ciliating=unname(C["ciliating"]),"Mature ciliated"=unname(C["cilia"]))
p5a<-ggplot(st,aes(num(dim1),num(dim2),color=group))+geom_point(size=2.4,alpha=.85)+scale_color_manual(values=spal)+labs(title="Epithelial state map",x="State PCA 1",y="State PCA 2",color=NULL)+theme_pub()+theme(legend.position="bottom")
p5b<-ggplot(st,aes(num(CILIA_CONSENSUS),num(Repair_Cilia_Imbalance),color=group))+geom_hline(yintercept=0,color="#aaa",linewidth=.3)+geom_vline(xintercept=0,color="#aaa",linewidth=.3)+geom_point(size=2.3)+scale_color_manual(values=spal,guide="none")+labs(title="Repair-high versus ciliated states",x="CILIA_CONSENSUS",y="Repair-Cilia score")+theme_pub()+theme(panel.grid=element_blank())
features<-c("KRT14_KRT17_REPAIR","ABNORMAL_REPAIR","MULTICILIOGENESIS","MATURE_CILIATED"); st$state_label<-make.unique(as.character(st$subcluster_ident)); mm<-do.call(rbind,lapply(features,function(f)data.frame(state=st$state_label,feature=f,value=scale(num(st[[f]]))[,1]))); mm$state<-factor(mm$state,levels=st$state_label[order(num(st$state_order))])
p5c<-ggplot(mm,aes(feature,state,fill=value))+geom_tile(color="white",linewidth=.25)+scale_fill_gradient2(low=C["cilia"],mid="white",high=C["repair"],midpoint=0,limits=c(-2,2),oob=scales::squish)+scale_x_discrete(labels=c(KRT14_KRT17_REPAIR="KRT14/17\nrepair",ABNORMAL_REPAIR="Abnormal\nrepair",MULTICILIOGENESIS="Multi-\nciliogenesis",MATURE_CILIATED="Mature\nciliated"))+labs(title="Selected state programs",x=NULL,y=NULL,fill="State z")+theme_pub()+theme(panel.grid=element_blank(),axis.text.y=element_text(size=6.9),axis.text.x=element_text(size=7.3),legend.text=element_text(size=7.2),plot.margin=margin(5,2,5,8))
comp$est<-num(comp$estimate_percentage_points);comp$lo<-num(comp$CI95_lower_percentage_points);comp$hi<-num(comp$CI95_upper_percentage_points);comp$state<-factor(comp$state,levels=comp$state)
p5d<-ggplot(comp,aes(est,state))+geom_vline(xintercept=0,linetype=2,color="#777",linewidth=.35)+geom_errorbarh(aes(xmin=lo,xmax=hi),height=.12,color=C["neutral"],linewidth=.5)+geom_point(color=C["neutral"],size=2)+labs(title="Donor composition",x="Heavy − never smoker\n(percentage points)",y=NULL)+theme_pub()+theme(panel.grid.major.y=element_blank(),axis.text.y=element_text(size=6.4),plot.title=element_text(size=9.3,face="bold"),plot.margin=margin(5,7,5,4))+annotate("text",x=Inf,y=-Inf,label="Not FDR-significant",hjust=1.05,vjust=-.7,size=2.2,family="Arial",color=C["neutral"])
p<-(p5a+p5b)/(p5c+p5d+plot_layout(widths=c(1.18,.82)))+plot_annotation(tag_levels="A")
save_plot(p,5,7.09,5.5)

# Figure 6
d<-rd("Fig6_ALI_trajectory"); emb<-d[d$record_type=="cell_embedding",]; bins<-d[d$record_type=="donor_equal_cell_bin",]; tests<-d[d$record_type=="trajectory_test",]
cluster_order<-c("basal.colonies","basal.subconfluent","basal.confluent","p.basal","d.basal","ciliating.early","ciliating.late","ciliated"); cpal<-setNames(c(C["repair"],C["repair_light"],"#C47B52",C["transition_dark"],C["transition"],C["ciliating"],"#2B8CBE",C["cilia"]),cluster_order)
curve_dat<-data.frame(dim1=num(emb$dim1),dim2=num(emb$dim2),ptbin=cut(num(emb$pseudotime),breaks=18)); curve<-aggregate(cbind(dim1,dim2)~ptbin,data=curve_dat,median)
p6a<-ggplot(emb,aes(num(dim1),num(dim2),color=cluster))+geom_point(size=.35,alpha=.28)+geom_path(data=curve,aes(dim1,dim2),inherit.aes=FALSE,color=C["ink"],linewidth=.75)+geom_segment(data=curve[nrow(curve)-2,],aes(x=dim1,y=dim2,xend=curve$dim1[nrow(curve)],yend=curve$dim2[nrow(curve)]),inherit.aes=FALSE,arrow=arrow(length=unit(1.6,"mm")),color=C["ink"],linewidth=.75)+scale_color_manual(values=cpal,guide="none")+labs(title="Slingshot lineage",x="PCA 1",y="PCA 2")+theme_pub()+theme(panel.grid=element_blank())+annotate("text",x=-Inf,y=-Inf,label="3 donors",hjust=-.05,vjust=-.7,size=2.2,family="Arial")
bins$pseudotime<-num(bins$pseudotime);bins$Repair_Cilia_Imbalance<-num(bins$Repair_Cilia_Imbalance)
ov<-aggregate(cbind(pseudotime,Repair_Cilia_Imbalance)~pseudotime_bin,bins,mean);tt<-tests[tests$program=="Repair_Cilia_Imbalance",][1,]; model_line<-data.frame(pseudotime=range(bins$pseudotime)); model_line$score<-mean(ov$Repair_Cilia_Imbalance[ov$pseudotime_bin<=2])+num(tt$slope_per_unit_pseudotime)*(model_line$pseudotime-min(model_line$pseudotime))
p6b<-ggplot(bins,aes(pseudotime,Repair_Cilia_Imbalance,group=donor))+geom_line(color=C["imbalance"],alpha=.30,linewidth=.45)+geom_point(color=C["imbalance"],alpha=.38,size=.7)+geom_line(data=model_line,aes(pseudotime,score,group=1),inherit.aes=FALSE,color=C["imbalance"],linewidth=1.2)+labs(title="Repair-Cilia along pseudotime",subtitle="Donor x equal-cell pseudotime bins (all bins: >=150 cells; 3 donors)",x="Slingshot pseudotime",y="Repair-Cilia score")+theme_pub()+theme(plot.subtitle=element_text(size=7.4,color=C["neutral"],margin=margin(b=3)))+annotate("text",x=Inf,y=Inf,label=sprintf("β = %.3f\n95%% CI %.3f to %.3f\nFDR = %.5f\nn = 3 donors",num(tt$slope_per_unit_pseudotime),num(tt$CI95_lower),num(tt$CI95_upper),num(tt$FDR)),hjust=1.03,vjust=1.06,size=2.55,lineheight=.9,family="Arial")
normtraj<-function(features,labels,cols,title){z<-do.call(rbind,lapply(seq_along(features),function(i){f<-features[i];o<-aggregate(cbind(pseudotime=num(bins$pseudotime),value=num(bins[[f]]))~pseudotime_bin,bins,mean);o$value<-(o$value-mean(o$value,na.rm=TRUE))/sd(o$value,na.rm=TRUE);o$program<-labels[i];o}));ggplot(z,aes(pseudotime,value,color=program))+geom_hline(yintercept=0,color="#aaa",linewidth=.25)+geom_line(linewidth=.95)+scale_color_manual(values=setNames(cols,labels))+labs(title=title,x="Slingshot pseudotime",y="Normalized score",color=NULL)+theme_pub()+theme(legend.position="bottom")}
p6c<-normtraj(c("KRT14_KRT17_REPAIR","KRT4_KRT13_TRANSITIONAL_FATE","MULTICILIOGENESIS","MATURE_CILIATED"),c("KRT14/17 repair","KRT4/13 transitional","Multiciliogenesis","Mature ciliated"),c(C["repair"],C["transition"],C["ciliating"],C["cilia"]),"Program trajectories")
p6d<-normtraj(c("FOXJ1","MCIDAS","GMNC"),c("FOXJ1","MCIDAS","GMNC"),c(C["cilia"],C["ciliating"],C["neutral"]),"Canonical regulators differ")+labs(subtitle=sprintf("FOXJ1 directional (FDR %.3g)\nMCIDAS weak (FDR %.3g); GMNC no clear rise (FDR %.3g)",num(tests$FDR[tests$program=="FOXJ1"]),num(tests$FDR[tests$program=="MCIDAS"]),num(tests$FDR[tests$program=="GMNC"])))+theme(plot.subtitle=element_text(size=7.7,lineheight=.9,color=C["neutral"],margin=margin(b=4)))
p<-(p6a+p6b)/(p6c+p6d)+plot_annotation(tag_levels="A") & theme(legend.text=element_text(size=6.4),legend.key.width=unit(3,"mm"))
save_plot(p,6,7.09,5.5)

# Figure 7
d<-rd("Fig7_maturation"); flow<-d[d$record_type=="workflow",]; proj<-d[d$record_type=="projection_effect",]; ss<-d[d$record_type=="sample_score",]; corr<-d[d$record_type=="correlation",]
flow$x<-num(flow$step); flow$short<-c("Normal ALI\nreference","Donor-aware\npseudo-\ntime","Frozen\n53-gene\nscore","PML\nprojection")
p7a<-ggplot(flow)+
  geom_rect(aes(xmin=x-.36,xmax=x+.36,ymin=.38,ymax=.76),fill="#F8FAFB",color=C["maturation"],linewidth=.62)+
  geom_segment(data=data.frame(x=1:3),aes(x=x+.39,xend=x+.61,y=.57,yend=.57),arrow=arrow(length=unit(1.45,"mm")),linewidth=.42,color=C["neutral"])+
  geom_text(aes(x=x,y=.57,label=short),size=2.36,family="Arial",lineheight=.84)+
  annotate("text",x=2.5,y=.20,label="Derived in normal ALI only; no PML labels used.",family="Arial",fontface="bold",size=2.52,color=C["maturation"])+
  coord_cartesian(xlim=c(.42,4.58),ylim=c(.10,.84),clip="off")+
  labs(title="Frozen normal-reference score",tag="A")+
  theme_void(base_family="Arial")+
  theme(plot.title=element_text(size=10,face="bold",margin=margin(b=5)),plot.tag.position=c(.01,.99),plot.tag=element_text(size=12.5,face="bold"),plot.margin=margin(7,6,5,8))
sev<-proj[proj$comparison=="per_one_histology_grade",]
p7b<-forest_plot(sev,"dataset",estimate="estimate_SD",lower="CI95_lower",upper="CI95_upper",color=C["maturation"],title="Maturation decreases with severity",xlabel="Maturation score per grade (95% CI)")+
  labs(tag="B")+theme(plot.tag.position=c(.01,.99),plot.tag=element_text(size=12.5,face="bold"),plot.margin=margin(7,8,5,8))
sub<-ss[ss$molecular_subtype %in% c("Normal","Proliferative"),]
pal7<-setNames(c(unname(C["neutral"]),unname(C["repair"])),c("Normal","Proliferative"))
ee<-proj[proj$dataset=="GSE109743" & proj$comparison=="Proliferative_vs_Normal",][1,]
stat7<-sprintf("Proliferative vs Normal: %.2f SD (95%% CI %.2f to %.2f); FDR %.2g",num(ee$estimate_SD),num(ee$CI95_lower),num(ee$CI95_upper),num(ee$FDR))
p7c<-box_points(sub,"molecular_subtype","ALI_MATURATION_SCORE",c("Normal","Proliferative"),pal7,"Subtype contrast in GSE109743","ALI maturation score")+
  labs(subtitle=stat7,tag="C")+
  theme(axis.text.x=element_text(angle=0,hjust=.5),plot.subtitle=element_text(size=7.1,color=C["neutral"],margin=margin(b=4)),plot.tag.position=c(.01,.99),plot.tag=element_text(size=12.5,face="bold"),plot.margin=margin(7,7,6,8))
p7d<-forest_plot(corr,"dataset",estimate="correlation",lower="CI95_lower",upper="CI95_upper",color=C["transition_dark"],title="Inverse score correlation across cohorts",xlabel="Spearman correlation (95% CI)")+
  labs(subtitle="Related, not independent constructs; cluster-bootstrap 95% CIs",caption="Numerical underflow is displayed as FDR <1e-300.",tag="D")+
  theme(plot.subtitle=element_text(size=7.1,color=C["neutral"],margin=margin(b=4)),plot.caption=element_text(size=6.5,hjust=0,color=C["neutral"],margin=margin(t=5)),plot.tag.position=c(.01,.99),plot.tag=element_text(size=12.5,face="bold"),plot.margin=margin(7,8,6,8))
p<-(p7a+p7b)/(p7c+p7d)+plot_layout(widths=c(.97,1.03),heights=c(.94,1.06))
save_plot(p,7,7.09,5.55)

# Figure 8
d<-rd("Fig8_exposure"); ex<-d[d$record_type=="exposure_effect",]; model<-d[d$record_type=="framework",]
rows<-c("OXIDATIVE_STRESS","DNA_DAMAGE_RESPONSE","ABNORMAL_REPAIR","KRT14_KRT17_REPAIR","CILIA_CONSENSUS","MULTICILIOGENESIS","MATURE_CILIATED","Repair_Cilia_Imbalance");cols<-c("PM2.5","Whole smoke","Stretch");ex$feature<-factor(ex$feature,levels=rev(rows));ex$exposure<-factor(ex$exposure,levels=cols);ex$display_effect<-num(ex$display_effect)
p8a<-ggplot(ex,aes(exposure,feature,fill=display_effect,alpha=evidence_strength,color=evidence_strength))+geom_tile(linewidth=.7)+geom_text(aes(label=ifelse(evidence_strength=="FDR<0.05","●","")),size=2.2,color=C["ink"],show.legend=FALSE)+scale_fill_gradient2(low=C["cilia"],mid="#F4F4F4",high=C["repair"],midpoint=0,limits=c(-1,1),na.value="#F7F7F7",name="Direction")+scale_alpha_manual(values=c("FDR<0.05"=.95,"P<0.05"=.60,"not significant"=.25,"not assessed"=.12),name="Support",guide=guide_legend(override.aes=list(fill="#777777",color="white")))+scale_color_manual(values=c("FDR<0.05"=C["ink"],"P<0.05"="white","not significant"="white","not assessed"="white"),guide="none")+scale_y_discrete(labels=function(x)gsub("_"," ",x))+labs(title="Controlled perturbations",x=NULL,y=NULL,caption="Effects normalized within exposure for display; ● = FDR<0.05")+theme_pub()+theme(panel.grid=element_blank(),axis.text.x=element_text(angle=15,hjust=1),plot.caption=element_text(size=7.1,hjust=0,margin=margin(t=5)))
model$y<-rev(seq_len(nrow(model)))
model$short<-c("Different airway insults","Distinct proximal responses","Repair-ciliogenesis\nbalance disturbed","Repair-dominant,\nless-mature state")
p8b<-ggplot(model)+geom_rect(aes(xmin=.03,xmax=.97,ymin=y-.34,ymax=y+.34),fill="white",color=c(C["stress"],C["neutral"],C["transition"],C["maturation"]),linewidth=.65)+geom_segment(data=data.frame(y=4:2),aes(x=.50,xend=.50,y=y-.36,yend=y-.64),arrow=arrow(length=unit(1.6,"mm")),color="#777",linewidth=.4)+geom_text(aes(x=.50,y=y,label=short),family="Arial",size=2.72,lineheight=.88)+annotate("text",x=.50,y=.28,label="Association with PML severity and\nlesion persistence (not causal)",family="Arial",size=2.32,lineheight=.9,color=C["neutral"])+coord_cartesian(xlim=c(0,1),ylim=c(.03,4.55),clip="off")+labs(title="Partial convergence")+theme_void(base_family="Arial")+theme(plot.title=element_text(size=10,face="bold",margin=margin(l=16)),plot.margin=margin(5,4,5,5))
p<-p8a+p8b+plot_layout(widths=c(1.18,.82))+plot_annotation(tag_levels="A") &
  theme(plot.tag.position=c(.018,.985),plot.tag=element_text(size=13,face="bold"))
save_plot(p,8,7.09,3.8)

cat("Rebuilt Figures 1–8 from", workbook, "\n")
