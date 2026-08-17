#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

set.seed(20260811)
project <- Sys.getenv("PROJECT_ROOT", unset = normalizePath(".", mustWork = FALSE))
root <- file.path(project, "02_scRNA/SECOND_VALIDATION/GSE233145")
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

welch_effect <- function(x1, x0, label1 = "COPD", label0 = "nonCLD") {
  x1 <- x1[is.finite(x1)]; x0 <- x0[is.finite(x0)]
  n1 <- length(x1); n0 <- length(x0); est <- mean(x1) - mean(x0)
  v1 <- var(x1); v0 <- var(x0); se <- sqrt(v1/n1 + v0/n0)
  df <- (v1/n1 + v0/n0)^2 / ((v1/n1)^2/(n1-1) + (v0/n0)^2/(n0-1))
  stat <- est/se; p <- 2*pt(abs(stat), df, lower.tail=FALSE); crit <- qt(.975, df)
  data.frame(group_1=label1, group_0=label0, n_group_1=n1, n_group_0=n0,
             estimate=est, std_error=se, CI95_lower=est-crit*se, CI95_upper=est+crit*se,
             statistic=stat, df=df, p_value=p)
}

meta_all <- fread(cmd = sprintf("gzip -cd %s", shQuote(file.path(root, "data/GSE233145_cells_metadata.txt.gz"))))
meta <- copy(meta_all)
meta <- meta[!is.na(cell_type) & cell_type != "nan"]
states <- sort(unique(meta$cell_type))
donor_time <- unique(meta[, .(patient, health_state, time_point)])
grid <- donor_time[, .(cell_type = states), by=.(patient, health_state, time_point)]
counts <- meta[, .(n_cells=.N), by=.(patient, health_state, time_point, cell_type)]
comp <- merge(grid, counts, all.x=TRUE, by=c("patient","health_state","time_point","cell_type"))
comp[is.na(n_cells), n_cells := 0L]
comp[, total_cells := sum(n_cells), by=.(patient,time_point)]
comp[, proportion := n_cells/total_cells]
basal <- comp[grepl("^Basal", cell_type), .(n_cells=sum(n_cells), total_cells=unique(total_cells)),
              by=.(patient,health_state,time_point)][, `:=`(cell_type="Basal_combined", proportion=n_cells/total_cells)]
comp <- rbind(comp, basal, fill=TRUE)
fwrite(comp, file.path(out, "donor_time_cell_state_proportions.csv"))

focus_states <- c("Basal_combined", "Ciliated", "Transitional Ciliated", "Secretory", "Suprabasal")
comp_effects <- rbindlist(unlist(lapply(c("day 0", "day 28"), function(day) lapply(focus_states, function(state) {
  z <- comp[time_point == day & cell_type == state]
  ans <- welch_effect(z[health_state == "COPD"]$proportion, z[health_state == "Donor"]$proportion)
  ans$time_point <- day; ans$cell_type <- state
  ans$estimate_percentage_points <- 100*ans$estimate
  ans$CI95_lower_percentage_points <- 100*ans$CI95_lower
  ans$CI95_upper_percentage_points <- 100*ans$CI95_upper
  ans
})), recursive=FALSE), fill=TRUE)
comp_effects[, FDR := p.adjust(p_value, "BH"), by=time_point]
comp_effects[, analysis_method := "Welch donor-level contrast; n=2 COPD-IV vs n=2 non-CLD"]
fwrite(comp_effects, file.path(out, "cell_state_composition_effects.csv"))

pb <- fread(file.path(root, "derived/target_donor_state_pseudobulk.csv"))
wide <- dcast(pb, patient + health_state + time_point + cell_type + source_object + n_cells + value_scale ~ gene,
              value.var="value")
target_genes <- intersect(unique(unlist(gene_sets)), names(wide))
for (g in target_genes) wide[[paste0(g,"_Z")]] <- as.numeric(scale(wide[[g]]))
for (nm in names(gene_sets)) {
  use <- paste0(intersect(gene_sets[[nm]], names(wide)), "_Z")
  wide[[nm]] <- if(length(use)) rowMeans(wide[, ..use], na.rm=TRUE) else NA_real_
}
wide[, ABNORMAL_REPAIR_Z := as.numeric(scale(ABNORMAL_REPAIR))]
wide[, CILIA_CONSENSUS_Z := as.numeric(scale(CILIA_CONSENSUS))]
wide[, Repair_Cilia_Imbalance := ABNORMAL_REPAIR_Z - CILIA_CONSENSUS_Z]
score_cols <- c(names(gene_sets), "Repair_Cilia_Imbalance")
fwrite(wide[, c("patient","health_state","time_point","cell_type","n_cells","value_scale",score_cols,
                intersect(c("KRT14","KRT17","FOXJ1","MCIDAS","GMNC"),names(wide))), with=FALSE],
       file.path(out, "donor_state_signature_scores.csv"))

coverage <- rbindlist(lapply(names(gene_sets), function(nm) data.frame(
  signature=nm, n_defined=length(gene_sets[[nm]]), n_detected=length(intersect(gene_sets[[nm]],names(wide))),
  detected_genes=paste(intersect(gene_sets[[nm]],names(wide)),collapse=";"),
  missing_genes=paste(setdiff(gene_sets[[nm]],names(wide)),collapse=";"))))
fwrite(coverage,file.path(out,"signature_coverage.csv"))

profiles <- wide[cell_type %in% c("Basal_1","Basal_2","Basal_combined","Suprabasal","Secretory","Transitional Ciliated","Ciliated"),
                 c(lapply(.SD,mean,na.rm=TRUE),list(n_donor_states=.N)), by=.(time_point,cell_type), .SDcols=score_cols]
fwrite(profiles,file.path(out,"cell_state_signature_profiles.csv"))

tests <- rbindlist(lapply(list(
  c("day 0","Basal_combined"), c("day 0","All_epithelial"),
  c("day 28","Basal_combined"), c("day 28","Transitional Ciliated"),
  c("day 28","Ciliated"), c("day 28","All_epithelial")
), function(spec) {
  z <- wide[time_point==spec[1] & cell_type==spec[2]]
  rbindlist(lapply(score_cols,function(nm){
    ans <- welch_effect(z[health_state=="COPD"][[nm]],z[health_state=="Donor"][[nm]])
    ans$time_point <- spec[1]; ans$cell_type <- spec[2]; ans$signature <- nm; ans
  }))
}))
tests[, FDR := p.adjust(p_value,"BH")]
tests[, analysis_method := "donor-state expression summary; Welch contrast; n=2/group"]
fwrite(tests,file.path(out,"within_state_signature_effects.csv"))

plot_comp <- comp_effects[time_point=="day 28"]
plot_comp[, label := factor(cell_type,levels=rev(cell_type))]
p1 <- ggplot(plot_comp,aes(estimate_percentage_points,label))+
  geom_vline(xintercept=0,linetype=2,colour="grey50")+
  geom_errorbarh(aes(xmin=CI95_lower_percentage_points,xmax=CI95_upper_percentage_points),height=.18)+
  geom_point(shape=21,size=3,fill="#4C78A8")+
  labs(title="Day 28 cell-state redistribution",x="COPD-IV minus non-CLD (percentage points; 95% CI)",y=NULL)+
  theme_minimal(base_size=11)
ggsave(file.path(out,"figures/GSE233145_composition_validation.png"),p1,width=8,height=4.8,dpi=240)

plot_sig <- tests[time_point=="day 28" & cell_type=="All_epithelial"]
plot_sig[, label := factor(signature,levels=rev(signature))]
p2 <- ggplot(plot_sig,aes(estimate,label))+
  geom_vline(xintercept=0,linetype=2,colour="grey50")+
  geom_errorbarh(aes(xmin=CI95_lower,xmax=CI95_upper),height=.18)+
  geom_point(shape=21,size=3,fill="#D55E00")+
  labs(title="Day 28 epithelial signature effects",x="COPD-IV minus non-CLD (95% CI)",y=NULL)+
  theme_minimal(base_size=11)
ggsave(file.path(out,"figures/GSE233145_signature_validation.png"),p2,width=8,height=5.2,dpi=240)

profile_focus <- profiles[time_point=="day 28" & cell_type %in% c("Basal_combined","Suprabasal","Transitional Ciliated","Ciliated")]
profile_long <- melt(profile_focus,
  id.vars=c("time_point","cell_type","n_donor_states"),
  measure.vars=c("ABNORMAL_REPAIR","KRT14_KRT17_REPAIR","CILIA_CONSENSUS","MULTICILIOGENESIS","MATURE_CILIATED","Repair_Cilia_Imbalance"),
  variable.name="signature",value.name="mean_score")
profile_long[, cell_type:=factor(cell_type,levels=c("Basal_combined","Suprabasal","Transitional Ciliated","Ciliated"))]
profile_long[, signature:=factor(signature,levels=c("ABNORMAL_REPAIR","KRT14_KRT17_REPAIR","CILIA_CONSENSUS","MULTICILIOGENESIS","MATURE_CILIATED","Repair_Cilia_Imbalance"))]
p3 <- ggplot(profile_long,aes(signature,cell_type,fill=mean_score))+
  geom_tile(colour="white",linewidth=.7)+
  geom_text(aes(label=sprintf("%.2f",mean_score)),size=3.2)+
  scale_fill_gradient2(low="#3B6FB6",mid="white",high="#D94B4B",midpoint=0)+
  labs(title="GSE233145 day-28 epithelial state localization",
       subtitle="Mean donor-state scores; fixed project signatures",
       x=NULL,y=NULL,fill="Score")+
  theme_minimal(base_size=11)+theme(axis.text.x=element_text(angle=40,hjust=1),panel.grid=element_blank())
ggsave(file.path(out,"figures/GSE233145_state_profile_validation.png"),p3,width=8.8,height=4.7,dpi=240)

key_comp <- comp_effects[time_point=="day 28" & cell_type %in% c("Basal_combined","Ciliated","Transitional Ciliated")]
key_sig <- tests[cell_type %in% c("Basal_combined","Transitional Ciliated","Ciliated","All_epithelial") &
                 signature %in% c("ABNORMAL_REPAIR","KRT14_KRT17_REPAIR","CILIA_CONSENSUS","MULTICILIOGENESIS","MATURE_CILIATED","Repair_Cilia_Imbalance")]
lines <- c("# Second-stage scRNA validation: GSE233145", "",
  sprintf("Author metadata contain %s total cells; %s cells have a retained author cell-state annotation. The design includes 2 COPD-IV and 2 non-CLD donors across seven ALI differentiation time points.", format(nrow(meta_all),big.mark=","), format(nrow(meta),big.mark=",")),
  "Expression validation uses the author-provided combined day-0 and day-28 H5AD objects. All inference is based on donor summaries; cells are not treated as independent replicates.", "",
  "## State localization", "",
  sprintf("- Day-28 Basal_combined: ABNORMAL_REPAIR %.3f; KRT14_KRT17_REPAIR %.3f; CILIA_CONSENSUS %.3f; MATURE_CILIATED %.3f; imbalance %.3f.",
    profiles[time_point=="day 28" & cell_type=="Basal_combined"]$ABNORMAL_REPAIR,
    profiles[time_point=="day 28" & cell_type=="Basal_combined"]$KRT14_KRT17_REPAIR,
    profiles[time_point=="day 28" & cell_type=="Basal_combined"]$CILIA_CONSENSUS,
    profiles[time_point=="day 28" & cell_type=="Basal_combined"]$MATURE_CILIATED,
    profiles[time_point=="day 28" & cell_type=="Basal_combined"]$Repair_Cilia_Imbalance),
  sprintf("- Day-28 Ciliated: ABNORMAL_REPAIR %.3f; KRT14_KRT17_REPAIR %.3f; CILIA_CONSENSUS %.3f; MATURE_CILIATED %.3f; imbalance %.3f.",
    profiles[time_point=="day 28" & cell_type=="Ciliated"]$ABNORMAL_REPAIR,
    profiles[time_point=="day 28" & cell_type=="Ciliated"]$KRT14_KRT17_REPAIR,
    profiles[time_point=="day 28" & cell_type=="Ciliated"]$CILIA_CONSENSUS,
    profiles[time_point=="day 28" & cell_type=="Ciliated"]$MATURE_CILIATED,
    profiles[time_point=="day 28" & cell_type=="Ciliated"]$Repair_Cilia_Imbalance), "",
  "## Day-28 composition", "")
for(i in seq_len(nrow(key_comp))){z<-key_comp[i];lines<-c(lines,sprintf("- %s: %.2f percentage points (95%% CI %.2f to %.2f), P=%s, FDR=%s.",z$cell_type,z$estimate_percentage_points,z$CI95_lower_percentage_points,z$CI95_upper_percentage_points,format(z$p_value,digits=3),format(z$FDR,digits=3)))}
lines <- c(lines,"","## Within-state/signature inference","")
for(i in seq_len(nrow(key_sig))){z<-key_sig[i];lines<-c(lines,sprintf("- %s | %s | %s: effect %.3f (95%% CI %.3f to %.3f), P=%s, FDR=%s.",z$time_point,z$cell_type,z$signature,z$estimate,z$CI95_lower,z$CI95_upper,format(z$p_value,digits=3),format(z$FDR,digits=3)))}
lines <- c(lines,"","## Assessment", "",
  "The independent dataset reproduces the repair-high/cilia-low basal-versus-mature-ciliated state axis. Disease-group contrasts are directionally supportive for fewer mature ciliated cells, accumulation of transitional-ciliated cells, and higher basal repair/KRT14-KRT17/imbalance, but none passes BH-FDR and some ciliated/all-epithelial intrinsic effects are discordant.",
  "Overall classification: directionally supportive.",
  "Because there are only two donors per group, effect directions and biological localization are more informative than nominal significance; negative or discordant effects are retained.")
writeLines(lines,file.path(project,"02_scRNA/SECOND_VALIDATION/second_stage_scRNA_summary.md"))
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
