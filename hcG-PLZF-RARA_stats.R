suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(grid))
suppressMessages(library(S4Vectors))
suppressMessages(library(patchwork))
suppressMessages(library(ggpubr))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(cowplot))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(scales))
suppressMessages(library(gProfileR))
suppressMessages(library(biomaRt))


## Directoies

inputSeurat <- "/shared/projects/scrna_gmp_rara/stage_babacar/scranseq/analysis/output/RNA/Old_New_integration/rna_new.rds"

outdir <- "/shared/projects/scrna_gmp_rara/stage_babacar/scranseq/analysis/output/RNA/Old_New_integration/Analysis/hcG-PLZF-RARA_stats"

# Assay choice
assay <- "RNA"
seurat <- readRDS(inputSeurat)
rna_new = seurat

outdir <- outdir

# colors as factor
colorTreatment <- c("#664CFF", "#FF8000")
colorCluster <- c("#E69F00", "#CC79A7", "#0072B2", "#009E73", "#D55E00", "#56B4E9")
colorPhases <- c(brewer.pal(9,"RdPu"))[c(3,6,9)]
colorBranches <- c("#332288", "#117733", "#DDCC77", "#AA4499", "#88CCEE")

DefaultAssay(seurat) <- assay


# compare mean expression of hcG-PLZF-RARA per cluster per treatment

p1 = VlnPlot(rna_new, feature = "PLZF-RARA-hcG", pt.size = 0, group.by = "predicted.id", split.by = "sampleName", 
             cols = rep(colorTreatment)) 
# Add pairwise comparisons p-value
p1 = p1 + stat_compare_means(method = "anova", label = "p.signif") 

png(paste(outdir,"/hcG-PLZF-RARA_MeanComparisonPerClusterPerTreatment.png", sep = "") ,res = 300,height=6,width =10,units = "in") # 
grid.arrange(p1)
dev.off()

# compare mean expression of hcG-PLZF-RARA per treatment

#my_comparisons <- list(c("RNA_Ctrl", "RNA_RA"))
p2 = VlnPlot(rna_new, feature = "PLZF-RARA-hcG", pt.size = 0, group.by = "sampleName", cols = rep(colorTreatment)) + NoLegend()
p2 = p2 + stat_compare_means(method = "t.test", label = "p.signif")

png(paste(outdir,"/hcG-PLZF-RARA_MeanComparisonPerTreatment.png", sep = "") ,res = 300,height=6,width =10,units = "in") # 
grid.arrange(p2)
dev.off()


# compare mean expression of hcG-PLZF-RARA per cluster

#my_comparisons <- list(c("Rep", "Neu1", "Neu2", "Neu3", "NeuRA1", "NeuRA2"))
p3 = VlnPlot(rna_new, feature = "PLZF-RARA-hcG", pt.size = 0, group.by = "predicted.id", cols = rep(colorCluster)) 
p3 = p3 + stat_compare_means(method = "anova", label = "p.signif")

png(paste(outdir,"/hcG-PLZF-RARA_MeanComparisonPerCluster.png", sep = "") ,res = 300,height=6,width =10,units = "in") # 
grid.arrange(p3)
dev.off()

# save plots
png(paste(outdir,"/hcG-PLZF-RARA_MeanComparison.png", sep = "") ,res = 300,height=6,width =10,units = "in") # 
grid.arrange(p1, p2, p3, nrow = 1)
dev.off()


## Per cluster
Idents(rna_new) = "predicted.id"
for (cl1 in levels(rna_new@meta.data$predicted.id)) {
  for (cl2 in levels(rna_new@meta.data$predicted.id)) {
    DE = FindMarkers(rna_new, features = "PLZF-RARA-hcG", ident.1 = cl1 , ident.2 = cl2, logfc.threshold = 0.0, min.pct = 0.0, verbose = FALSE)
    if (dim(DE)[1] == 0) {
      break }
    #else if (cl1 != "NeuRA2" && cl2 != "NeuRA2"&& cl1 !=  cl2 && abs(DE$p_val < 0.05)) 
    else {#if (abs(DE$p_val < 0.05)) {
      cat(cl1,"  ", cl2); cat("\n")
      print(DE); cat("\n\n") }
  }
}

## Per cluster per treatment
list_data = list()
Idents(rna_new) = "Treatment_clust"
for (cluster in levels(rna_new@meta.data$predicted.id)) {
  list_data[[cluster]] = FindMarkers(rna_new, features = "PLZF-RARA-hcG", ident.1 = paste(cluster, "RNA_Ctrl", sep = "_"), ident.2 = paste(cluster, "RNA_RA", sep = "_"), logfc.threshold = 0.0, min.pct = 0.1, verbose = FALSE)
  list_data[[cluster]]$cluster = cluster
  #print(list_data[[cluster]]$cluster)
}

list_to_df = do.call("rbind", list_data)

write.table(x = list_to_df, file = paste0(outdir, "MarkersCtrlVsTreated_table.tsv"), quote = FALSE)


## Per treatment
rna_new@meta.data$sampleName = factor(rna_new@meta.data$sampleName, levels = c("RNA_Ctrl", "RNA_RA"))
FindMarkers(rna_new, features = "PLZF-RARA-hcG", group.by = "sampleName", ident.1 = "RNA_Ctrl", ident.2 = "RNA_RA", logfc.threshold = 0.0, min.pct = 0.0, verbose = FALSE)


