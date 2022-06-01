## -----------------------------------------------------------------------------
## Libraries
## -----------------------------------------------------------------------------

suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(gProfileR)
library(biomaRt)
library(getopt)

## -----------------------------------------------------------------------------
## Load function
## -----------------------------------------------------------------------------

source("/shared/projects/scrna_gmp_rara/stage_babacar/scranseq/analysis/R_src/funForSeuratAnalysis.R")

## -----------------------------------------------------------------------------
## Command line args
## -----------------------------------------------------------------------------

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  "inputSeurat",  "i",1, "character", "input seurat object with cluster column names numclust",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  "resolution", "r", 1, "character", "clustering resolution to use",
  "SigList", "s", 1, "character", "Path to the directory with all the signatures (.txt)",
  "clustName", "c", 1, "character", "List of cluster name sep by +",
  "logfc_threshold", "f", 1, "numeric", "Set a logFC threshold",
  "pval", "p", 1, "numeric", "Set a pval adj threshold",
  "assay", "a", 1, "character", "Chose default assay"
), byrow=TRUE, ncol=5);

opt = getopt(spec)


inputNew = "/shared/projects/scrna_gmp_rara/stage_babacar/scranseq/analysis/output/RNA/Seurat4_integration/combined.rds"
outdir = "/shared/projects/scrna_gmp_rara/stage_babacar/scranseq/analysis/output/RNA/IntegrationAnalysis/"
resolution <- "integrated_snn_res.0.3"
assay <- "RNA"
clustName <- "Cl1+Cl2+Cl3+Cl4+Cl5+Cl6+Cl7"
logfc_threshold <- 0.25
pval <- 0.05

integrated = readRDS(inputNew)
DefaultAssay(integrated) = assay
Idents(integrated) = resolution

# Prepare cluster name
clustName <- strsplit(clustName, split = "\\+")[[1]]

theme_set(theme_classic())
dir.create(outdir, recursive = T)

## ---------------------------------------------------------------------
## Color scheme
## ---------------------------------------------------------------------

colorTreatment <- c("#664CFF", "#FF8000")
colorCluster <- c("#E69F00", "#CC79A7", "#0072B2", "#009E73", "#D55E00", "#56B4E9", "#FF0000")
colorPhases <- c(brewer.pal(9,"RdPu"))[c(3,6,9)]
colorBranches <- c("#332288", "#117733", "#DDCC77", "#AA4499", "#88CCEE")


##-------------------------------------------------#
##       Data preparation : Clust name etc
##-------------------------------------------------#

names(clustName) = levels(integrated)
integrated = RenameIdents(integrated, clustName)
integrated@meta.data$FinalCluster = Idents(integrated)
integrated@meta.data$FinalCluster = factor(integrated@meta.data$FinalCluster, levels = c("Cl1", "Cl2", "Cl3", "Cl4", "Cl5", "Cl6", "Cl7"))
Idents(integrated) = "FinalCluster"

# add Treatment_clust as metadata
integrated$Treatment_clust = paste(integrated$FinalCluster, integrated$sampleName, sep = "_")
integrated$Treatment_clust = factor(integrated$Treatment_clust, levels = c("Cl1_RNA_Ctrl", "Cl1_RNA_RA", "Cl2_RNA_Ctrl", "Cl2_RNA_RA", "Cl3_RNA_Ctrl", "Cl3_RNA_RA", "Cl4_RNA_Ctrl", "Cl4_RNA_RA", "Cl5_RNA_Ctrl", "Cl5_RNA_RA", "Cl6_RNA_Ctrl", "Cl6_RNA_RA", "Cl7_RNA_Ctrl", "Cl7_RNA_RA")) 


# plot labels 
p <- DimPlot(integrated, group.by = "FinalCluster", label = TRUE, repel = TRUE, cols = colorCluster) + ggtitle("Labels") + 
  NoLegend() #+ scale_colour_hue(drop = FALSE)

png(paste(outdir,"/ClustName.png", sep = ""),res = 300,height=6,width =10,units = "in") # 
grid.arrange(p)
dev.off()


##--------------------------------------------------#
##   Looking for cluster markers with new names
##--------------------------------------------------#

Idents(integrated) <- "FinalCluster"
markers_cluster <- FindAllMarkers(integrated, logfc.threshold = logfc_threshold)
markers_cluster_sig <- markers_cluster[markers_cluster$p_val_adj < pval,]

write.table(x = markers_cluster_sig, file = paste0(outdir, "Markers_table.tsv"), quote = F, col.names = T, row.names = F, sep = "\t")


##-----------------------------------------------------#
##   Looking for geneSet enrichment on cluster markers
##-----------------------------------------------------#

enrich_clust <- list()
for(cluster in unique(markers_cluster$cluster)){
  print(cluster)
  #enrich_clust[[cluster]] <- gprofileR2::gost(markers_cluster[markers_cluster$cluster %in% cluster, "gene"],
  #organism = "mmusculus",
  #custom_bg = rownames(seurat),
  #ordered_query = "none")
  
  markers_cluster_gost = markers_cluster[markers_cluster$cluster %in% cluster, "gene"]
  enrich_clust[[cluster]] <- gProfileR::gprofiler(query = markers_cluster_gost,
                                                  organism = "mmusculus",
                                                  custom_bg = rownames(seurat),
                                                  ordered_query = "none")
  enrich_clust[[cluster]]$cluster <- cluster
}
enrich_clust.df <- do.call("rbind", enrich_clust)

write.table(enrich_clust.df, file = paste0(outdir, "enrichmentTable.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)



top_markers_cl1 = head(markers_cluster_sig[markers_cluster$cluster == "Cl1", ]$gene, 10)
top_markers_cl2 = head(markers_cluster_sig[markers_cluster$cluster == "Cl2", ]$gene, 10)
top_markers_cl3 = head(markers_cluster_sig[markers_cluster$cluster == "Cl3", ]$gene, 10)
top_markers_cl4 = head(markers_cluster_sig[markers_cluster$cluster == "Cl4", ]$gene, 10)
top_markers_cl5 = head(markers_cluster_sig[markers_cluster$cluster == "Cl5", ]$gene, 10)
top_markers_cl6 = head(markers_cluster_sig[markers_cluster$cluster == "Cl6", ]$gene, 10)
top_markers_cl7 = head(markers_cluster_sig[markers_cluster$cluster == "Cl7", ]$gene, 10)

write.table(markers_cluster_sig[markers_cluster$cluster == "Cl1", ]$gene, file = paste0(outdir, "markers_cl1.csv"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(markers_cluster_sig[markers_cluster$cluster == "Cl2", ]$gene, file = paste0(outdir, "markers_cl2.csv"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(markers_cluster_sig[markers_cluster$cluster == "Cl3", ]$gene, file = paste0(outdir, "markers_cl3.csv"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(markers_cluster_sig[markers_cluster$cluster == "Cl4", ]$gene, file = paste0(outdir, "markers_cl4.csv"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(markers_cluster_sig[markers_cluster$cluster == "Cl5", ]$gene, file = paste0(outdir, "markers_cl5.csv"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(markers_cluster_sig[markers_cluster$cluster == "Cl6", ]$gene, file = paste0(outdir, "markers_cl6.csv"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(markers_cluster_sig[markers_cluster$cluster == "Cl7", ]$gene, file = paste0(outdir, "markers_cl7.csv"), row.names = FALSE, col.names = FALSE, quote = FALSE)

##-----------------------------------------------------#
##   Looking for geneSet enrichment on cluster markers
##-----------------------------------------------------#

enrich_clust <- list()
for(cluster in unique(markers_cluster$cluster)){
  print(cluster)
  enrich_clust[[cluster]] <- gProfileR::gprofiler(markers_cluster[markers_cluster$cluster %in% cluster, "gene"],
                                                  organism = "mmusculus",
                                                  custom_bg = rownames(integrated),
                                                  ordered_query = "none")
  enrich_clust[[cluster]]$cluster <- cluster
}
enrich_clust.df <- do.call("rbind", enrich_clust)

write.table(enrich_clust.df, file = paste0(outdir, "enrichmentTable.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)


#-------------------------------------------------#
#              Basic UMAP plotting
#-------------------------------------------------#
dir.create(paste0(outdir, "UMAP/"), recursive = T)

UMAP_Clust <- DimPlot(integrated, pt.size = 0.4, cols = colorCluster, group.by = "FinalCluster") + theme(axis.title = element_blank(), axis.text = element_blank())
UMAP_treatment <- DimPlot(integrated, split.by = "sampleName", group.by = "sampleName", pt.size = 0.4, cols = colorTreatment) + theme(axis.title = element_blank(), axis.text = element_blank())

# UMAP treatment split
Idents(integrated) <- "sampleName"
seurat_control <- subset(integrated, idents = "RNA_Ctrl")
seurat_treated <- subset(integrated, idents = "RNA_RA")
Idents(integrated) <- "FinalCluster"

UMAP_control <- DimPlot(seurat_control, split.by = "sampleName", group.by = "sampleName", pt.size = 0.4, cols = colorTreatment[1]) + theme(axis.title = element_blank(), axis.text = element_blank(), legend.position = "none")
UMAP_treated <- DimPlot(seurat_treated, split.by = "sampleName", group.by = "sampleName", pt.size = 0.4, cols = colorTreatment[2]) + theme(axis.title = element_blank(), axis.text = element_blank(), legend.position = "none")

UMAP_CellCycle <- DimPlot(integrated, pt.size = 0.4, cols = colorPhases, group.by = "phases") + theme(axis.title = element_blank(), axis.text = element_blank())
UMAP_CellCycle_split <- DimPlot(integrated, split.by = "sampleName", group.by = "phases", pt.size = 0.4, cols = colorPhases) + theme(axis.title = element_blank(), axis.text = element_blank())

ggsave(paste0(outdir,"UMAP/", 'UMAP_CellCycle', '.png'), plot = UMAP_CellCycle, device = 'png', path = NULL,
       scale = 1, width = 15, height = 15, units = 'cm', dpi = 300)

ggsave(paste0(outdir,"UMAP/", 'UMAP_split_CellCycle', '.png'), plot = UMAP_CellCycle_split, device = 'png', path = NULL,
       scale = 1, width = 20, height = 15, units = 'cm', dpi = 300)
ggsave(paste0(outdir,"UMAP/", 'UMAP_cluster', '.png'), plot = UMAP_Clust, device = 'png', path = NULL,
       scale = 1, width = 15, height = 15, units = 'cm', dpi = 300)

ggsave(paste0(outdir, "UMAP/", 'UMAP_split_treatment', '.png'), plot = UMAP_treatment, device = 'png', path = NULL,
       scale = 1, width = 20, height = 15, units = 'cm', dpi = 300)
ggsave(paste0(outdir, "UMAP/", 'UMAP_treatment_Control', '.png'), plot = UMAP_control, device = 'png', path = NULL,
       scale = 1, width = 15, height = 15, units = 'cm', dpi = 300)
ggsave(paste0(outdir, "UMAP/", 'UMAP_treatment_Treated', '.png'), plot = UMAP_treated, device = 'png', path = NULL,
       scale = 1, width = 15, height = 15, units = 'cm', dpi = 300)


#-------------------------------------------------#
# Generate the barplot for cell distrib in cluster
#-------------------------------------------------#
#BP of cell cycle distrib per cluster
summary_phases <- ddply(integrated@meta.data,~FinalCluster + phases + sampleName, nrow)

bp_phasesClust_control <- ggplot(data.frame(summary_phases[summary_phases$sampleName %in% "RNA_Ctrl",]), aes(fill = phases,y = V1, x=FinalCluster)) +
  geom_bar( stat="identity", position="fill")+
  scale_fill_manual( values= colorPhases)+
  scale_y_continuous(name = "Sample (%)", labels = c(0,25,50,75,100))+
  ylab(label = "")+xlab(label = "") + coord_flip()+
  theme(legend.title=element_blank()) + ggtitle("Control")

bp_phasesClust_treated <- ggplot(data.frame(summary_phases[summary_phases$sampleName %in% "RNA_RA",]), aes(fill = phases,y = V1, x=FinalCluster)) +
  geom_bar( stat="identity", position="fill")+
  scale_fill_manual( values= colorPhases)+
  scale_y_continuous(name = "Sample (%)", labels = c(0,25,50,75,100))+
  ylab(label = "")+xlab(label = "") + coord_flip()+
  theme(legend.title=element_blank()) + ggtitle("Treated")

cellCycle_legend <- get_legend(bp_phasesClust_treated)

bp_grid_cellCycle <- plot_grid(bp_phasesClust_control + theme(legend.position = "none"), bp_phasesClust_treated + theme(legend.position = "none"))
bp_cellCycle_leg_grid <- plot_grid(bp_grid_cellCycle, cellCycle_legend,  rel_widths = c(3,0.4))

ggsave(paste0(outdir, "BarPlot/", 'BP_cellCycle_split', '.png'), plot = bp_cellCycle_leg_grid, device = 'png', path = NULL,
       scale = 1, width = 20, height = 20, units = 'cm', dpi = 300)

CellPerClust <- ddply(integrated@meta.data,~FinalCluster + sampleName, nrow)

bp_rawCellCount_control <- ggplot(data.frame(CellPerClust[CellPerClust$sampleName %in% "RNA_Ctrl",]), aes(y = V1, x=FinalCluster)) +
  geom_bar( stat="identity") + coord_flip() + ylab("nCells") + expand_limits(y = c(NULL, 3000))
bp_rawCellCount_treated <- ggplot(data.frame(CellPerClust[CellPerClust$sampleName %in% "RNA_RA",]), aes(y = V1, x=FinalCluster)) +
  geom_bar( stat="identity") + coord_flip() + ylab("nCells") + expand_limits(y = c(NULL, 3000))

ggsave(paste0(outdir, "BarPlot/", 'BP_raw_cellCount_control', '.png'), plot = bp_rawCellCount_control, device = 'png', path = NULL,
       scale = 1, width = 15, height = 20, units = 'cm', dpi = 300)
ggsave(paste0(outdir, "BarPlot/", 'BP_raw_cellCount_treated', '.png'), plot = bp_rawCellCount_treated, device = 'png', path = NULL,
       scale = 1, width = 15, height = 20, units = 'cm', dpi = 300)

bp_control_raw_cycle <- plot_grid(bp_phasesClust_control + theme(legend.position = "none"), bp_rawCellCount_control + theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + ggtitle(" "), rel_widths = c(2,1))
bp_treated_raw_cycle <- plot_grid(bp_phasesClust_treated + theme(legend.position = "none"), bp_rawCellCount_treated + theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + ggtitle(" "), rel_widths = c(2,1))

bp_all_grid_noleg <- plot_grid(bp_control_raw_cycle, bp_treated_raw_cycle, nrow = 2)
bp_all_grid_leg <- plot_grid(bp_all_grid_noleg, cellCycle_legend, rel_widths = c(3,0.4))

ggsave(paste0(outdir, "BarPlot/", 'BP_grid_cellCycle_cellCount', '.png'), plot = bp_treated_raw_cycle, device = 'png', path = NULL,
       scale = 1, width = 20, height = 20, units = 'cm', dpi = 300)

# BP proportion treatment per cluster
clustersampleName <- ddply(integrated@meta.data,~FinalCluster + sampleName,nrow)
propExpect <- table(integrated@meta.data$sampleName)/length(integrated@meta.data$sampleName)[]

#enrich_clust <- getEnrichPopClust(hspc.combined = seurat, Xname = "Control", Yname = "Treated", colorX = #colorTreatment[1], colorY = colorTreatment[2], metaCol = "sampleName", clustCol = "predicted.id")

enrich_clust <- getEnrichPopClust(hspc.combined = integrated, Xname = "RNA_Ctrl", Yname = "RNA_RA", colorX = colorTreatment[1], colorY = colorTreatment[2], metaCol = "sampleName", clustCol = "FinalCluster")


clustersampleName$sampleName <- factor(clustersampleName$sampleName, levels = c("RNA_RA", "RNA_Ctrl"))
clustersampleName$condition <- factor(clustersampleName$sampleName, levels = c("RNA_Ctrl", "RNA_RA"))

bp_treatmentClust1 <- ggplot(data.frame(clustersampleName), aes(fill = sampleName,y = V1, x=FinalCluster)) +
  geom_bar( stat="identity", position="fill")+
  scale_fill_manual( values= rev(colorTreatment))+
  scale_y_continuous(name = "Sample (%)", labels = c(0,25,50,75,100))+
  ylab(label = "")+xlab(label = "") + theme_classic() + coord_flip()+
  theme(legend.title=element_blank(), axis.text.y = element_text(colour = enrich_clust[,"color"])) +
  geom_hline(yintercept = propExpect[1])

ggsave(paste0(outdir, "BarPlot/",'BP_treatment_fill', '.png'), plot = bp_treatmentClust1, device = 'png', path = NULL,
       scale = 1, width = 15, height = 20, units = 'cm', dpi = 300)


bp_treatmentClust2 <- ggplot(data.frame(clustersampleName), aes(fill = condition,y = V1, x=FinalCluster)) +
  #geom_bar( stat="identity", position="fill")+
  geom_bar( stat="identity", position = position_dodge2())+
  scale_fill_manual( values= colorTreatment)+
  ylab(label = "nCells")+xlab(label = "") +
  theme_classic() + coord_flip() +
  theme(legend.title=element_blank(), axis.text.y = element_text(colour = enrich_clust[,"color"])) 
#geom_hline(yintercept = propExpect[1])

ggsave(paste0(outdir, "BarPlot/",'BP_treatment_position_dodge', '.png'), plot = bp_treatmentClust2, device = 'png', path = NULL,
       scale = 1, width = 15, height = 20, units = 'cm', dpi = 300)

# BP cell number distribution per cluster 
CellPerClust <- ddply(integrated@meta.data,~FinalCluster, nrow)

bp_rawTreatmentClust <- ggplot(data.frame(CellPerClust), aes(y = V1, x=FinalCluster)) +
  geom_bar( stat="identity") + theme_classic() + coord_flip() + ylab("nCells") + xlab("")

ggsave(paste0(outdir, "BarPlot/",'BP_raw_treatment', '.png'), plot = bp_rawTreatmentClust, device = 'png', path = NULL,
       scale = 1, width = 15, height = 20, units = 'cm', dpi = 300)

# Arranging BP grid
treatment_legend <- get_legend(bp_treatmentClust2)

bp_grid1 <- plot_grid(bp_treatmentClust1 + theme(legend.position = "none", axis.line.y = element_blank(), 
                                                 axis.ticks.y = element_blank()), 
                      bp_treatmentClust2  + theme(legend.position = "none", axis.text.y = element_blank(), 
                                                  axis.title.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank()),
                      treatment_legend, rel_widths = c(2,2,0.5), nrow = 1)

ggsave(paste0(outdir, "BarPlot/",'BP_grid1_treatment_dodge', '.png'), plot = bp_grid1, device = 'png', path = NULL,
       scale = 1, width = 15, height = 20, units = 'cm', dpi = 300)


bp_grid2 <- plot_grid(bp_treatmentClust1 + theme(legend.position = "none", axis.line.y = element_blank(), 
                                                 axis.ticks.y = element_blank()), 
                      bp_rawTreatmentClust + theme(axis.text.y = element_blank(), axis.title.y = element_blank(), 
                                                   axis.line.y = element_blank(), 
                                                   axis.ticks.y = element_blank()), 
                      treatment_legend, rel_widths = c(2,1,0.5), nrow = 1)

ggsave(paste0(outdir, "BarPlot/",'BP_grid2_treatment_dodge', '.png'), plot = bp_grid2, device = 'png', path = NULL,
       scale = 1, width = 15, height = 20, units = 'cm', dpi = 300)


############################################################################################################ 
############################### Nb cells per cluster ###############################


nb_cells = summary(integrated$FinalCluster)

NbCellsPerCluster = data.frame(nb_cells)
NbCellsPerCluster$Cluster = rownames(NbCellsPerCluster)
NbCellsPerCluster$Cluster =  factor(NbCellsPerCluster$Cluster, levels = c("cl2", "cl6", "cl1", "cl4", "cl7", "cl3", "cl5"))

rownames(NbCellsPerCluster) = NULL
colnames(NbCellsPerClusterOldvsNew) = c("Old", "New", "Cluster")


p= ggplot(data = NbCellsPerCluster , aes(x = Cluster, y = nb_cells)) +
  geom_bar(stat = "identity", fill="#664CFF") +
  geom_text(aes(label=nb_cells), vjust=1.6, color="white", size=3.5) +
  theme_classic() +
  ggtitle("Nb cells per cluster with new fusion gene on reference") +
  ylab("Nb Cells") + xlab("")

############################################################################################################ 
############################### PLZF-RARA-hcG expression ###############################
############################################################################################################

#PLZF-RARA-hcG expression per cluster
p1 = VlnPlot(integrated, feature = "PLZF-RARA-hcG", pt.size = 0, group.by = "FinalCluster", cols = colorCluster) #+ ggtitle("PLZF-RARA-hcG expression")
p2 = DotPlot(integrated, feature = "PLZF-RARA-hcG", group.by = "FinalCluster") #+ ggtitle("PLZF-RARA-hcG expression")

png(paste(outdir,"/PLZF-RARA-hcG_expression.png", sep = ""),res = 300,height=6,width =10,units = "in") # 
grid.arrange(p1,p2,nrow = 1)
dev.off()


# PLZF-RARA-hcG expression per cluster per treatment
p1 = VlnPlot(integrated, feature = "PLZF-RARA-hcG", pt.size = 0, group.by = "Treatment_clust") + NoLegend()
png(paste(outdir,"/PLZF-RARA-hcG_expressionPerClusterPerTreatment.png", sep = ""),res = 300,height=6,width =10,units = "in") # 
grid.arrange(p1)
dev.off()


p = VlnPlot(integrated, feature = "PLZF-RARA-hcG", pt.size = 0, group.by = "FinalCluster", split.by = "sampleName", cols = rep(colorTreatment))
png(paste(outdir,"/PLZF-RARA-hcG_expressionPerClusterSplitPerSampleName.png", sep = ""),res = 300,height=6,width =10,units = "in") # 
grid.arrange(p)
dev.off()

# PLZF-RARA-hcG expression per cluster on CT
rna_new_CT = integrated[, integrated$Treatment_clust %in% c("Cl1_RNA_Ctrl", "Cl2_RNA_Ctrl", "Cl3_RNA_Ctrl", "Cl4_RNA_Ctrl", "Cl6_RNA_Ctrl", "Cl7_RNA_Ctrl")]
p = VlnPlot(rna_new_CT, feature = "PLZF-RARA-hcG", pt.size = 0, split.by = "sampleName", cols = rep(colorTreatment[1]))

png(paste(outdir,"/PLZF-RARA-hcG_expressionPerClusterCT.png", sep = ""),res = 300,height=6,width =10,units = "in") # 
grid.arrange(p)
dev.off()


## statitics associated with VlnPlots of PLZF-RARA-hcG

#PLZF-RARA-hcG expression per cluster  
Idents(integrated) = "FinalCluster"
PLZF_RARA_hcG_statsPerClust = FindAllMarkers(integrated, features = "PLZF-RARA-hcG", logfc.threshold = 0.0)
write.table(PLZF_RARA_hcG_statsPerClust, file = paste0(outdir, "PLZF_RARA_hcG_statsPerClust.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)


#PLZF-RARA-hcG expression per cluster per treatment 
Idents(rna_new) = "Treatment_clust"
PLZF_RARA_hcG_statsPerClustPerSample = FindAllMarkers(rna_new, features = "PLZF-RARA-hcG", logfc.threshold = 0.0)
write.table(PLZF_RARA_hcG_statsPerClustPerSample, file = paste0(outdir, "PLZF_RARA_hcG_statsPerClustPerSample.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

# PLZF-RARA-hcG expression per cluster on CT
Idents(rna_new_CT) = "Treatment_clust"
PLZF_RARA_hcG_CT_stats = FindAllMarkers(rna_new_CT, features = "PLZF-RARA-hcG", logfc.threshold = 0.0)
write.table(PLZF_RARA_hcG_CT_stats, file = paste0(outdir, "PLZF-RARA-hcG-CT_stats.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

############################################################################################################

Idents(integrated) = "FinalCluster"
for (cl1 in levels(integrated@meta.data$FinalCluster)) {
  for (cl2 in levels(integrated@meta.data$FinalCluster)) {
    DE = FindMarkers(integrated, features = "PLZF-RARA-hcG", ident.1 = cl1, ident.2 = cl2, logfc.threshold = 0.0, verbose = FALSE)
    if (dim(DE)[1] == 0) {
      break }
    else if (cl1 !=  cl2) {
      cat(cl1,"  ", cl2); cat("\n")
      print(DE); cat("\n\n") }
  }
}

Idents(integrated) = "Treatment_clust"

control = c("Cl1_RNA_Ctrl", "Cl2_RNA_Ctrl", "Cl3_RNA_Ctrl", "Cl4_RNA_Ctrl", "Cl5_RNA_Ctrl", "Cl6_RNA_Ctrl", "Cl7_RNA_Ctrl")
control = factor(control, levels = c("Cl1_RNA_Ctrl", "Cl2_RNA_Ctrl", "Cl3_RNA_Ctrl", "Cl4_RNA_Ctrl", "Cl5_RNA_Ctrl", "Cl6_RNA_Ctrl", "Cl7_RNA_Ctrl"))

treated = c("Cl1_RNA_RA", "Cl2_RNA_RA", "Cl3_RNA_RA", "Cl4_RNA_RA", "Cl5_RNA_RA", "Cl6_RNA_RA", "Cl7_RNA_RA")
treated = c("Cl1_RNA_RA", "Cl2_RNA_RA", "Cl3_RNA_RA", "Cl4_RNA_RA", "Cl5_RNA_RA", "Cl6_RNA_RA", "Cl7_RNA_RA")

for (cl1 in levels(integrated)) {
  for (cl2 in levels(treated)) {
    DE = FindMarkers(integrated, features = "PLZF-RARA-hcG", ident.1 = cl1, ident.2 = cl2, logfc.threshold = 0.0, verbose = FALSE)
    if (dim(DE)[1] == 0) {
      break }
    else if (cl1 !=  cl2 && DE$p_val < 0.05) {
      cat(cl1,"  ", cl2); cat("\n")
      print(DE); cat("\n\n") }
  }
}


list_data = list()
Idents(integrated) = "Treatment_clust"
for (cluster in levels(integrated@meta.data$FinalCluster)) {
  list_data[[cluster]] = FindMarkers(integrated, features = "PLZF-RARA-hcG", ident.1 = paste(cluster, "RNA_Ctrl", sep = "_"), ident.2 = paste(cluster, "RNA_RA", sep = "_"), logfc.threshold = 0.0, min.pct = 0.0, verbose = FALSE)
  list_data[[cluster]]$cluster = cluster
  print(list_data[[cluster]]$cluster)
}

list_to_df = do.call("rbind", list_data)

write.table(x = list_to_df, file = paste0(outdir, "MarkersCtrlVsTreated_table.tsv"), quote = FALSE)

p1 = VlnPlot(integrated, feature = "PLZF-RARA-hcG", pt.size = 0, group.by = "FinalCluster", split.by = "sampleName", 
             cols = rep(colorTreatment)) 
# Add pairwise comparisons p-value
p1 = p1 + stat_compare_means(method = "anova", label = "p.signif")

png(paste(outdir,"/hcG-PLZF-RARA_MeanComparisonPerClusterPerTreatment.png", sep = "") ,res = 300,height=6,width =10,units = "in") # 
grid.arrange(p1)
dev.off()

## save
Idents(seurat) <- "FinalCluster"
DefaultAssay(seurat) <- DefAssay
saveRDS(integrated, paste0(outdir, "seurat_annotated", ".rds"))

