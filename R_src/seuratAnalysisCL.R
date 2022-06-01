### Supplemental analysis and cluster annotation


## -----------------------------------------------------------------------------
## Libraries
## -----------------------------------------------------------------------------

suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
library(plyr)
library(dplyr)
library(stringr)
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

source("/shared/ifbstor1/projects/scrna_gmp_rara/stage_babacar/scranseq/analysis/R_src/funForSeuratAnalysis.R")

## -----------------------------------------------------------------------------
## Command line args
## -----------------------------------------------------------------------------

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  "inputSeurat",  "i",1, "character", "input seurat object with cluster column names numclust",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  "resolution", "r", 1, "character", "clustering resolution to use",
  #"SigList", "s", 1, "character", "Path to the directory with all the signatures (.txt)",
  #"clustName", "c", 1, "character", "List of cluster name sep by +",
  "logfc_threshold", "f", 1, "numeric", "Set a logFC threshold",
  "pval", "p", 1, "numeric", "Set a pval adj threshold",
  "assay", "a", 1, "character", "Chose default assay"
), byrow=TRUE, ncol=5);

opt = getopt(spec)

# if help was asked, print a friendly message
# and exit with a non-zero error code

args <- commandArgs()
if ( !is.null(opt$help) | is.null(opt$inputSeurat)) {
  cat("Create influence graph from regulon table")
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

## ---------------------------------------------------------------------
## Read inputs
## ---------------------------------------------------------------------

# # Testing opt
inputSeurat <- "/shared/projects/scrna_gmp_rara/stage_babacar/scranseq/analysis/output/RNA/Old_New_integration/rna_new.rds"
outdir <- "/shared/projects/scrna_gmp_rara/stage_babacar/scranseq/analysis/output/RNA/Old_New_integration/Analysis/"
assay <- "RNA"


# Read seurat object
seurat <- readRDS(inputSeurat)

outdir <- outdir

# Set logFC threshold
if (is.null(opt$logfc_threshold)) {
  opt$logfc_threshold <- 0.25
}

# Set logFC threshold
if (is.null(opt$pval)) {
  opt$pval <- 0.05
}

# Define ref clustering
if(is.null(opt$resolution)){
  opt$resolution <- "integrated_snn_res.0.3"
}
#targetRes <- opt$resolution

# Define default assay
if(is.null(opt$assay)){
  opt$assay <- "RNA"
}

DefAssay <- assay

# Prepare cluster name
#clustName <- strsplit(opt$clustName, split = "\\+")[[1]]


theme_set(theme_classic())
dir.create(outdir, recursive = T)

## ---------------------------------------------------------------------
## Color scheme
## ---------------------------------------------------------------------

colorTreatment <- c("#664CFF", "#FF8000")
colorCluster <- c("#E69F00", "#CC79A7", "#0072B2", "#009E73", "#D55E00", "#56B4E9")
colorPhases <- c(brewer.pal(9,"RdPu"))[c(3,6,9)]
colorBranches <- c("#332288", "#117733", "#DDCC77", "#AA4499", "#88CCEE")

##-------------------------------------------------#
##       Data preparation : Clust name etc
##-------------------------------------------------#

DefaultAssay(seurat) <- DefAssay

seurat@meta.data$phases <- factor(seurat@meta.data$phases, levels = c("G1_G0", "S", "G2_M"))

Idents(seurat) <- "predicted.id"
markers_cluster <- FindAllMarkers(seurat, logfc.threshold = opt$logfc_threshold)
markers_cluster_sig <- markers_cluster[markers_cluster$p_val_adj < opt$pval,]

# Rep markers
markers_Rep_sig = markers_cluster_sig[markers_cluster$cluster == "Rep", ]$gene

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


##-----------------------------------------------------#
##   Looking for gene enrichment on cluster markers
##-----------------------------------------------------#

top_markers_Rep = head(markers_cluster_sig[markers_cluster$cluster == "Rep", ]$gene, 10)
top_markers_Neu1 = head(markers_cluster_sig[markers_cluster$cluster == "Neu1", ]$gene, 10)
top_markers_Neu2 = head(markers_cluster_sig[markers_cluster$cluster == "Neu2", ]$gene, 10)
top_markers_Neu3 = head(markers_cluster_sig[markers_cluster$cluster == "Neu3", ]$gene, 10)
top_markers_NeuRA1 = head(markers_cluster_sig[markers_cluster$cluster == "NeuRA1", ]$gene, 10)
top_markers_NeuRA2 = head(markers_cluster_sig[markers_cluster$cluster == "NeuRA2", ]$gene, 10)

write.table(markers_cluster_sig[markers_cluster$cluster == "Rep", ]$gene, file = paste0(outdir, "markers_Rep.csv"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(markers_cluster_sig[markers_cluster$cluster == "Neu1", ]$gene, file = paste0(outdir, "markers_Neu1.csv"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(markers_cluster_sig[markers_cluster$cluster == "Neu2", ]$gene, file = paste0(outdir, "markers_Neu2.csv"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(markers_cluster_sig[markers_cluster$cluster == "Neu3", ]$gene, file = paste0(outdir, "markers_Neu3.csv"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(markers_cluster_sig[markers_cluster$cluster == "NeuRA1", ]$gene, file = paste0(outdir, "markers_NeuRA1.csv"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(markers_cluster_sig[markers_cluster$cluster == "NeuRA2", ]$gene, file = paste0(outdir, "markers_NeuRA2.csv"), row.names = FALSE, col.names = FALSE, quote = FALSE)


#-------------------------------------------------#
#              Basic UMAP plotting
#-------------------------------------------------#
dir.create(paste0(outdir, "UMAP/"), recursive = T)

UMAP_Clust <- DimPlot(seurat, pt.size = 0.4, cols = colorCluster, group.by = "predicted.id") + theme(axis.title = element_blank(), axis.text = element_blank())
UMAP_treatment <- DimPlot(seurat, split.by = "sampleName", group.by = "sampleName", pt.size = 0.4, cols = colorTreatment) + theme(axis.title = element_blank(), axis.text = element_blank())

# UMAP treatment split
Idents(seurat) <- "sampleName"
seurat_control <- subset(seurat, idents = "RNA_Ctrl")
seurat_treated <- subset(seurat, idents = "RNA_RA")
Idents(seurat) <- "predicted.id"

UMAP_control <- DimPlot(seurat_control, split.by = "sampleName", group.by = "sampleName", pt.size = 0.4, cols = colorTreatment[1]) + theme(axis.title = element_blank(), axis.text = element_blank(), legend.position = "none")
UMAP_treated <- DimPlot(seurat_treated, split.by = "sampleName", group.by = "sampleName", pt.size = 0.4, cols = colorTreatment[2]) + theme(axis.title = element_blank(), axis.text = element_blank(), legend.position = "none")

UMAP_CellCycle <- DimPlot(seurat, pt.size = 0.4, cols = colorPhases, group.by = "phases") + theme(axis.title = element_blank(), axis.text = element_blank())
UMAP_CellCycle_split <- DimPlot(seurat, split.by = "sampleName", group.by = "phases", pt.size = 0.4, cols = colorPhases) + theme(axis.title = element_blank(), axis.text = element_blank())

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
dir.create(paste0(outdir, "BarPlot/"), recursive = T)
#BP of cell cycle distrib per cluster
summary_phases <- ddply(seurat@meta.data,~predicted.id + phases + sampleName, nrow)

bp_phasesClust_control <- ggplot(data.frame(summary_phases[summary_phases$sampleName %in% "RNA_Ctrl",]), aes(fill = phases,y = V1, x=predicted.id)) +
  geom_bar( stat="identity", position="fill")+
  scale_fill_manual( values= colorPhases)+
  scale_y_continuous(name = "Sample (%)", labels = c(0,25,50,75,100))+
  ylab(label = "")+xlab(label = "") + coord_flip()+
  theme(legend.title=element_blank()) + ggtitle("Control")

bp_phasesClust_treated <- ggplot(data.frame(summary_phases[summary_phases$sampleName %in% "RNA_RA",]), aes(fill = phases,y = V1, x=predicted.id)) +
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

CellPerClust <- ddply(seurat@meta.data,~predicted.id + sampleName, nrow)

bp_rawCellCount_control <- ggplot(data.frame(CellPerClust[CellPerClust$sampleName %in% "RNA_Ctrl",]), aes(y = V1, x=predicted.id)) +
  geom_bar( stat="identity") + coord_flip() + ylab("nCells") + expand_limits(y = c(NULL, 3000))
bp_rawCellCount_treated <- ggplot(data.frame(CellPerClust[CellPerClust$sampleName %in% "RNA_RA",]), aes(y = V1, x=predicted.id)) +
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
clustersampleName <- ddply(seurat@meta.data,~predicted.id + sampleName,nrow)
propExpect <- table(seurat@meta.data$sampleName)/length(seurat@meta.data$sampleName)[]

enrich_clust <- getEnrichPopClust(hspc.combined = seurat, Xname = "RNA_Ctrl", Yname = "RNA_RA", colorX = colorTreatment[1], colorY = colorTreatment[2], metaCol = "sampleName", clustCol = "predicted.id")

clustersampleName$sampleName <- factor(clustersampleName$sampleName, levels = c("RNA_RA", "RNA_Ctrl"))

bp_treatmentClust <- ggplot(data.frame(clustersampleName), aes(fill = sampleName,y = V1, x=predicted.id)) +
  geom_bar( stat="identity", position="fill")+
  scale_fill_manual( values= rev(colorTreatment))+
  scale_y_continuous(name = "Sample (%)", labels = c(0,25,50,75,100))+
  ylab(label = "")+xlab(label = "") + coord_flip()+ theme_minimal() +
  theme(legend.title=element_blank(), axis.text.y = element_text(colour = enrich_clust[,"color"])) +
  geom_hline(yintercept = propExpect[1])

ggsave(paste0(outdir, "BarPlot/",'BP_treatment_fill', '.png'), plot = bp_treatmentClust, device = 'png', path = NULL,
       scale = 1, width = 15, height = 20, units = 'cm', dpi = 300)

# BP cell number distribution per cluster 
CellPerClust <- ddply(seurat@meta.data,~predicted.id, nrow)

bp_rawTreatmentClust <- ggplot(data.frame(CellPerClust), aes(y = V1, x=predicted.id)) +
  geom_bar( stat="identity") + coord_flip() + ylab("nCells")

ggsave(paste0(outdir, "BarPlot/",'BP_raw_treatment_fill', '.png'), plot = bp_rawTreatmentClust, device = 'png', path = NULL,
       scale = 1, width = 15, height = 20, units = 'cm', dpi = 300)

# Arranging BP grid
treatment_legend <- get_legend(bp_treatmentClust)

bp_grid <- plot_grid(bp_treatmentClust + theme(legend.position = "none"), 
                     bp_rawTreatmentClust + theme(axis.text.y = element_blank(), axis.title.y = element_blank()), 
                     treatment_legend, rel_widths = c(2,1,0.5), nrow = 1)

ggsave(paste0(outdir, "BarPlot/",'BP_grid_treatment_fill', '.png'), plot = bp_grid, device = 'png', path = NULL,
       scale = 1, width = 15, height = 20, units = 'cm', dpi = 300)



Idents(seurat) <- "predicted.id"
DefaultAssay(seurat) <- DefAssay

saveRDS(seurat, paste0(outdir, "seurat_annotated", ".rds"))


#################################################PLZF-RARA-hcG vs all analysis#################################################

# BP proportion PLZF-RARA-hcG per cluster per treatment

PLZF_RARA_hcG = subset(seurat, features = "PLZF-RARA-hcG")
PLZF_RARA_hcG_perSample <- ddply(PLZF_RARA_hcG@meta.data,~predicted.id + sampleName,nrow)

propExpect <- table(seurat@meta.data$sampleName)/length(seurat@meta.data$sampleName)[]

enrich_clust <- getEnrichPopClust(hspc.combined = seurat, Xname = "Control", Yname = "Treated", colorX = colorTreatment[1], colorY = colorTreatment[2], metaCol = "sampleName", clustCol = "predicted.id")

PLZF_RARA_hcG$condition <- factor(PLZF_RARA_hcG_perSample$sampleName, levels = c("RNA_Ctrl", "RNA_RA"))

bp_treatmentClust <- ggplot(data.frame(PLZF_RARA_hcG_perSample), aes(fill = condition,y = V1, x=predicted.id)) +
  geom_bar( stat="identity", position="fill")+
  scale_fill_manual( values= colorTreatment)+
  scale_y_continuous(name = "Sample (%)", labels = c(0,25,50,75,100))+
  ylab(label = "")+xlab(label = "") + coord_flip()+
  theme(legend.title=element_blank(), axis.text.y = element_text(colour = enrich_clust[,"color"])) +
  geom_hline(yintercept = propExpect[1])

ggsave(paste0(outdir, "BarPlot/",'BP_treatment', '.png'), plot = bp_treatmentClust, device = 'png', path = NULL,
       scale = 1, width = 15, height = 20, units = 'cm', dpi = 300)



# BP cell number distribution per cluster : PLZF_RARA_hcG
PLZF_RARA_hcG_CellPerClust <- ddply(PLZF_RARA_hcG@meta.data,~predicted.id, nrow)

bp_PLZF_RARA_hcG_CellPerClust <- ggplot(data.frame(PLZF_RARA_hcG_CellPerClust), aes(y = V1, x=predicted.id)) +
  geom_bar( stat="identity") + ylab("nCells") + xlab("Cluster") +  #+ coord_flip()
  geom_text(aes(label=V1), vjust=1.6, color="white", size=3.5) +
  theme_minimal() 

ggsave(paste0(outdir, "BarPlot/",'BP_PLZF_RARA_hcG_cluster', '.png'), plot = bp_PLZF_RARA_hcG_CellPerClust, device = 'png', path = NULL,
       scale = 1, width = 15, height = 20, units = 'cm', dpi = 300)


# BP cell number distribution per cluster per treatment : PLZF_RARA_hcG
PLZF_RARA_hcG = PLZF_RARA_hcG[, PLZF_RARA_hcG$nCount_RNA > 0]
CellPerClust <- ddply(PLZF_RARA_hcG@meta.data,~predicted.id, nrow)

PLZF_RARA_hcG_CellPerClustPerTreatment <- ddply(PLZF_RARA_hcG@meta.data,~Treatment_clust, nrow)
PLZF_RARA_hcG_CellPerClustPerTreatment$Treatment_clust = str_replace(PLZF_RARA_hcG_CellPerClustPerTreatment$Treatment_clust, "RNA_", "")
PLZF_RARA_hcG_CellPerClustPerTreatment$Treatment_clust = str_replace(PLZF_RARA_hcG_CellPerClustPerTreatment$Treatment_clust, "Ctrl", "CT")

PLZF_RARA_hcG_CellPerClustPerTreatment$Treatment_clust = factor(PLZF_RARA_hcG_CellPerClustPerTreatment$Treatment_clust, 
                                                                levels = c("Rep_CT", "Rep_RA", "Neu1_CT", "Neu1_RA", "Neu2_CT", "Neu2_RA", "Neu3_CT", "Neu3_RA", "NeuRA1_CT", "NeuRA1_RA", "NeuRA2_CT", "NeuRA2_RA"))


PLZF_RARA_hcG_RepCTProp = (PLZF_RARA_hcG_CellPerClustPerTreatment$V1[1] / CellPerClust$V1[1]) *100
PLZF_RARA_hcG_RepRAProp = (PLZF_RARA_hcG_CellPerClustPerTreatment$V1[2] / CellPerClust$V1[1]) *100

PLZF_RARA_hcG_Neu1CTProp = (PLZF_RARA_hcG_CellPerClustPerTreatment$V1[3] / CellPerClust$V1[2]) *100
PLZF_RARA_hcG_Neu1RAProp = (PLZF_RARA_hcG_CellPerClustPerTreatment$V1[4] / CellPerClust$V1[2]) *100

PLZF_RARA_hcG_Neu2CTProp = (PLZF_RARA_hcG_CellPerClustPerTreatment$V1[5] / CellPerClust$V1[3]) *100
PLZF_RARA_hcG_Neu2RAProp = (PLZF_RARA_hcG_CellPerClustPerTreatment$V1[6] / CellPerClust$V1[3]) *100

PLZF_RARA_hcG_Neu3CTProp = (PLZF_RARA_hcG_CellPerClustPerTreatment$V1[7] / CellPerClust$V1[4]) *100
PLZF_RARA_hcG_Neu3RAProp = (PLZF_RARA_hcG_CellPerClustPerTreatment$V1[8] / CellPerClust$V1[4]) *100

PLZF_RARA_hcG_NeuRA1CTProp = (PLZF_RARA_hcG_CellPerClustPerTreatment$V1[9] / CellPerClust$V1[5]) *100
PLZF_RARA_hcG_NeuRA1RAProp = (PLZF_RARA_hcG_CellPerClustPerTreatment$V1[10] / CellPerClust$V1[5]) *100

PLZF_RARA_hcG_NeuRA2CTProp = (PLZF_RARA_hcG_CellPerClustPerTreatment$V1[11] / CellPerClust$V1[6]) *100
PLZF_RARA_hcG_NeuRA2RAProp = (PLZF_RARA_hcG_CellPerClustPerTreatment$V1[12] / CellPerClust$V1[6]) *100


PLZF_RARA_hcG_CellPerClustPerTreatment$Prop = 
  c(PLZF_RARA_hcG_RepCTProp, PLZF_RARA_hcG_RepRAProp, PLZF_RARA_hcG_Neu1CTProp, PLZF_RARA_hcG_Neu1RAProp,
                    PLZF_RARA_hcG_Neu2CTProp, PLZF_RARA_hcG_Neu2RAProp, PLZF_RARA_hcG_Neu3CTProp, PLZF_RARA_hcG_Neu3RAProp, 
                    PLZF_RARA_hcG_NeuRA1CTProp, PLZF_RARA_hcG_NeuRA1RAProp, PLZF_RARA_hcG_NeuRA2CTProp, PLZF_RARA_hcG_NeuRA1RAProp)

bp_PLZF_RARA_hcG_CellPerClustPerTreatment <- ggplot(data.frame(PLZF_RARA_hcG_CellPerClustPerTreatment), aes(y = V1, x=Treatment_clust)) +
  geom_bar( stat="identity", fill = rep(colorTreatment, 6)) + ylab("nCells") + xlab("Cluster") + coord_flip() +
  geom_text(aes(label=V1), vjust = 1.3, color="black", size=2.5) + # vjust=1.6
  theme_minimal() 

ggsave(paste0(outdir, "BarPlot/",'NbCell-hcG-PLZF-RARA-PerClusterPerTreatement', '.png'), plot = bp_PLZF_RARA_hcG_CellPerClustPerTreatment, device = 'png', path = NULL,
       scale = 1, width = 15, height = 20, units = 'cm', dpi = 300)

