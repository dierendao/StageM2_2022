
CorPCPlot <- function(object, assay = NULL,reduction = "lsi", n = 10, metaCol, ...) 
{
  #assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  dr <- object[[reduction]]
  embed <- Embeddings(object = dr)
  counts <- object[[metaCol]]
  embed <- embed[rownames(x = counts), ]
  #n <- SetIfNull(x = n, y = ncol(x = embed))
  embed <- embed[, seq_len(length.out = n)]
  depth.cor <- as.data.frame(cor(x = embed, y = counts))
  depth.cor$counts <- depth.cor[, 1]
  depth.cor$Component <- seq_len(length.out = nrow(x = depth.cor))
  p <- ggplot(depth.cor, aes(Component, counts)) + geom_point() + 
    scale_x_continuous(n.breaks = n, limits = c(1, n)) + 
    ylab("Correlation") + ylim(c(-1, 1)) + theme_light() + 
    ggtitle("Correlation between depth and reduced dimension components", 
            subtitle = paste0("Assay: ", assay, "\t", "Reduction: ", 
                              reduction))
  return(p)
}


rna <- readRDS("output/RNA/PLZF_RARA_CT/Seurat4/seurat.rds")
atac <- readRDS("output/ATAC/smp/Ctrl/atac.rds")
DimPlot(atac)


##Data preprocessing

## Here, we process the gene activity matrix in order to find anchors between cells in the scATAC-seq dataset and the scRNA-seq dataset.

DefaultAssay(atac) <- "activities"
atac <- FindVariableFeatures(atac)
atac <- NormalizeData(atac)
atac <- ScaleData(atac)

transfer.anchors <- FindTransferAnchors(reference = rna, query = atac, features = VariableFeatures(object = rna), 
                                        reference.assay = "RNA", query.assay = "activities", reduction = "cca")

## To improve discard first integratedLSI ?
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna$numclust, dims = 1:30,
                                     weight.reduction = atac[["lsi"]])
atac <- AddMetaData(atac, metadata = celltype.predictions)

hist(atac$prediction.score.max)
abline(v = 0.5, col = "red")

atac.filtered <- subset(atac, subset = prediction.score.max > 0.5)
atac.filtered$predicted.id <- factor(atac.filtered$predicted.id, levels = levels(rna))  # to make the colors match
p1 <- DimPlot(atac.filtered, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + 
  NoLegend() + scale_colour_hue(drop = FALSE)
p2 <- DimPlot(rna, group.by = "numclust", label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells") + 
  NoLegend()
p1 + p2


# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(rna)
refdata <- GetAssayData(rna, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
## To improve discard first integratedLSI ?
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = atac[["lsi"]])

# this line adds the imputed data matrix to the atac object
atac[["RNA"]] <- imputation
coembed <- merge(x = rna, y = atac)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)
coembed$celltype <- ifelse(!is.na(coembed$numclust), coembed$numclust, coembed$predicted.id)


p1 <- DimPlot(coembed, group.by = "sampleName")
p2 <- DimPlot(coembed, group.by = "celltype", label = TRUE, repel = TRUE)
p1 + p2


