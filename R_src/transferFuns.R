
integrateByConditionOneAtacRnaAll <- function(atac,
                                              rna0,
                                              clustCol ="numclust",
                                              refAssayLabel = "RNA",
                                              refAssayTransfer = "RNA",
                                              condCol,
                                              condition,
                                              dimsLSI = c(1:30),
                                              dimsCCA = c(1:30),
                                              outdir = "./") {

  
  Idents(rna0) <- condCol
  rna <- subset(rna0,idents = condition)
  
  ##Data preprocessing
  
  ## Here, we process the gene activity matrix in order to find anchors between cells in the scATAC-seq dataset and the scRNA-seq dataset.
  rna$numclust <- rna@meta.data[,clustCol]
  DefaultAssay(atac) <- "activities"
  atac <- FindVariableFeatures(atac)
  atac <- NormalizeData(atac)
  atac <- ScaleData(atac)
  
  transfer.anchors <- FindTransferAnchors(reference = rna, query = atac, features = VariableFeatures(object = rna), 
                                          reference.assay = refAssayLabel, query.assay = "activities", reduction = "cca")
  
  ## To improve discard first integratedLSI ?
  celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna$numclust, dims = dimsLSI,
                                       weight.reduction = atac[["lsi"]])
  atac <- AddMetaData(atac, metadata = celltype.predictions)
  
  png(paste0(outdir,"prediction.score.max.png"))
  hist(atac$prediction.score.max)
  abline(v = 0.5, col = "red")
  dev.off()
  
  ## We filter atac cells with a weak prediction score for rna cluster transfer
  atac.filtered <- subset(atac, subset = prediction.score.max > 0.4)
  atac.filtered$predicted.id <- factor(atac.filtered$predicted.id, levels = levels(rna$numclust))  # to make the colors match
  
  p1 <- DimPlot(atac.filtered, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + 
    NoLegend() + scale_colour_hue(drop = FALSE)
  p2 <- DimPlot(rna, group.by = "numclust", label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells") + 
    NoLegend()
  
  png(paste0(outdir,"predictedATACvsRNAclust.png"),units = "in",res = 300,height=6,width =10)
  grid.arrange(p1,p2,nrow = 1)
  dev.off()
  
  ### Input with all rna data as ref
  # note that we restrict the imputation to variable genes from scRNA-seq 
  genes.use <- VariableFeatures(rna0) # all genes of the slot if integrated 
  refdata <- GetAssayData(rna0, assay = refAssayTransfer, slot = "data")[genes.use, ]
  
  # refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
  # (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
  transfer.anchors.2 <- FindTransferAnchors(reference = rna0, 
                                            query = atac.filtered, 
                                            features = VariableFeatures(object = rna0), 
                                            reference.assay = refAssayTransfer, 
                                            query.assay = "activities", 
                                            reduction = "cca",
                                            dims = dimsCCA)
  
  imputation <- TransferData(anchorset = transfer.anchors.2, 
                             refdata = refdata, 
                             weight.reduction = atac.filtered[["lsi"]],
                             dims = dimsLSI)
  
  # this line adds the imputed data matrix to the atac object
  atac.filtered[[refAssayTransfer]] <- imputation
  
  return(atac.filtered)
}
