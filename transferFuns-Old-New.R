integrateByclustColOldNew <- function(rna_new,
                                              rna0,
                                              clustCol ="numclust",
                                              refAssayLabel = "RNA",
                                              refAssayTransfer = "RNA",
                                              condCol,
                                              #condition,
                                              dimsLSI = c(1:30),
                                              dimsCCA = c(1:30),
                                              outdir = "./") {
  
  
  Idents(rna0) <- condCol
  rna_old <- subset(rna0,idents = condition)
  
  ##Data preprocessing
  
  ## Here, we process the gene activity matrix in order to find anchors between cells in the scATAC-seq dataset and the scRNA-seq dataset.
  rna_old$numclust <- rna_old@meta.data[,clustCol]
  DefaultAssay(atac) <- "activities"
  atac <- FindVariableFeatures(atac)
  atac <- NormalizeData(atac)
  atac <- ScaleData(atac)
  
  transfer.anchors <- FindTransferAnchors(reference = rna_old, query = rna_new, features = VariableFeatures(object = rna_old))
  
  ## To improve discard first integratedLSI ?
  celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna_old$numclust)
  rna_new <- AddMetaData(atac, metadata = celltype.predictions)
  
  png(paste0(outdir,"prediction.score.max.png"))
  hist(rna_new$prediction.score.max)
  abline(v = 0.5, col = "red")
  dev.off()
  
  ## We filter atac cells with a weak prediction score for rna cluster transfer
  rna_new.filtered <- subset(rna_new, subset = prediction.score.max > 0.4)
  rna_new.filtered$predicted.id <- factor(rna_new.filtered$predicted.id, levels = levels(rna_old$numclust))  # to make the colors match
  
  p1 <- DimPlot(rna_new.filtered, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("New scRNA-seq cells") + 
    NoLegend() + scale_colour_hue(drop = FALSE)
  p2 <- DimPlot(rna_old, group.by = "numclust", label = TRUE, repel = TRUE) + ggtitle("Old scRNA-seq cells") + 
    NoLegend()
  
  png(paste0(opt$outdir,"predictedNewvsOldclust.png"),units = "in",res = 300,height=6,width =10)
  grid.arrange(p1,p2,nrow = 1)
  dev.off()
  
  ### Input with all rna data as ref
  # note that we restrict the imputation to variable genes from scRNA-seq 
  genes.use <- VariableFeatures(rna0) # all genes of the slot if integrated 
  refdata <- GetAssayData(rna0, assay = refAssayTransfer, slot = "data")[genes.use, ]
  
  # refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
  # (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
  transfer.anchors.2 <- FindTransferAnchors(reference = rna0, 
                                            query = rna_new.filtered, 
                                            features = VariableFeatures(object = rna0), 
                                            reference.assay = refAssayTransfer, 
                                            query.assay = "activities", 
                                            reduction = "cca",
                                            dims = dimsCCA)
  
  imputation <- TransferData(anchorset = transfer.anchors.2, 
                             refdata = refdata, 
                             weight.reduction = rna_new.filtered[["lsi"]],
                             dims = dimsLSI)
  
  # this line adds the imputed data matrix to the atac object
  atac.filtered[[refAssayTransfer]] <- imputation
  
  return(atac.filtered)
}
