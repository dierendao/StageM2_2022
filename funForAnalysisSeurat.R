#Visualize position of a cluster on the UMAP
clusterDistrib_seurat <- function(seurat = seurat, numclust = numclust, metaData = "numclust", flip = FALSE){
  print(numclust)
  seurat@meta.data$clusterDis <- FALSE
  seurat@meta.data[which(seurat@meta.data[,metaData] == numclust),"clusterDis"] <- TRUE
  print(nrow(seurat@meta.data[which(seurat@meta.data$clusterDis == TRUE),]))
  seurat@meta.data$clusterDis <- factor(seurat@meta.data$clusterDis, levels = c(TRUE, FALSE))
  if(flip == TRUE){
    plot <- DimPlot(seurat, group.by = "clusterDis", pt.size = 0.01) + labs(title = numclust) + NoLegend() +
      theme(axis.title = element_text(size = 6), axis.text = element_text(size = 6), title = element_text(size = 7)) +
      scale_colour_manual(values = c("red", "grey")) +
      scale_alpha_manual(values = c(1,0)) + coord_flip() + scale_y_reverse()
  }else{
    plot <- DimPlot(seurat, group.by = "clusterDis", pt.size = 0.01) + labs(title = numclust) + NoLegend() +
      theme(axis.title = element_text(size = 6), axis.text = element_text(size = 6), title = element_text(size = 7)) +
      scale_colour_manual(values = c("red", "grey")) +
      scale_alpha_manual(values = c(1,0))
  }
  return(plot)
}

clusterDistrib_seurat_cond <- function(seurat = seurat, numclust = numclust, metaData = "numclust", flip = FALSE, color){
  print(numclust)
  seurat@meta.data$clusterDis <- "Other"
  seurat@meta.data[which(seurat@meta.data[,metaData] == numclust),"clusterDis"] <- numclust
  seurat@meta.data$clusterCondDis <- paste(seurat$clusterDis, seurat$condition, sep = "_")
  # colorClusterCond <- c(colorGenotype, alpha(colorGenotype, 0.4), rep("grey", 4))
  colorClusterCond <- c(color, rep("grey", 4))
  # print(nrow(seurat@meta.data[which(seurat@meta.data$clusterDis == TRUE),]))
  seurat@meta.data$clusterCondDis <- factor(seurat@meta.data$clusterCondDis, 
                                            levels = c(paste(numclust, "Young_Wild", sep = "_"), paste(numclust, "Young_Mutant", sep = "_"), 
                                                       paste(numclust, "Aged_Wild", sep = "_"), paste(numclust, "Aged_Mutant", sep = "_"), 
                                                       "Other_Young_Wild", "Other_Young_Mutant", "Other_Aged_Wild", "Other_Aged_Mutant"))
  names(colorClusterCond) <- levels(seurat$clusterCondDis)
  if(flip == TRUE){
    plot <- DimPlot(seurat, group.by = "clusterCondDis", pt.size = 0.1) + labs(title = numclust) +
      theme(axis.title = element_text(size = 6), axis.text = element_text(size = 6), title = element_text(size = 7)) +
      scale_colour_manual(values = colorClusterCond) +
      coord_flip() + scale_y_reverse()
  }else{
    plot <- DimPlot(seurat, group.by = "clusterCondDis", pt.size = 0.1) + labs(title = numclust) +
      theme(axis.title = element_text(size = 6), axis.text = element_text(size = 6), title = element_text(size = 7)) +
      scale_colour_manual(values = colorClusterCond) 
  }
  return(plot)
}

#Barplot to show distribution of samples along the cluster
getSamplePropPerClustBarplot <- function(hspc.combined) {
  clustersampleName <- ddply(hspc.combined@meta.data,~FinalCluster + sampleName,nrow)
  nCellSample <- ddply(hspc.combined@meta.data,~sampleName, nrow)
  names(nCellSample)[2] <- "TotalCell"
  test <- join(x = clustersampleName, y = nCellSample, by = "sampleName")
  
  test$propSample <- test$V1/test$TotalCell
  
  propExpect <- table(hspc.combined@meta.data$sampleName)/length(hspc.combined@meta.data$sampleName)[]
  propYoungExp <- propExpect[[unique(hspc.combined@meta.data$sampleName)[1]]]
  
  sampleName <- ggplot(data.frame(test), aes(fill = sampleName,y = propSample, x=FinalCluster)) +
    geom_bar( stat="identity", position="fill")+
    scale_fill_manual( values= rev(hue_pal()(length(unique(hspc.combined@meta.data$sampleName)))))+
    scale_y_continuous(name = "Sample (%)", labels = c(0,25,50,75,100))+
    ylab(label = "")+xlab(label = "") + coord_flip()+
    theme(legend.title=element_blank(), axis.text.y = element_text(size = 16)) 
  return(sampleName)
}

getGenotypePropPerClustBarplot <- function(hspc.combined) {
  clusterGenotype <- ddply(hspc.combined@meta.data,~FinalCluster + Genotype,nrow)
  
  propExpect <- table(hspc.combined@meta.data$Genotype)/length(hspc.combined@meta.data$Genotype)[]
  # propWildExp <- propExpect[[unique(hspc.combined@meta.data$Genotype)[2]]]
  propWildExp <- propExpect["Wild"][[1]]
  #clusterAGE$numclust <- factor(v$clusterNature , levels = c(""))
  clusterGenotype$Genotype <- factor(clusterGenotype$Genotype , levels = c("Mutant","Wild")) 
  
  
  Genotype <- ggplot(data.frame(clusterGenotype), aes(fill = Genotype,y = V1, x=FinalCluster)) +
    geom_bar( stat="identity", position="fill")+
    scale_fill_manual( values= rev(hue_pal()(length(unique(hspc.combined@meta.data$Genotype)))))+
    scale_y_continuous(name = "Genotype (%)", labels = c(0,25,50,75,100))+
    ylab(label = "")+xlab(label = "") + coord_flip() + geom_hline(yintercept = propWildExp)+
    theme(legend.title=element_blank(), axis.text.y = element_text(size = 16))
  return(Genotype)
  
}

getAGEPropPerClustBarplot_2 <- function(hspc.combined) {
  clusterAge <- ddply(hspc.combined@meta.data,~FinalCluster + AGE,nrow)
  
  propExpect <- table(hspc.combined@meta.data$AGE)/length(hspc.combined@meta.data$AGE)[]
  propYoungExp <- propExpect[[unique(hspc.combined@meta.data$AGE)[1]]]
  
  #clusterAGE$numclust <- factor(v$clusterNature , levels = c(""))
  clusterAge$AGE <- factor(clusterAge$AGE , levels = c("Aged","Young"))
  
  
  AGE <- ggplot(data.frame(clusterAge), aes(fill = AGE,y = V1, x=FinalCluster,levels = "Young","Aged")) +
    geom_bar( stat="identity", position="fill")+
    scale_fill_manual( values= rev(hue_pal()(length(unique(hspc.combined@meta.data$AGE)))))+
    scale_y_continuous(name = "Age (%)", labels = c(0,25,50,75,100))+
    ylab(label = "")+xlab(label = "") + coord_flip() + geom_hline(yintercept = 1-propYoungExp)+
    theme(legend.title=element_blank(), axis.text.y = element_text(size = 16))
  return(AGE)
  
}

#Bind sampleName and barcode (needed to annotate CR2 data)
tagCell <- function(cds){
  pData(cds)$taggedBarcode <- paste(pData(cds)$sampleName, pData(cds)$barcode, sep = "_")
  return(cds)
}

#Add information if the cell is found with CellRanger2
getCR2diff <- function(seurat = seurat, allCellCR2){
  uniqueCell <- setdiff(rownames(seurat@meta.data), allCellCR2)
  
  seurat@meta.data$CellRanger2 <- "Both"
  seurat@meta.data[uniqueCell,"CellRanger2"] <- "CellRanger3"
  seurat@meta.data[which(is.element(seurat@meta.data$sampleName, set = c("young_Zbtb_B", "old_Zbtb_B"))),"CellRanger2"] <- "NewData"
  
  clusterCellRanger2 <- ddply(seurat@meta.data,~FinalCluster + CellRanger2,nrow)
  
  originName <- ggplot(data.frame(clusterCellRanger2), aes(fill = CellRanger2,y = V1, x=FinalCluster)) +
    geom_bar( stat="identity", position="fill")+
    scale_fill_manual( values= rev(hue_pal()(length(unique(seurat@meta.data$CellRanger2)))))+
    scale_y_continuous(name = "Sample (%)", labels = c(0,25,50,75,100))+
    ylab(label = "")+xlab(label = "") + coord_flip()+
    theme(legend.title=element_blank()) 
  return(originName)
}

#Hypergeometric test and plotting the result with a bp function, use a funtion from funForSeurat.R
getEnrichPopClust <- function(hspc.combined, Xname, Yname, colorX, colorY, metaCol = "AGE"){
  conditionEnrich <- getEnrichAge(hspc.combined = hspc.combined,clustCol = "numclust", metaCol = metaCol)
  conditionEnrich <- as.data.frame(t(conditionEnrich))
  conditionEnrich$color <- "black"
  conditionEnrich[which(as.numeric(as.vector(conditionEnrich$phyper)) < 0.05 & conditionEnrich$enrich == Xname),"color"] <- colorX
  conditionEnrich[which(as.numeric(as.vector(conditionEnrich$phyper)) < 0.05 & conditionEnrich$enrich == Yname),"color"] <- colorY
  return(conditionEnrich)
}

# Allow you to find the markers between 2 conditions for each 
MarkersAnalysis <- function(seurat = seurat, analysis){
  seurat@meta.data$Analysis <- paste(seurat@meta.data[,analysis], seurat@meta.data$FinalCluster, sep = "_")
  condAnalysis <- unique(seurat@meta.data[,analysis])
  
  Idents(seurat) <- "Analysis"
  
  allMarkers <- data.frame()
  for(clustName in unique(seurat@meta.data$FinalCluster)){
    print(clustName)
    markersClust <- FindMarkers(object = seurat, ident.1 = paste(condAnalysis[1], clustName, sep = "_"), 
                                ident.2 = paste(condAnalysis[2], clustName, sep = "_"), logfc.threshold = 0.25)
    
    markersClust[,paste(condAnalysis[1], condAnalysis[2], sep = "_vs_")] <- paste(condAnalysis[1], clustName, "up", sep ="_")
    markersClust[which(markersClust$avg_logFC < 0),paste(condAnalysis[1], condAnalysis[2], sep = "_vs_")] <- paste(condAnalysis[1], clustName, "down",sep ="_")
    markersClust$Gene <- rownames(markersClust)
    markersClust <- markersClust[order(markersClust$avg_logFC),]
    markersClust <- markersClust[which(markersClust$p_val_adj < 0.05),]
    allMarkers <- rbind(allMarkers, markersClust)
  }
  return(allMarkers)
}

getPredictedPerSample <- function(seurat){
  propPredicted <- ddply(seurat@meta.data,~sampleName + predicted, nrow)
  predicted <- ggplot(data.frame(propPredicted), aes(fill = predicted,y = V1, x=sampleName)) +
    geom_bar( stat="identity", position="fill")+
    scale_fill_manual( values= hue_pal()(length(unique(seurat@meta.data$predicted))))+
    scale_y_continuous(name = "Cell type (%)", labels = c(0,25,50,75,100))+
    ylab(label = "")+xlab(label = "") +
    theme(legend.title=element_blank(), axis.text.x = element_text(angle = 75, hjust = 1, size = 16), axis.text.y = element_text(size = 12))
  
  return(predicted)
}

getCellTypePerSample <- function(seurat){
  propPredicted <- ddply(seurat@meta.data,~sampleName + cellType_Actinn, nrow)
  predicted <- ggplot(data.frame(propPredicted), aes(fill = cellType_Actinn,y = V1, x=sampleName)) +
    geom_bar( stat="identity", position="fill")+
    scale_fill_manual( values= hue_pal()(length(unique(seurat@meta.data$cellType_Actinn))))+
    scale_y_continuous(name = "Cell type (%)", labels = c(0,25,50,75,100))+
    ylab(label = "")+xlab(label = "") +
    theme(legend.title=element_blank(), axis.text.x = element_text(angle = 75, hjust = 1, size = 16), axis.text.y = element_text(size = 12))
  
  return(predicted)
}

getPhasesPerSample <- function(seurat){
  propPhases <- ddply(seurat@meta.data,~sampleName + phases, nrow)
  phases <- ggplot(data.frame(propPhases), aes(fill = phases,y = V1, x=sampleName)) +
    geom_bar( stat="identity", position="fill")+
    scale_fill_manual( values= hue_pal()(length(unique(seurat@meta.data$phases))))+
    scale_y_continuous(name = "Cell type (%)", labels = c(0,25,50,75,100))+
    ylab(label = "")+xlab(label = "") +
    theme(legend.title=element_blank(), axis.text.x = element_text(angle = 75, hjust = 1, size = 16),  axis.text.y = element_text(size = 12))
  
  return(phases)
}

# getPhasePropPerClustBarplotSplit <- function(hspc.combined, color = hue_pal()(length(unique(hspc.combined@meta.data$phases)))) {
#   clusterpredicted <- ddply(hspc.combined@meta.data,~numclust + phases + condition,nrow)
#   
#   split_plot <- list()
#   for(i in unique(clusterpredicted$condition)){
#     print(i)
#     predicted <- ggplot(data.frame(clusterpredicted[which(clusterpredicted$condition == i),]), aes(fill = phases,y = V1, x=numclust)) +
#       geom_bar( stat="identity", position="fill")+
#       scale_fill_manual( values= color)+
#       scale_y_continuous(name = "Cell type (%)", labels = c(0,25,50,75,100))+
#       ylab(label = "")+xlab(label = "") + coord_flip()+
#       labs(title = i) + theme(plot.title = element_text(size = 10), axis.text.x = element_text(size = 12))
#     split_plot[[i]] <- predicted
#   }
#   return(split_plot)
# }
getPhasePropPerClustBarplotSplit <- function(hspc.combined, color = hue_pal()(length(unique(hspc.combined@meta.data$phases))),
                                                  condition = "Young_Wild", sig_val = "black") {
  clusterpredicted <- ddply(hspc.combined@meta.data,~FinalCluster + phases + condition,nrow)
  
  split_plot <- list()
  for(i in condition){
    print(i)
    predicted <- ggplot(data.frame(clusterpredicted[which(clusterpredicted$condition == i),]), aes(fill = phases,y = V1, x=FinalCluster)) +
      geom_bar( stat="identity", position="fill")+
      scale_fill_manual( values= color)+
      scale_y_continuous(name = "Cell type (%)", labels = c(0,25,50,75,100))+
      ylab(label = "")+xlab(label = "") + coord_flip()+ 
      labs(title = i) + theme(plot.title = element_text(size = 10), axis.text.x = element_text(size = 12), 
                              axis.text.y = element_text(colour = sig_val, size = 16))
    split_plot[[i]] <- predicted
  }
  return(split_plot)
}

Fisher_for_cellCycle <- function(seurat, condition){
  summary_data <- ddply(seurat@meta.data,~condition + FinalCluster + phases, nrow)
  
  results_list <- list()
  for(i in levels(seurat$FinalCluster)){
    print(i)
    summ_clust <- summary_data[which(summary_data$FinalCluster == i & summary_data$condition %in% condition),]
    contingence <- data.frame(dcast(summ_clust, phases ~ condition, fill = 0, value.var = "V1"), row.names = 1)
    if(nrow(contingence)>1){
      fisher <- fisher.test(contingence)
      results_list[i] <- fisher$p.value
    } else{
      results_list[i] <- 1
    }
  }
  result <- as.data.frame(do.call("rbind",results_list))
  result$color <- "black"
  result[which(result$V1 < 0.05),"color"] <- "#DD3497"
  
  return(result)
}

Fisher_for_cellType <- function(seurat, condition){
  summary_data <- ddply(seurat@meta.data,~condition + FinalCluster + cellType_Actinn, nrow)
  
  results_list <- list()
  for(i in levels(seurat$FinalCluster)){
    print(i)
    summ_clust <- summary_data[which(summary_data$FinalCluster == i & summary_data$condition %in% condition),]
    contingence <- data.frame(dcast(summ_clust, cellType_Actinn ~ condition, fill = 0, value.var = "V1"), row.names = 1)
    if(nrow(contingence)>1){
      fisher <- fisher.test(contingence)
      results_list[i] <- fisher$p.value
    } else{
      results_list[i] <- 1
    }
  }
  result <- as.data.frame(do.call("rbind",results_list))
  result$color <- "black"
  result[which(result$V1 < 0.05),"color"] <- "red"
  
  return(result)
}

Fisher_for_predicted <- function(seurat, condition){
  summary_data <- ddply(seurat@meta.data,~condition + FinalCluster + predicted, nrow)
  
  results_list <- list()
  for(i in levels(seurat$FinalCluster)){
    print(i)
    summ_clust <- summary_data[which(summary_data$FinalCluster == i & summary_data$condition %in% condition),]
    contingence <- data.frame(dcast(summ_clust, predicted ~ condition, fill = 0, value.var = "V1"), row.names = 1)
    if(nrow(contingence)>1){
      fisher <- fisher.test(contingence)
      results_list[i] <- fisher$p.value
    } else{
      results_list[i] <- 1
    }
  }
  result <- as.data.frame(do.call("rbind",results_list))
  result$color <- "black"
  result[which(result$V1 < 0.05),"color"] <- "red"
  
  return(result)
}


getPredictedPropPerClustBarplotSplit <- function(hspc.combined, colors = hue_pal()(length(unique(hspc.combined@meta.data$predicted))), 
                                                 condition = "Young_Wild", sig_val = "black") {
  clusterpredicted <- ddply(hspc.combined@meta.data,~FinalCluster + predicted + condition,nrow)
  
  split_plot <- list()
  for(i in condition){
    print(i)
    predicted <- ggplot(data.frame(clusterpredicted[which(clusterpredicted$condition == i),]), aes(fill = predicted,y = V1, x=FinalCluster)) +
      geom_bar( stat="identity", position="fill")+
      scale_fill_manual( values= colors)+
      scale_y_continuous(name = "Cell type (%)", labels = c(0,25,50,75,100))+
      ylab(label = "")+xlab(label = "") + coord_flip()+
      labs(title = i) + theme(plot.title = element_text(size = 10), axis.text.x = element_text(size = 12), 
                              axis.text.y = element_text(colour = sig_val, size = 16))
    split_plot[[i]] <- predicted
  }
  return(split_plot)
}

getCellTypePropPerClustBarplotSplit <- function(hspc.combined, colors = hue_pal()(length(unique(hspc.combined@meta.data$predicted))), 
                                                 condition = "Young_Wild", sig_val = "black") {
  clusterpredicted <- ddply(hspc.combined@meta.data,~FinalCluster + cellType_Actinn + condition,nrow)
  
  split_plot <- list()
  for(i in condition){
    print(i)
    cellType <- ggplot(data.frame(clusterpredicted[which(clusterpredicted$condition == i),]), aes(fill = cellType_Actinn,y = V1, x=FinalCluster)) +
      geom_bar( stat="identity", position="fill")+
      scale_fill_manual( values= colors)+
      scale_y_continuous(name = "Cell type (%)", labels = c(0,25,50,75,100))+
      ylab(label = "")+xlab(label = "") + coord_flip()+
      labs(title = i) + theme(plot.title = element_text(size = 10),  
                              axis.text.y = element_text(colour = sig_val, size = 16))
    split_plot[[i]] <- cellType
  }
  return(split_plot)
}


getPhasePropPerPredictedBarplotSplit <- function(hspc.combined, colors = hue_pal()(length(unique(hspc.combined@meta.data$phases))), 
                                                 condition = "Young_Wild") {
  clusterpredicted <- ddply(hspc.combined@meta.data,~phases + predicted + condition,nrow)
  
  split_plot <- list()
  for(i in condition){
    print(i)
    predicted <- ggplot(data.frame(clusterpredicted[which(clusterpredicted$condition == i),]), aes(fill = phases,y = V1, x=predicted)) +
      geom_bar( stat="identity", position="fill")+
      scale_fill_manual( values= colors)+
      scale_y_continuous(name = "Cell type (%)", labels = c(0,25,50,75,100))+
      ylab(label = "")+xlab(label = "") + 
      labs(title = i) + theme(plot.title = element_text(size = 10))
    split_plot[[i]] <- predicted
  }
  return(split_plot)
}

# getCellTypePropPerClustBarplotSplit <- function(hspc.combined, colors = hue_pal()(length(unique(hspc.combined@meta.data$cellType_Actinn)))) {
#   clusterpredicted <- ddply(hspc.combined@meta.data,~FinalCluster + cellType_Actinn + condition,nrow)
#   
#   split_plot <- list()
#   for(i in unique(clusterpredicted$condition)){
#     print(i)
#     predicted <- ggplot(data.frame(clusterpredicted[which(clusterpredicted$condition == i),]), aes(fill = cellType_Actinn,y = V1, x=FinalCluster)) +
#       geom_bar( stat="identity", position="fill")+
#       scale_fill_manual( values= colors)+
#       scale_y_continuous(name = "Cell type (%)", labels = c(0,25,50,75,100))+
#       ylab(label = "")+xlab(label = "") + coord_flip()+
#       labs(title = i) + theme(plot.title = element_text(size = 10))
#     split_plot[[i]] <- predicted
#   }
#   return(split_plot)
# }

getCellTypePerCondition <- function(seurat, colors = hue_pal()(length(unique(seurat$cellType_Actinn)))){
  cellType <- ddply(seurat@meta.data,~condition + cellType_Actinn, nrow)
  barplot <- ggplot(data = data.frame(cellType), aes(fill=cellType_Actinn, y=V1, x=condition)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = colors) + 
    scale_y_continuous(name = "Cell type (%)", labels = c(0,25,50,75,100)) +
    ggtitle("Actinn") + theme(axis.text.x = element_text(size = 16, angle = -15), axis.title.x = element_blank(), plot.title = element_text(size = 25))
  
  return(barplot)
}
