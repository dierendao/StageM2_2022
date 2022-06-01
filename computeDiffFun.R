#Functions to compute difference of feature in seurat object (because Seurat compute logFC)

getExtAvgDiff <- function(row) {
  res <- NA
  if(row[1] > 0 & row[2] > 0) {
    res <- min(c(row))
  }
  if(row[1] < 0 & row[2] < 0) {
    res <- max(c(row))
  }
  return(res)
}

featureDiff <- function(seurat,cells.1,cells.2,feature) {
  data <- GetAssayData(seurat,slot = "data")
  total.diff <- mean(data[feature,cells.1]) - mean(data[feature,cells.2])
  return(total.diff)
} 


getTrueDiff <- function (seurat,table,colIdent = "numclust",suffix = "",colTest="AGE") {
  
  if(!is.null(levels(seurat[[colTest]]))) {
    factors <- levels(seurat@meta.data[,colTest])
  } else {
    factors <- unique(seurat@meta.data[,colTest])
  }
  factors <- as.vector(factors)

  table[,paste("avg_diff",suffix,sep = "")] <- NA
  
  for (rm in rownames(table)) {
    feature <- table[rm,"Gene"]
    cells.1 <- colnames(seurat)[which(seurat[[colIdent]] == table[rm,colIdent] & seurat[[colTest]] == factors[1])]
    cells.2 <- colnames(seurat)[which(seurat[[colIdent]] == table[rm,colIdent] & seurat[[colTest]] == factors[2])]
    table[rm,paste("avg_diff",suffix,sep = "")] <- featureDiff(seurat,cells.1,cells.2,feature)
  }
  
  return(table)
}


## Test functions

FindMarkerPerClustGroupVar <- function(cluster,
                                       hspc.combined,
                                       condition ="AGE",        ## only var with two factors allowed
                                       grouping.var = "platform",
                                       max.cells.per.ident = Inf,
                                       filterOnPadj = T,
                                       logfc.threshold = 0.25,
                                       min.pct = 0.1,
                                       keepDiverging = F,
                                       test.use = "wilcox",
                                       identCol = "numclust",
                                       pseudocount.use = 1,
                                       computeTrueDiff= F) {
  
  clusterCondition <- paste0("cluster.",condition)
  if(!is.null(levels(seurat@meta.data[,condition]))) {
    factors <- levels(seurat@meta.data[,condition])
  } else {
    factors <- unique(seurat@meta.data[,condition])
  }
  
  hspc.combined@meta.data[,clusterCondition] <- paste(hspc.combined@meta.data[,identCol],
                                                      hspc.combined@meta.data[,condition], sep = "_")
  
  test <- table(hspc.combined@meta.data[,clusterCondition], hspc.combined@meta.data[,grouping.var])
  groups = unique(hspc.combined@meta.data[,grouping.var])
  
  ident1 = paste0(cluster,"_",factors[1])
  ident2 = paste0(cluster,"_",factors[2])
  
  
  if( (test[ident1,groups[1]]>3 & test[ident1,groups[2]]>3) &
      (test[ident2,groups[1]]>3&test[ident2,groups[2]]>3)) {
    
    Idents(object = hspc.combined) <- clusterCondition
    
    
    
    
    markers <- FindConservedMarkers(hspc.combined,
                                    assay = DefaultAssay(hspc.combined),
                                    grouping.var = grouping.var,
                                    test.use=test.use,
                                    ident.1 = ident1,
                                    ident.2 = ident2,
                                    pseudocount.use = pseudocount.use,
                                    min.pct = min.pct,
                                    logfc.threshold = logfc.threshold,
                                    max.cells.per.ident = max.cells.per.ident)
    
    print(dim(markers))
    #Add a column the min abs(logFC) 
    markers$min_avg_logFC <- NA
    
    if(!keepDiverging) {
      #print("keep only markers with conserved logfc sign")
      markers <- markers[which(markers[,paste0(groups[1],"_avg_logFC")]/markers[,paste0(group[2],"_avg_logFC")] > 0),]
    } else {
      #print("combined p_val of diverging markers between two batch are set to 1 and their min_avg_logfc to 0")
      markers[which(markers[,paste0(groups[1],"_avg_logFC")]/markers[,paste0(group[2],"_avg_logFC")] > 0),"min_avg_logFC"] <- 0
      markers[which(markers[,paste0(groups[1],"_avg_logFC")]/markers[,paste0(group[2],"_avg_logFC")] > 0),"minimump_p_val"] <- 1
    }
    
    for (g in rownames(markers)) {
      if(markers[g,paste0(groups[1],"_avg_logFC")] < 0) {
        markers[g,"min_avg_logFC"] <- max(markers[g,c(paste0(groups[1],"_avg_logFC"),paste0(groups[2],"_avg_logFC"))])
      } else {
        markers[g,"min_avg_logFC"] <- min(markers[g,c(paste0(groups[1],"_avg_logFC"),paste0(groups[2],"_avg_logFC"))])
      }
    }
    
    markers$Cluster <- paste0(ident1,"_up")
    markers$Cluster[which(markers$min_avg_logFC < 0)] <-  paste0(ident1,"_down")
    markers$Gene <- rownames(markers)
    markers <- markers[order(markers$min_avg_logFC),]
    
    if(computeTrueDiff) {
      markers$numclust <- str_split_fixed(markers$Cluster,pattern = "_",n=3)[,1]
      
      
      markers <- getTrueDiff(seurat[,which(seurat[[grouping.var]] == group[1])],
                                  markers,
                                  colIdent = "numclust",
                                  colTest = condition,
                                  suffix = paste0("_",group[1]))
      
      markers <- getTrueDiff(seurat[,which(seurat[[grouping.var]] == group[2])],
                                  markers,
                                  colIdent = "numclust",
                                  colTest = condition,
                                  suffix = paste0("_",group[2]))
      
      
      markers$min_avg_diff <- apply(markers[,c(paste0("avg_diff_",group[1]),paste0("avg_diff_",group[2]))],
                                    1,
                                    FUN=getExtAvgDiff) 
    }
    
    
    if (filterOnPadj) {
      markers <- markers[which(markers$A_p_val_adj < 0.05 & markers$B_p_val_adj < 0.05),]
    }
  }else {
    markers <- NULL 
  }
  return(markers)
  
}

FindMarkerPerClust <- function(cluster,
                               hspc.combined,
                               condition ="AGE",        ## only var with two factors allowed
                               max.cells.per.ident = Inf,
                               logfc.threshold = 0.25,
                               min.pct = 0.1,
                               test.use = "wilcox",
                               identCol = "numclust",
                               pseudocount.use = 1,
                               computeTrueDiff = F) {
  
  clusterCondition <- paste0("cluster.",condition)
  if(!is.null(levels(seurat@meta.data[,condition]))) {
    factors <- levels(seurat@meta.data[,condition])
  } else {
    factors <- unique(seurat@meta.data[,condition])
  }
  
  hspc.combined@meta.data[,clusterCondition] <- paste(hspc.combined@meta.data[,identCol],
                                                      hspc.combined@meta.data[,condition], sep = "_")
  
  test <- table(hspc.combined@meta.data[,clusterCondition])
  
  ident1 = paste0(cluster,"_",factors[1])
  ident2 = paste0(cluster,"_",factors[2])
  
  
  if( test[ident1]>3 & test[ident2]>3) {
    
    Idents(object = hspc.combined) <- clusterCondition
    
    markers <- FindMarkers(hspc.combined,
                           assay = DefaultAssay(hspc.combined),
                           test.use=test.use,
                           ident.1 = ident1,
                           ident.2 = ident2,
                           pseudocount.use = pseudocount.use,
                           min.pct = min.pct,
                           logfc.threshold = logfc.threshold,
                           max.cells.per.ident = max.cells.per.ident)
    
    markers$Cluster <- paste0(ident1,"_up")
    markers$Cluster[which(markers$avg_logFC < 0)] <-  paste0(ident1,"_down")
    markers$Gene <- rownames(markers)
    markers <- markers[order(markers$avg_logFC),]
    
    if(computeTrueDiff) {
      markers$numclust <- str_split_fixed(markers$Cluster,pattern = "_",n=3)[,1]
      
      
      markers <- getTrueDiff(seurat,
                             markers,
                             colIdent = "numclust",
                             colTest = condition
                             )
      
     }
    
  } else {
    markers <- NULL 
  }
  return(markers)
  
}

