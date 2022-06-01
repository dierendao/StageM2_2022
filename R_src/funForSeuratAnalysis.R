# Do hypergeometric test of to population to find if ine proportion is significantly increased
getEnrichAge <- function(hspc.combined,clustCol ='clusterName',metaCol = "age") {
  
  table <- table(hspc.combined@meta.data[,metaCol],hspc.combined@meta.data[,clustCol])
  
  #Remove null column in case of reclustering has been made
  
  table <- table[,as.vector(which(colSums(table)>0))]
  
  tablePercent <- prop.table(table,2)
  
  propExpect <- table(hspc.combined@meta.data[,metaCol])/length(hspc.combined@meta.data[,metaCol])
  propExpectAge_1<- propExpect[[unique(hspc.combined@meta.data[,metaCol])[1]]]
  propExpectAge_2<- propExpect[[unique(hspc.combined@meta.data[,metaCol])[2]]]
  phyper <- rep(NA,length(colnames(table)))
  enrich <- rep(NA,length(colnames(table)))
  tablePercent <- rbind(tablePercent,enrich,phyper)
  
  
  for (age in unique(hspc.combined@meta.data[,metaCol])) {
    for (cluster in colnames(table)) {
      if(tablePercent[age,cluster] > propExpect[[age]]) {
        cells_pull_marked <- table[age,as.character(cluster)]
        cells_pull <- as.numeric(colSums(table)[as.character(cluster)])
        cells_marked_all <- rowSums(table)[age]
        all_cells <- length(hspc.combined@meta.data[,metaCol])
        
        
        
        p.value <-  phyper(q=cells_pull_marked -1, 
                           m=cells_marked_all,
                           n=all_cells - cells_marked_all, k= cells_pull, lower.tail=FALSE)
        
        tablePercent["enrich",cluster] <- age
        
        tablePercent["phyper",cluster] <- p.value
        
      }
    }
  }
  return(tablePercent)
  
}

# Hypergeometric test and plotting the result with a bp function, use a getEnrichAge
getEnrichPopClust <- function(hspc.combined, Xname, Yname, colorX, colorY, metaCol = "AGE", clustCol = "numclust"){
  conditionEnrich <- getEnrichAge(hspc.combined = hspc.combined,clustCol = clustCol, metaCol = metaCol)
  conditionEnrich <- as.data.frame(t(conditionEnrich))
  conditionEnrich$color <- "black"
  conditionEnrich[which(as.numeric(as.vector(conditionEnrich$phyper)) < 0.01 & conditionEnrich$enrich == Xname),"color"] <- colorX
  conditionEnrich[which(as.numeric(as.vector(conditionEnrich$phyper)) < 0.01 & conditionEnrich$enrich == Yname),"color"] <- colorY
  return(conditionEnrich)
}

# Function to calculate the score diff of a signature between 2 cells pop
featureDiff <- function(seurat,cells.1,cells.2,feature) {
  data <- GetAssayData(seurat,slot = "data")
  total.diff <- mean(data[feature,cells.1]) - mean(data[feature,cells.2])
  return(total.diff)
} 

convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
