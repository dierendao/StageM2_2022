loadCellRangerMatrix <- function(matrix_dir,sample_name = NULL) {
  if (!is.null(sample_name)) {
    sample_name <- paste0(sample_name,"_")
  }
  barcode.path <- paste0(matrix_dir,"/",sample_name,"barcodes.tsv")
  features.path <- paste0(matrix_dir,"/",sample_name,"genes.tsv")
  matrix.path <- paste0(matrix_dir,"/",sample_name,"matrix.mtx")
  mat <- readMM(file = matrix.path)
  feature.names = read.delim(features.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(mat) = barcode.names$V1
  rownames(mat) = feature.names$V1
  colnames(feature.names) <- c("id", "symbol")
  rownames(feature.names) <- feature.names$id
  colnames(barcode.names) <- c("barcode")
  rownames(barcode.names) <- barcode.names$barcode
  results <- list()
  results$exprs <- mat
  results$fd <- feature.names
  results$pd <- barcode.names
  
  return(results)
}

loadCellRangerMatrix_cellranger3 <- function(matrix_dir,sample_name =NULL) {
  if (!is.null(sample_name)) {
    sample_name <- paste0(sample_name,"_")
  }
  barcode.path <- paste0(matrix_dir,"/",sample_name,"barcodes.tsv")
  features.path <- paste0(matrix_dir,"/",sample_name,"features.tsv")
  matrix.path <- paste0(matrix_dir,"/",sample_name,"matrix.mtx")
  mat <- readMM(file = matrix.path)
  feature.names = read.delim(features.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(mat) = barcode.names$V1
  rownames(mat) = feature.names$V1
  colnames(feature.names) <- c("id", "symbol")
  rownames(feature.names) <- feature.names$id
  colnames(barcode.names) <- c("barcode")
  rownames(barcode.names) <- barcode.names$barcode
  results <- list()
  results$exprs <- mat
  results$fd <- feature.names
  results$pd <- barcode.names

  return(results)
}
