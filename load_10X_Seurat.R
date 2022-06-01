##------------------------------------------------------------------------------
## B. NDAO
##------------------------------------------------------------------------------

## SET DIRECTORY
#setwd("/shared/projects/scrna_gmp_rara/stage_babacar/scranseq/analysis/")

## -----------------------------------------------------------------------------
## Libraries
## -----------------------------------------------------------------------------


suppressMessages(library("getopt"))
suppressMessages(library("Matrix"))
suppressMessages(library("Seurat"))

source("/shared/projects/scrna_gmp_rara/stage_babacar/scranseq/analysis/R_src/funForLoading.R")


## -----------------------------------------------------------------------------
## Command line args
## -----------------------------------------------------------------------------

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputMatrixDir', 'b', 1, "character", "REQUIRED:10X expression data (matrix files directory path)",
  'outfile',     'f',1, "character", 'Output file path (default .10X_data.rds)',
  'cellranger',   'c',1, "character", "cellranger version used to generate matrix files (default 2)",
  "subSample",    "p",1, "character", "proportion of cell to subsample (by default all the cells are used (e.g. 0.5))",
  "sampleInfo",   "s", 1, "character", "sample information that will be added to pData with the following format: age=2_months,runDate=10_12_2017..",
  "sampleName",   "n", 1, "character", "if matrix file have the sampl name as prefix (name_matrix.mtx...) provide it here"
  
), byrow=TRUE, ncol=5);

opt = getopt(spec)

# if help was asked, print a friendly message
# and exit with a non-zero error code
# For test
# opt <- list()
opt$inputMatrixDir <- "/shared/projects/scrna_gmp_rara/stage_babacar/scranseq/analysis/input/RAW/PLZF_RARA_RA/"
opt$cellranger <- 4
opt$sampleInfo <- "age=Adult,runDate=11/2020,sampleName=RNA_RA"
opt$outfile <- "/shared/projects/scrna_gmp_rara/stage_babacar/scranseq/analysis/input/10X_data_RNA_RA.rds"

args <- commandArgs()
if ( !is.null(opt$help) | is.null(opt$inputMatrixDir)) {
  cat("loading 10X genomics single cells data and save it as a seurat R object.")
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

## -----------------------------------------------------------------------------
## Processing data
## -----------------------------------------------------------------------------

# output

if(is.null(opt$outfile)) {
  opt$outfile <- "10X_data.rds"
  outdir = './'
} else {
  outdir <- paste(strsplit(opt$outfile,"/")[[1]][-length(strsplit(opt$outfile,"/")[[1]])],collapse="/")
  dir.create(outdir, showWarnings=FALSE)
}


if(is.null(opt$cellranger)) {
  opt$cellranger <- "2"
}


if (opt$cellranger == "3" | opt$cellranger == "4") {
  gbm <- loadCellRangerMatrix_cellranger3(opt$inputMatrixDir,sample_name=opt$sampleName)
} else
{
  gbm <- loadCellRangerMatrix(opt$inputMatrixDir,sample_name=opt$sampleName)
}


fd <- gbm$fd
pd <- gbm$pd

# Column 'symbol' is the one (from cellRanger workflow) that corresponds to featureData's gene short names.

colnames(fd)[which(colnames(fd)=="symbol")]<- "gene_short_name"

for (c in c(1:length(colnames(fd)))){
  if (is.na(colnames(fd)[c])) {
    colnames(fd)[c] <- paste0("col_",c)
  }
}


seurat <- CreateSeuratObject(counts = gbm$exprs,
                             assay = "RNA",
                             meta.data = gbm$pd)

# ## Subsample ## DEPRECATED
# if(is.null(opt$subSample)==F) {
#   opt$subSample <- as.numeric(opt$subSample)
#   print(opt$subSample)
#   
#   cellSubset <- sample(rownames(pData(gbm_cds)),size = opt$subSample*length(rownames(pData(gbm_cds))))
#   gbm_cds <- gbm_cds[,cellSubset] 
# }


#Add sample infos

print(opt$sampleInfo)
sampleInfos <- strsplit(opt$sampleInfo,split=",")[[1]]
for (i in sampleInfos) {
  print(i)
  info <- strsplit(i,split = "=")[[1]][1]
  print(info)
  value <- strsplit(i,split = "=")[[1]][2]
  print(value)
  seurat@meta.data[,info] <- value
}


output <- list(seurat = seurat,featureData = fd)

saveRDS(output,file = opt$outfile)


