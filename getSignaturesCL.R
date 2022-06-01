##------------------------------------------------------------------------------
## L. Herault
##------------------------------------------------------------------------------

## -----------------------------------------------------------------------------
## Libraries
## -----------------------------------------------------------------------------

suppressMessages(library(getopt))
suppressMessages(library(readxl))


source("R_src/getSigPerCells.R")



# Store previously published signature in an R object

## -----------------------------------------------------------------------------
## Command line args
## -----------------------------------------------------------------------------

spec = matrix(c(
  'help',        'h', 0, "logical",	"Help about the program",
  'input_tabfile', 't', 1, "character","tab file with gene lists per column (Geiger signatures tab file)",
  'input_microArraySig',"m", 1, "character", "directory  with micro array cell type signatures (Chambers) or directly xlsx file from cell paper",
  'outdir',      'o', 1, "character", 	"Output directory. Default to current working directory."),
  byrow=TRUE, ncol=5)

opt = getopt(spec)

# if help was asked, print a friendly message
# and exit with a non-zero error code
args <- commandArgs()
if ( !is.null(opt$help) | is.null(opt$input_tabfile) | is.null(opt$input_microArraySig)) {
  cat("get signatures with microo array data and/or tab files from previous studies\n")
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}



sessionInfo()

# output
if(is.null(opt$outdir))
  opt$outdir <- getwd()

dir.create(opt$outdir, showWarnings=FALSE)

#Miccro array sig

if (dir.exists(opt$input_microArraySig)) {
  arraySigFiles <- list.files(opt$input_microArraySig,full.names = F)
  print(arraySigFiles)
  arraySig <- getMicroArraySigList(arraySigFiles)
  
  
  
  allSignatures <- c(arraySig)
  
  allSignatures <- lapply(allSignatures,'[[',2)
} else {
  
  # Miccro array sig directly from xlsx file
  allsheets <- read_excel_allsheets(opt$input_microArraySig)
  
  allSignatures <- getMicroArraySigListXls(allsheets)
  
}

#Tab file signatures

if(endsWith(opt$input_tabfile,suffix = ".xls")) { #if downloaded file
  tab_sig <- read_xls(opt$input_tabfile)
  tab_sig <- as.data.frame(tab_sig)
  tab_sig <- tab_sig[-2,]
  colnames(tab_sig) <- tab_sig[1,]
  tab_sig <- tab_sig[-1,]
} else {
  tab_sig <- read.table(opt$input_tabfile,sep = "\t",header=T)
  
}


nonEmptyString <- function(stringVector){
  result <- which(stringVector != "")
  return(stringVector[result])
}

list_sig <- as.list(tab_sig)

list_sig <- lapply(list_sig,nonEmptyString)


allSignatures <- c(allSignatures,list_sig)



saveRDS(allSignatures,paste(opt$outdir,"publicSignatures.rds",sep = ""))

