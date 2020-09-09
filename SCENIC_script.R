require(SCENIC)
## Set up SCENIC parameters
org="hgnc" #mgi or hgnc, or dmel
dbDir="cisTarget_databases" # RcisTarget databases location
myDatasetTitle="SCENIC Human Liver" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)
# Modify if needed
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
# Databases:
# scenicOptions@settings$dbs <- c("mm9-5kb-mc8nr"="mm9-tss-centered-5kb-10species.mc8nr.feather")
# scenicOptions@settings$db_mcVersion <- "v8"

# Save to use at a later time...
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 


## Read in raw data
obj <- readRDS("/cluster/projects/macparland/TA/LiverMap2.0/RawData/Analysis/EmptyDrops/SoupX/C41_emptyDrops_table_SeurObj_SoupX.rds")

exprMat <- obj@assays$RNA@counts

genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)

exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)

rm(exprMat)

runCorrelation(exprMat_filtered, scenicOptions)

### This is extremely computationally expensive....

# Optional: add log (if it is not logged/normalized already)
exprMat_filtered <- log2(exprMat_filtered+1) 

# Run GENIE3
runGenie3(exprMat_filtered, scenicOptions)

## Reload data ##
loom <- open_loom(loomPath, mode="r")
exprMat <- get_dgem(loom)
close_loom(loom)
# Optional: log expression (for TF expression plot, it does not affect any other calculation)
logMat <- log2(exprMat+1)
dim(exprMat)

library(SCENIC)
scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123

# For a very quick run: 
# coexMethod=c("top5perTarget")
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # For toy run
# save...

runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) #** Only for toy run!!
runSCENIC_3_scoreCells(scenicOptions, logMat)


### Concerns:
# 1) How to deal with batch effects between individuals / samples?
# 2) How to balance different cell-type / sub-type resolution?
# 3) How to calculate correlations efficiently?
