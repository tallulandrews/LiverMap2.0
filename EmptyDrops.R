##### NOTE ####
#
# This script is specifically designed for datasets where an excessive number of cells were 
# called using cellranger's pipeline. DO NOT USE on datasets where cell ranger identified 
# fewer than 20,000 cells!!!
#
##############

require(methods)
require(Matrix)
require("DropletUtils")
require(dplyr)
require(Seurat)
script_dir = "/cluster/home/tandrews/scripts/LiverMap2.0"

i <- as.numeric(as.character(commandArgs(trailingOnly=TRUE)))
#i = 1



#mt_filter <- Params[i,"MTfilter"]
#nc_filter <- Params[i,"nCellfilter"]
#this_seed <- Params[i, "Seed"]
#name <- paste("EmptyDrops", as.character(Params[i, "Name"]), sep="_")
#folder <- as.character(Params[i, "Directory_UHN"])
folder <- "../McGilvray_Sonya__C-62_Enriched_5pr/raw_feature_bc_matrix/"
name <- "EmptyDrops_C62_5pr";

set.seed(101)

n_max_cells <- 20000
n_min_cells <- 100


mydata <- Read10X(data.dir = paste(folder, sep="/"))
mydata <- mydata[,Matrix::colSums(mydata) > 0]
mydata <- mydata[Matrix::rowSums(mydata) > 0,]
print(paste("Stats out:", dim(mydata), "input dimensions of", name))
br.out <- barcodeRanks(mydata)

plausible_cell_threshold <- max(br.out$total[br.out$rank > n_max_cells]);
mandatory_cell_threshold <- max(br.out$total[br.out$rank < n_min_cells]);

#plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
#o <- order(br.out$rank)
#lines(br.out$rank[o], br.out$fitted[o], col="red")

#abline(h=br.out$knee, col="dodgerblue", lty=2)
#abline(h=br.out$inflection, col="forestgreen", lty=2)
#legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
#    legend=c("knee", "inflection"))


set.seed(100)
e.out <- emptyDrops(mydata, lower=max(plausible_cell_threshold, min(br.out$knee, br.out$inflection)), niters=100000, ignore=3, retain=mandatory_cell_threshold)

e.out <- e.out[!is.na(e.out$PValue),]
e.out$q.value <- e.out$PValue;
e.out$q.value[e.out$Limited] <- 0;
e.out$q.value <- p.adjust(e.out$q.value, method="fdr");

mydata <- mydata[,match(rownames(e.out),colnames(mydata))]
e.out$ngenes <- Matrix::colSums(mydata>0);
e.out$q.value[e.out$Limited] <- 0;


is.cell <- e.out$q.value <= 0.01
sum(is.cell, na.rm=TRUE)

png(paste(name, "emptyDrops.png", sep="_"), width=6, height=6, units="in", res=150)
plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
    xlab="Total UMI count", ylab="-Log Probability")
dev.off()

outdata <- mydata[,is.cell]
print(paste("Stats out:", dim(outdata), "input dimensions of", name))
saveRDS(outdata, paste(name, "emptyDrops_table.rds", sep="_"));
