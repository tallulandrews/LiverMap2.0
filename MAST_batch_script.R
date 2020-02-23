args <- commandArgs(trailingOnly=TRUE);
this_contrast <- args[1];

my_rowMeans <- function(x) {
	if (!is.null(ncol(x))) {
		if (ncol(x) > 1) {
			return(Matrix::rowMeans(x))
		}
	}
	return(x);
}
my_rowSums <- function(x) {
	if (!is.null(ncol(x))) {
		if (ncol(x) > 1) {
			return(Matrix::rowSums(x))
		}
	}
	return(x);
}

## Significant with MAST
set.seed(8817)

require("MAST")

mod <- readRDS("harmony_MASTmodel.rds")

set.seed(3891)

res1 <- summary(mod, doLRT=this_contrast)
res <- res1$datatable[res1$datatable$contrast==this_contrast & res1$datatable$component %in% c("C","D"),]
res$fdr <- p.adjust(unlist(res[,4]), method="fdr")

con <- res[res$component =="C",]
dis <- res[res$component =="D",]

sig <- (con$fdr < 0.01 | dis$fdr < 0.01) & sign(con$coef) == sign(dis$coef)
out <- cbind(con[sig,], dis[sig,])
out$overall <- out[,7]+out[,16]
out <- out[order(out$overall, decreasing=T),]
write.table(out, file=paste(this_contrast,"_MAST_output.txt", sep=""), col.names=T, row.names=F)
 
