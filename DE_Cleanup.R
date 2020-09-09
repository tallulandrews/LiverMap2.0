files <- Sys.glob("*SEX-F*DE*.csv")

my_de_list <- vector();

for(f in files) {
	de_mat <- read.delim(f, sep=",", header=T)
	if (ncol(de_mat) ==8) {
		rownames(de_mat) <- as.character(de_mat[,1])
		de_mat <- de_mat[,-1]
	}
	colnames(de_mat) <- c("FC", "Mean", "pval", "qval", "dir", "desc", "known")
	de_mat <- de_mat[!grepl(":", rownames(de_mat)),]
	#de_mat$qval <- signif(de_mat$qval, digits=1)
	#de_mat$score <- round(-log10(de_mat$qval))*de_mat$dir
	#de_mat$score_undir <- round(-log10(de_mat$qval))
	de_mat$score <- de_mat$FC
	de_mat$score[de_mat$qval > 0.05] <- 0;
	de_mat$rank <- rank(de_mat$score)
	out <- de_mat$score
	names(out) <- rownames(de_mat);
	
	if (length(my_de_list) == 0) {
		my_de_list <- matrix(out, ncol=1)
		rownames(my_de_list) <- names(out)
		colnames(my_de_list) <- f
	} else {
		all_genes <- sort(unique(c(names(out), rownames(my_de_list))))
		out <- out[match(all_genes, names(out))]
		out[is.na(out)] <- 0
		my_de_list <- my_de_list[match(all_genes, rownames(my_de_list)),]
		my_de_list[is.na(my_de_list)] <- 0
		my_de_list <- cbind(my_de_list, out);
		rownames(my_de_list) <- all_genes;
	}
}

max_score <- max(abs(my_de_list[is.finite(my_de_list)]))
my_de_list[my_de_list > max_score] <- max_score+1
my_de_list[my_de_list < -1*max_score] <- -1*max_score-1

rank_mat <- apply(my_de_list, 2, rank)

scaled_rank_mat <- t(apply(rank_mat, 1, function(x){x/mean(x)}))
	
scaled_score_mat <- t(apply(my_de_list, 1, function(x){
					y<-x[x!=0]; if (length(y)==0) {return(rep(1, length(x)))} 
							else { x/(mean(abs(y)))}}))

tmp <- scaled_score_mat[,4]; tmp <- sort(tmp); tmp

#GSEA or other Pathway tools.

#gsea_input <- data.frame(NAME=rownames(scaled_score_mat), score=scaled_score_mat[,4])
#GSEA(gsea_input, c(1,2), gs.db="../../ExternalData/ReactomePathways.gmt")

require(fgsea)
gmt <- read.delim("../../ExternalData/ReactomePathways.gmt", sep="\t", header=F)

out <- apply(gmt, 1, function(x){name=paste(x[2], gsub(" ", "_", x[1]), sep="_"); 
				  genes <- x[3:length(x)]; genes <- genes[genes != ""]; 
				  l <- list(); l[[name]] <- genes; return(l)})
pathways <- list()
for (p in out) {
	tmp <- unlist(p)
	names(tmp) <- NULL
	pathways[[names(p)]] <- tmp
}

expr_score <- scaled_score_mat[,4]
expr_score <- expr_score[expr_score != 0]
out <- fgsea(pathways=pathways, stats=expr_score, minSize=15, maxSize=500)

up <- out[out$ES > 0]
up <- up[order(up$pval),]
head(up)

require(ggplot2)
plotEnrichment(pathways[["R-HSA-449147_Signaling_by_Interleukins"]], expr_score)+ labs(title="Signaling by Interleukins")