
caudate_freq <- table(merge_categories(obj_list$Caudate@meta.data$scmap_cell_anno))

flush_freq <- table(merge_categories(obj_list$Flush@meta.data$scmap_cell_anno))

type_freq_table <- cbind(caudate_freq/sum(caudate_freq)*100, flush_freq/sum(flush_freq)*100)
type_freq_table <- type_freq_table[apply(type_freq_table,1,max)>1,]
colnames(type_freq_table) <- c("Caudate", "Flush")
type_freq_table <- type_freq_table[order(apply(type_freq_table,1,max)),]

# Significance:
# Start from biggest difference to smallest
# if Sig remove all those cells from both datasets before continuing
# -> avoid sig because of proportionality b/c of other group that is sig

diff <- abs(type_freq_table[,1]-type_freq_table[,2])
p <- list()
for (type in names(diff)[order(diff, decreasing=T)]) {
	a = caudate_freq[type]
	b = sum(caudate_freq[names(caudate_freq) != type])
	d = flush_freq[type]
	e = sum(flush_freq[names(flush_freq) != type])
	
	out <- fisher.test(cbind(c(a,b),c(d,e)))
	p[[type]] <- out$p.value
	if (p[[type]] < 0.05/length(diff)){
		caudate_freq <- caudate_freq[names(caudate_freq) != type]
		flush_freq <- flush_freq[names(flush_freq) != type]
	}

}


png("CF_Cell_type_freqs.png", width=7, height=6, units="in", res=300)
par(mar=c(3,6,1,1))
loc <- barplot(t(type_freq_table), beside=T, las=1, horiz=T, xlab="% of cells", col=c("firebrick", "grey85"))
mtext("% of cells", 1,2)
sig <- names(p)[p< 0.01]
colnames(loc) <- rownames(type_freq_table)
xes <- colMeans(loc)
yes <- apply(type_freq_table, 1, max)
p_vals <- signif(unlist(p)[match(names(xes), names(p))], digits=1)
labs <- as.character(p_vals)
labs[p_vals > 0.01] <- ""
text(yes, xes, labs, pos=4)
legend("bottomright", c("Flush", "Caudate"), bty="n", fill=c("grey85", "firebrick"))
dev.off()


