
metadata <- read.delim("Metadata20LiverMapPlusParams.csv", sep=",", header=T)
metadata <- metadata[1:(nrow(metadata)-2),]
# Sex Pie Graph
png("Sex_piegraph.png", width=5, height=5, units="in", res=300)
metadata$sex <- as.character(metadata$sex)
metadata$sex[!metadata$sex %in% c("M","F")] <- "Unk"
pie(table(factor(metadata$sex)), labels=labels)
labels <- paste(names(table(factor(metadata$sex))), " (", table(factor(metadata$sex)), ")", sep="")
pie(table(factor(metadata$sex)), labels=labels, col=c("salmon", "cornflowerblue","grey50"), main="Sex")
dev.off()


par(mar=c(4,1,0,1))
# Age dotplot
png("Age_dotgraph.png", width=5, height=3, units="in", res=300)
stripchart(metadata$age, method = "stack", offset = .5, at = .15, pch = 19,main = "", xlab = "Donor Age (years)", frame.plot=FALSE, xaxt="n")
axis(1, at=seq(from=0, to=70, by=10))
abline(v=c(35, 60), col="grey35", lty=3)
legend("topleft", bty="n", c("Young"))
legend("top", bty="n", c("Middle"))
legend("topright", bty="n", c("Old"))
dev.off()

# BMI dotplot
png("BMI_dotgraph.png", width=5, height=4, units="in", res=300)
set.seed(3927)
stripchart(metadata$BMI, method="jitter", jitter=0.01, ylim=c(1,1.5), pch=19, main="", xlab="Donor BMI", frame.plot=FALSE, offset=.1, xaxt="n")
axis(1, at=seq(from=0, to=70, by=10))
abline(v=c(15, 18, 25, 30, 35, 40), col="grey35", lty=3)
text(21.5,1.1, "Normal", pos=3)
text(27.5,1.1, "Overweight", pos=3)
text(32.5,1.1, "Obese", pos=3)
dev.off()

