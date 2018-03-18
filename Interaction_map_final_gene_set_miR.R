s1 <- read.table("Final_gene_set.bed", sep = "\t", stringsAsFactor = FALSE)
s1 <- s1[, c(4)]
s1 <- data.frame(s1)
names(s1)[1] <- paste("Target")
m1 <- read.csv("EZH2_activated_miRs_DIANA_targets.csv", stringsAsFactor = FALSE)
m1 <- m1[, -c(1)]
names(m1)[1] <- paste("miR")
names(m1)[2] <- paste("Target")
m1$count <- 1
sm <- merge(s1, m1, by = "Target", all.x = TRUE)
sm$miR <- ifelse(is.na(sm$miR), 'No_miR', as.character(sm$miR))
sm$count[is.na(sm$count)] <- 0
p1 <- reshape(sm, idvar = "miR", timevar = "Target", direction = "wide")
p1[is.na(p1)] <- 0
library(RColorBrewer)
library(gplots)
hmcol = colorRampPalette(brewer.pal(9, "Reds"))(100)
pdf("EZH2_activated_miRs_DIANA_targets_final_gene_set.pdf")
write.table(p1, "EZH2_activated_miRs_DIANA_targets_final_gene_set_reshaped.txt", sep = "\t", row.names = FALSE, quote = FALSE)
a <- read.table("EZH2_activated_miRs_DIANA_targets_final_gene_set_reshaped.txt", sep = "\t", stringsAsFactor = FALSE, header = TRUE)
a1 <- data.frame(a, row.names = 1)
a1 <- a1[setdiff(rownames(a1),"No_miR"),]
zzz <- as.matrix(a1)
heatmap.2(zzz, trace="none", col = hmcol, dendrogram= "none", margin=c(10, 10))
dev.off()

