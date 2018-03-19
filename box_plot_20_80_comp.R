a <- read.table("DESeq2_EZH2KO_WT_enrich_padj005_protcod_norm_counts_protcod_logfold_above_below_05_plus1_ratios_EZH2_enriched.txt", sep = "\t", stringsAsFactor = FALSE, header = TRUE)
a1 <- a[, c(1,2)]
a1 <- data.frame(a1)
names(a1)[1] <- paste("Genes")
b <- read.table("DESeq2_EZH2KO_WT_enrich_padj005_protcod_WTVP_ind_EZH2_enriched.txt", sep = "\t", stringsAsFactor = FALSE, header = TRUE)
b1 <- b[, c(1,2)]
names(b1)[1] <- paste("Genes")
ab <- merge(a1, b1, by = "Genes", all.x = TRUE)
ab1 <- ab[is.na(ab$X.1.y),]
p <- read.table("Normalized_counts_protcod.txt", sep = "\t", stringsAsFactor = FALSE)
write.csv(p, "test.csv")
p <- read.csv("test.csv", stringsAsFactor = FALSE)
p1 <- p[, c(1,10,11)]
p1$WT_input1 <- p1$WT_input1 + 1
p1$WT_input2 <- p1$WT_input2 + 1
p1$avg <- (p1$WT_input1 + p1$WT_input2)/2
p2 <- p1$avg
p2 <- data.frame(p2)
p2 <- p1[, c(1,4)]
k1 <- ab1
nrow(k1)
ab2 <- merge(a1, b1, by = "Genes")
nrow(ab2)
k2 <- ab2
k1 <- k1$Genes
k1 <- data.frame(k1)
k2 <- k2$Genes
k2 <- data.frame(k2)
k2$type <- paste("Double")
k1$type <- paste("Single")
names(p2)[1] <- paste("Genes")
names(k1)[1] <- paste("Genes")
names(k2)[1] <- paste("Genes")
kp1 <- merge(k1, p2, by = "Genes")
kp2 <- merge(k2, p2, by = "Genes")
write.csv(kp1, "Single_silenced_genes.csv")
write.csv(kp2, "Double_silenced_genes.csv")
g <- read.table("gencode.v27.annotation.protcod.bed", sep = "\t", stringsAsFactor = FALSE)
g$V5 <- g$V3 - g$V2
g1 <- g[, c(4,5)]
names(g1)[1] <- paste("Genes")
names(g1)[2] <- paste("kb")
library(plyr)
g2 <- ddply(g1, .(Genes), summarize, kb=mean(kb))
a <- read.csv("Single_silenced_genes.csv", stringsAsFactor = FALSE, header = TRUE)
k <- read.table("Normalized_counts_VP55.txt", sep = "\t", stringsAsFactor = FALSE)
write.csv(k, "test.csv")
k <- read.csv("test.csv", stringsAsFactor = FALSE)
k$WT_VP_C11 <- k$WT_VP_C11 + 1
k$WT_VP_C12 <- k$WT_VP_C12 + 1
k$vp <- (k$WT_VP_C11 + k$WT_VP_C12)/2
k1 <- k[, c(1, 10)]
a1 <- a[, -c(1)]
names(k1)[1] <- paste("Genes")
ka <- merge(k1, a, by = "Genes")
kaa <- merge(ka, g2, by = "Genes")
kaa$avg <- kaa$avg/kaa$kb
kaa$vp <- kaa$vp/kaa$kb
ka1 <- kaa[, c(4,5)]
ka1$type <- paste("A")
ka2 <- kaa[, c(4,2)]
ka2$type <- paste("B")
b <- read.csv("Double_silenced_genes.csv", stringsAsFactor = FALSE, header = TRUE)
b1 <- b[, -c(1)]
kb <- merge(k1, b, by = "Genes")
kbb <- merge(kb, g2, by = "Genes")
kbb$avg <- kbb$avg/kbb$kb
kbb$vp <- kbb$vp/kbb$kb
kb1 <- kbb[, c(4,5)]
kb2 <- kbb[, c(4,2)]
kb1$type <- paste("C")
kb2$type <- paste("D")
names(ka2)[2] <- paste("avg")
names(kb2)[2] <- paste("avg")
kab <- rbind(ka1, ka2, kb1, kb2)
nrow(kab)
kab$log <- log(kab$avg, 2)
library(ggplot2)
write.csv(kab, "Box_plot_20_80_comp.csv")
pdf("Box_plot_20_80_comp.pdf")
p <- ggplot(kab, aes(kab$type, kab$log, fill = type))
p + geom_boxplot(outlier.colour = "white", notch = TRUE) + theme_bw() + theme(axis.text.x=element_text(angle=60, hjust=1), panel.grid.major = element_blank()) + stat_boxplot(geom ='errorbar', width = 0.2)
dev.off()


