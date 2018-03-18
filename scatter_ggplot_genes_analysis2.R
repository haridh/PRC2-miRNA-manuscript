library(plyr)
library(ggplot2)
a <- read.csv("DESeq2_WT_VP_vs_NO_padj005_protcod_norm_counts.csv", stringsAsFactor = FALSE)
a1 <- a[, c(2,4)]
a1 <- data.frame(a1)
b <- read.table("AGO2_iCLIP_BR1_counts.bed", sep = "\t", stringsAsFactor = FALSE)
b1 <- b[, c(13,25)]
b <- read.table("AGO2_iCLIP_BR2_counts.bed", sep = "\t", stringsAsFactor = FALSE)
b2 <- b[, c(13,25)]
b3 <- merge(b1, b2, by = "V13")
b3$counts <- (b3$V25.x + b3$V25.y)/2
b4 <- b3[, c(1,4)]
b4 <- unique(b4)
b5 <- ddply(b4, "V13", summarize, counts=sum(counts))
names(b5)[1] <- paste("Row.names")
ab <- merge(a1, b5, by = "Row.names")
ab <- subset(ab, ab$log2FoldChange >= 0.5)

ab$counts <- ab$counts + 1
k1 <- read.table("AGO2_common_peaks_noreps.lift.hg38.genes.bed", sep = "\t", stringsAsFactor = FALSE)
k1 <- k1[, c(4,5)]
k1 <- data.frame(k1)
names(k1)[1] <- paste("Row.names")
kab <- merge(ab, k1, by = "Row.names", all.x = TRUE)
kab <- unique(kab)
kab$label <- ifelse(!is.na(kab[,4]), "AGO2_Enriched", "Non-AGO2_enriched")
kab <- kab[, -c(4)]
k2 <- read.table("DESeq2_EZH2KO_WT_enrich_padj005_protcod_WTVP_ind_EZH2_enriched.txt", sep = "\t", stringsAsFactor = FALSE, header = TRUE)
k2 <- subset(k2, k2$WT_KO_inp1 < 1)
names(k2)[1] <- paste("Row.names")
k2 <- k2[, c(1,2)]
kkab <- merge(kab, k2, all.x = TRUE)
kkab$label <- ifelse(!is.na(kkab[,5]) & kkab$label == "AGO2_Enriched", "AGO2_Enriched_EZH2_repressed", ifelse(!is.na(kkab[,5]) & kkab$label == "Non-AGO2_enriched", "Non-AGO2_Enriched_EZH2_repressed", kkab$label))
m3 <- data.frame(kkab, row.names = 1)
m3$log_count <- log(m3$counts, 2)
m3$X.1 <- ifelse(m3$label == "AGO2_Enriched_EZH2_repressed", "D", "A")
m3$X.1 <- ifelse(m3$label == "Non-AGO2_Enriched_EZH2_repressed", "C", m3$X.1)
m3$X.1 <- ifelse(m3$label == "AGO2_Enriched", "B", m3$X.1)
library(ggplot2)
names(m3)[4] <- paste("tags")
m3 <- m3[with(m3, order(m3$tags)), ]
pdf("Scatter_clip_input_genes_log_change_based.pdf")
ggplot(m3, aes(x=log2FoldChange, y=log_count)) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_point(aes(shape=tags, colour = tags)) + scale_shape_manual(values = c("B" = 16,"C"= 1, "D"=16,  "A"= 1)) + scale_colour_manual(values = c("B" = "black","C"= "skyblue1", "D"="blue",  "A"= "black"))
dev.off()
write.csv(m3, "Scatter_clip_input_genes_log_change_based.csv")



a <- read.csv("DESeq2_WT_VP_vs_NO_padj005_protcod_norm_counts.csv", stringsAsFactor = FALSE)
a$avg <- (a$WT_NO_C11 + a$WT_NO_C12)/2
a$avg <- a$avg + 1
a$log_avg <- log(a$avg, 2)
a <- subset(a, a$log2FoldChange >= 0.5)
a1 <- a[, c(2,18)]

a1 <- data.frame(a1)
b <- read.table("AGO2_iCLIP_BR1_counts.bed", sep = "\t", stringsAsFactor = FALSE)
b1 <- b[, c(13,25)]
b <- read.table("AGO2_iCLIP_BR2_counts.bed", sep = "\t", stringsAsFactor = FALSE)
b2 <- b[, c(13,25)]
b3 <- merge(b1, b2, by = "V13")
b3$counts <- (b3$V25.x + b3$V25.y)/2
b4 <- b3[, c(1,4)]
b4 <- unique(b4)
b5 <- ddply(b4, "V13", summarize, counts=sum(counts))
names(b5)[1] <- paste("Row.names")
ab <- merge(a1, b5, by = "Row.names")
ab$counts <- ab$counts + 1
k1 <- read.table("AGO2_common_peaks_noreps.lift.hg38.genes.bed", sep = "\t", stringsAsFactor = FALSE)
k1 <- k1[, c(4,5)]
k1 <- data.frame(k1)
names(k1)[1] <- paste("Row.names")
kab <- merge(ab, k1, by = "Row.names", all.x = TRUE)
kab <- unique(kab)
kab$label <- ifelse(!is.na(kab[,4]), "AGO2_Enriched", "Non-AGO2_enriched")
kab <- kab[, -c(4)]
k2 <- read.table("DESeq2_EZH2KO_WT_enrich_padj005_protcod_WTVP_ind_EZH2_enriched.txt", sep = "\t", stringsAsFactor = FALSE, header     = TRUE)
k2 <- subset(k2, k2$WT_KO_inp1 < 1)
names(k2)[1] <- paste("Row.names")
k2 <- k2[, c(1,2)]
kkab <- merge(kab, k2, all.x = TRUE)
kkab$label <- ifelse(!is.na(kkab[,5]) & kkab$label == "AGO2_Enriched", "AGO2_Enriched_EZH2_repressed", ifelse(!is.na(kkab[,5]) & kkab$label == "Non-AGO2_enriched", "Non-AGO2_Enriched_EZH2_repressed", kkab$label))
m3 <- data.frame(kkab)
m3$log_count <- log(m3$counts, 2)
m3$X.1 <- ifelse(m3$label == "AGO2_Enriched_EZH2_repressed", "D", "A")
m3$X.1 <- ifelse(m3$label == "Non-AGO2_Enriched_EZH2_repressed", "C", m3$X.1)
m3$X.1 <- ifelse(m3$label == "AGO2_Enriched", "B", m3$X.1)
names(m3)[5] <- paste("tags")
m3 <- m3[with(m3, order(m3$tags)), ]
pdf("Scatter_clip_input_genes_exp_based.pdf")
ggplot(m3, aes(x=log_avg, y=log_count)) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_point(aes(colour = tags, shape = tags))+ scale_shape_manual(values = c("B" = 16,"C"= 1, "D"=16,  "A"= 1)) + scale_colour_manual(values = c("B" = "black","C"= "skyblue1", "D"="blue", "A" = "black"))
dev.off()
write.csv(m3,"Scatter_clip_input_genes_exp_based.csv")
