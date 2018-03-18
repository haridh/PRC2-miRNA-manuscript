a <- read.csv("DESeq2_WT_VP_vs_NO_padj005_protcod_norm_counts.csv", stringsAsFactor = FALSE, header = TRUE)
a1 <- subset(a, a$log2FoldChange >=0.5)
a1 <- a1$Row.names
a1 <- data.frame(a1)
names(a1)[1] <- paste("genes")
b <- read.csv("EZH2_peak_genes_order.csv", stringsAsFactor = FALSE)
b1 <- b$NAME
b1 <- data.frame(b1)
names(b1)[1] <- paste("genes")
ab <- merge(a1, b1, by = "genes")
ab <- unique(ab)
c <- read.table("DESeq2_EZH2KO_WT_enrich_padj005_protcod_norm_counts_protcod_logfold_above_below_05_plus1_ratios.txt", header = TRUE, sep = "\t", stringsAsFactor = FALSE)
c1 <- c[, c(1:4)]
c1 <- data.frame(c1)
names(c1)[1] <- paste("genes")
abc <- merge(ab, c1, by = "genes")
abc <- unique(abc)
abc$avg <- (abc$WT_KO_inp1 + abc$WT_KO_inp2)/2
abc <- abc[with(abc, order(avg)), ]
abc1 <- abc[, -c(5)]
write.table(abc1, "DESeq2_EZH2KO_WT_enrich_padj005_protcod_WTVP_ind_EZH2_enriched.txt", sep = "\t", row.names = FALSE, quote = FALSE)
k <- subset(abc, abc$avg < 1)
k1 <- k$genes
k1 <- data.frame(k1)
names(k1)[1] <- paste("Row.names")
p <- read.table("DESeq2_KO_WT_WTVP_NO_intersect_protcod_all_ratios.txt", sep = "\t", stringsAsFactor = FALSE, header = TRUE)
pk <- merge(p, k1, by = "Row.names")
pk$avg1 <- log(pk$WT_KO_ind_1,2)
pk$avg2 <- log(pk$WT_KO_ind_2,2)
pk$avg <- (pk$avg1 + pk$avg2)/2
pk$test <- pk$avg1/pk$avg2
pk <- subset(pk, pk$test > 0)
h <- subset(pk, (pk$avg >= 0.5))
print(nrow(h))
pkm <- subset(pk, (pk$avg >= 0.5)| (pk$avg <= -0.5))
pkm <- pkm[with(pkm, order(-avg)), ]
pkm1 <- pkm[, -c(11:14)]
write.table(pkm1, "DESeq2_EZH2KO_WT_enrich_padj005_protcod_WTVP_ind_EZH2_enriched_all_ratios.txt", sep = "\t", row.names = FALSE, quote = FALSE)


