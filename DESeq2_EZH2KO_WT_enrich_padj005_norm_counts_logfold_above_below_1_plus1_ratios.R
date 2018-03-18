a <- read.csv("DESeq2_EZH2KO_WT_enrich_padj005_protcod_norm_counts.csv", stringsAsFactor = FALSE)
a <- a[, -c(1)]
a <- data.frame(a, row.names = 1)
b <- subset(a, (a$log2FoldChange <= -0.5 | a$log2FoldChange >= 0.5))
b <- b[,-c(1:6)]
b <- b + 1
write.csv(b, "DESeq2_EZH2KO_WT_enrich_padj005_protcod_norm_counts_protcod_logfold_above_below_05_plus1.csv")
b <- read.csv("DESeq2_EZH2KO_WT_enrich_padj005_protcod_norm_counts_protcod_logfold_above_below_05_plus1.csv", stringsAsFactor = FALSE)
b <- data.frame(b)
b <- unique(b)
b$WT_KO_inp1 <- b$WT_input1/b$KO_input1
b$WT_KO_inp2 <- b$WT_input2/b$KO_input2
b$WT_SUZKO1 <- b$WT_input1/b$SUZKO_input1
b$WT_SUZKO2 <- b$WT_input2/b$SUZKO_input2
c <- b[, c(1,1,12,13,14,15)]
c$avg <- (c$WT_KO_inp1+c$WT_KO_inp2)/2
c <- c[with(c, order(avg)), ]
c1 <- c[, -c(7)]
write.table(c1, "DESeq2_EZH2KO_WT_enrich_padj005_protcod_norm_counts_protcod_logfold_above_below_05_plus1_ratios.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
