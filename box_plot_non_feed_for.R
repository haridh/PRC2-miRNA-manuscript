a <- read.table("DESeq2_EZH2KO_WT_enrich_padj005_protcod_WTVP_ind_EZH2_enriched_all_ratios.txt", sep = "\t", stringsAsFactor = FALSE, header = TRUE)
a$avg <- (a$WT_KO_ind_1 + a$WT_KO_ind_1)/2
a1 <- subset(a, a$avg > 1)
a3 <- a1[, c(1)]
a3 <- data.frame(a3)
names(a3)[1] <- paste("Genes")
a3$temp <- 1:nrow(a3)
b <- read.table("DESeq2_EZH2KO_WT_enrich_padj005_protcod_WTVP_ind_EZH2_enriched.txt", sep = "\t", stringsAsFactor = FALSE, header = TRUE)
b <- subset(b, b$WT_KO_inp1 < 1)
b1 <- b[, c(1)]
b1 <- data.frame(b1)
names(b1)[1] <- paste("Genes")
print(nrow(a3))
print(nrow(b1))
ab <- merge(b1, a3, by = "Genes", all.x = TRUE)
ab1 <- ab[is.na(ab$temp),]
ab2 <- ab1[, c(1)]
ab2 <- data.frame(ab2, row.names = 1)
m <- read.table("Normalized_counts_protcod.txt", sep = "\t", stringsAsFactor = FALSE)
ma <- merge(ab2, m, by = "row.names")
k <- read.table("Normalized_counts_VP55.txt", sep = "\t", stringsAsFactor = FALSE)
names(ma)[1] <- paste("NAME")
write.csv(ma, "test.csv")
ma <- read.csv("test.csv", stringsAsFactor = FALSE)
ma <- ma[, -c(1)]
ma <- data.frame(ma, row.names = 1)
mak <- merge(ma, k, by = "row.names")
mak <- data.frame(mak, row.names = 1)
mak <- mak + 1
mk1 <- mak[, c(3,4,9,10,13,14,17,18)]
mk1$KO <- (mk1$KO_input1 + mk1$KO_input2)/2
mk1$WT <- (mk1$WT_input1 + mk1$WT_input2)/2
mk1$KO_VP <- (mk1$KO_VP_C1 + mk1$KO_VP_C3)/2
mk1$WT_VP <- (mk1$WT_VP_C11 + mk1$WT_VP_C12)/2
mk2 <- mk1[, c(9:12)]
p1 <- mk2[, c(1)]
p1 <- data.frame(p1)
p1$id <- paste("C")
p2 <- mk2[, c(2)]
p2 <- data.frame(p2)
p2$id <- paste("A")
p3 <- mk2[, c(3)]
p3 <- data.frame(p3)
p3$id <- paste("D")
p4 <- mk2[, c(4)]
p4 <- data.frame(p4)
p4$id <- paste("B")
names(p1)[1] <- paste("norm_reads")
names(p2)[1] <- paste("norm_reads")
names(p3)[1] <- paste("norm_reads")
names(p4)[1] <- paste("norm_reads")
p <- rbind(p1, p2, p3, p4)
p$log_norm_exp <- log(p$norm_reads, 2)
library(ggplot2)
s <- p
write.csv(s, "Box_plot_EZH2_miR_repressed_genes_non_feed.csv")
pdf("Box_plot_EZH2_miR_repressed_genes_non_feed.pdf")
p <- ggplot(s, aes(s$id, s$log_norm_exp, fill = id))
p + geom_boxplot(outlier.colour = "white", notch = TRUE) + theme_bw() + theme(axis.text.x=element_text(angle=60, hjust=1), panel.grid.major = element_blank()) + stat_boxplot(geom ='errorbar', width = 0.2) + coord_cartesian(ylim = c(0,20))
dev.off()
