a <- read.csv("DESeq2_EZH2KO_WT_enrich_padj005_protcod_norm_counts.csv", stringsAsFactor = FALSE)
b <- read.csv("DESeq2_WT_VP_vs_NO_padj005_protcod_norm_counts.csv", stringsAsFactor = FALSE)
a <- a[, -c(1)]
b <- b[, -c(1)]
ab <- merge(a, b, by = "Row.names")
ab <- unique(ab)
write.csv(ab, "DESeq2_KO_WT_WTVP_NO_VP_protcod_intersect.csv")
m <- read.csv("DESeq2_KO_VP_vs_NO_padj005_protcod_norm_counts.csv", stringsAsFactor = FALSE)
m <- m[, -c(1)]
ma <- merge(a, m, by = "Row.names")
ma <- unique(ma)
write.csv(ma, "DESeq2_KO_WT_KOVP_NO_VP_protcod_intersect.csv")
k <- read.table("Normalized_counts_VP55.txt", sep = "\t", stringsAsFactor = FALSE)
k <- k+1
k$KO_ind_1 <- k$KO_VP_C1/k$KO_NO_C1
k$KO_ind_2 <- k$KO_VP_C3/k$KO_NO_C3
k$WT_ind_1 <- k$WT_VP_C11/k$WT_NO_C11
k$WT_ind_2 <- k$WT_VP_C12/k$WT_NO_C12
k1 <- k[, c(9:12)]
a1 <- subset(a, a$log2FoldChange <= -0.5 | a$log2FoldChange >= 0.5)
a2 <- a1[, c(1,10,11,16,17)]
a2 <- data.frame(a2, row.names = 1)
a2 <- a2 +1
a2$WT_KO1 <- a2$WT_input1/a2$KO_input1
a2$WT_KO2 <- a2$WT_input2/a2$KO_input2
a3 <- a2[, c(5,6)]
ak <- merge(a3, k1, by = "row.names")
ak <- unique(ak)
write.csv(ak, "DESeq2_EZH2_reg_KO_WT_VP_KO_WT_protcod_norm_cts_ratios.csv")
j <- read.csv("DESeq2_EZH2_reg_KO_WT_VP_KO_WT_protcod_norm_cts_ratios.csv", stringsAsFactor = FALSE, header = TRUE)
j$WT_KO_ind_1 <- j$WT_ind_1/j$KO_ind_1
j$WT_KO_ind_2 <- j$WT_ind_2/j$KO_ind_2
j2 <- j[, c(2,2,3,4,7,8,5,6,9,10)]
write.table(j2, "DESeq2_EZH2_reg_KO_WT_VP_KO_WT_protcod_norm_cts_ALL_ratios.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
names(j2)[1] <- paste("Row.names")
vk <- read.csv("DESeq2_KO_VP_vs_NO_padj005_protcod_norm_counts.csv", stringsAsFactor = FALSE)
vk0 <- subset(vk, vk$log2FoldChange <= -0.5 | vk$log2FoldChange >= 0.5)
vk0 <- vk0[, -c(1)]
vk0 <- vk0[, c(1,8,9,10,11)]
vk0 <- data.frame(vk0, row.names = 1)
vk0 <- vk0+1
vk0$KO_ind_1 <- vk0$KO_VP_C1/vk0$KO_NO_C1
vk0$KO_ind_2 <- vk0$KO_VP_C3/vk0$KO_NO_C3
vk1 <- vk0[, c(5,6)]
h <- a1[, c(1,10,11, 16,17)]
h <- data.frame(h, row.names = 1)
h <- h+1
h$WT_KO1 <- h$WT_input1/h$KO_input1
h$WT_KO2 <- h$WT_input2/h$KO_input2
h1 <- h[, c(5,6)]
hvk <- merge(h1, vk1, by = "row.names")
hvk <- hvk[, c(1,1,2,3,4,5)]
hvk$avg1 <- (hvk$WT_KO1 + hvk$WT_KO2)/2
hvk$avg2 <- (hvk$KO_ind_1 + hvk$KO_ind_2)/2
hvk1 <- subset(hvk, hvk$avg1 < 1)
hvk2 <- subset(hvk, hvk$avg1 >=1)
hvk1 <- hvk1[with(hvk1, order(-avg2)), ]
hvk2 <- hvk2[with(hvk2, order(-avg2)), ]
hvk3 <- rbind(hvk1, hvk2)
hvk3 <- hvk3[, -c(7:8)]
write.table(hvk3, "DESeq2_KO_WT_KOVP_NO_intersect_protcod.txt", sep = "\t", row.names = FALSE, quote = FALSE)
vk <- read.csv("DESeq2_WT_VP_vs_NO_padj005_protcod_norm_counts.csv", stringsAsFactor = FALSE)
vk0 <- subset(vk, vk$log2FoldChange <= -0.5 | vk$log2FoldChange >= 0.5)
vk0 <- vk0[, -c(1)]
vk0 <- vk0[, c(1,12,13,14,15)]
vk0 <- data.frame(vk0, row.names = 1)
vk0 <- vk0+1
vk0$WT_ind_1 <- vk0$WT_VP_C11/vk0$WT_NO_C11
vk0$WT_ind_2 <- vk0$WT_VP_C12/vk0$WT_NO_C12
vk1 <- vk0[, c(5,6)]
h <- a1[, c(1,10,11, 16,17)]
h <- data.frame(h, row.names = 1)
h <- h+1
h$WT_KO1 <- h$WT_input1/h$KO_input1
h$WT_KO2 <- h$WT_input2/h$KO_input2
h1 <- h[, c(5,6)]
hvk <- merge(h1, vk1, by = "row.names")
hvk <- hvk[, c(1,1,2,3,4,5)]
names(hvk)[1] <- paste("Row.names")
hvkj <- merge(hvk, j2, by = "Row.names")
fin <- hvkj[, c(1,1,3,9,10,11,12,13,14,15)]
fin <- data.frame(fin)
fin$avg1 <- (fin$WT_KO1.x + fin$WT_KO2.y)/2
fin$avg2 <- (fin$WT_KO_ind_1 + fin$WT_KO_ind_2)/2
fin1 <- subset(fin, fin$avg1 < 1)
fin1$avg3 <- (fin1$WT_ind_1.y + fin1$WT_ind_2.y)/2
fin11 <- subset(fin1, fin1$avg3 >=1)
fin12 <- subset(fin1, fin1$avg3 <1)
fin11 <- fin11[with(fin11, order(-avg2)), ]
fin12 <- fin12[with(fin12, order(-avg2)), ]
fin2 <- rbind(fin11, fin12)
fin2 <- fin2[, -c(11:13)]


hvk$avg1 <- (hvk$WT_KO1 + hvk$WT_KO2)/2
hvk$avg2 <- (hvk$WT_ind_1 + hvk$WT_ind_2)/2
hvk1 <- subset(hvk, hvk$avg1 < 1)
hvk2 <- subset(hvk, hvk$avg1 >=1)
hvk1 <- hvk1[with(hvk1, order(-avg2)), ]
hvk2 <- hvk2[with(hvk2, order(-avg2)), ]
hvk3 <- rbind(hvk1, hvk2)
hvk3 <- hvk3[, -c(7:8)]
write.table(fin2, "DESeq2_KO_WT_WTVP_NO_intersect_protcod_all_ratios.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(hvk3, "DESeq2_KO_WT_WTVP_NO_intersect_protcod.txt", sep = "\t", row.names = FALSE, quote = FALSE)
