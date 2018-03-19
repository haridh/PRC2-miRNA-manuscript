library("ggplot2")
library("plyr")
library("dplyr")
library("zoo")
a <- read.table("T98_EZH2_HsAbcam1.b50.score.width.tss.pergene.txt", sep = "\t", stringsAsFactor = FALSE, header = TRUE)
a$sum <- rowSums(a[,2:201])
a1 <- a[, c(1,202)]
s1 <- sum(a1$sum)
a1$sum <- a1$sum/s1
b <- read.table("T98_H3k27me3_HsCmb.b50.score.width.tss.pergene.txt", sep = "\t", stringsAsFactor = FALSE, header = TRUE)
b$sum <- rowSums(b[,2:201])
b1 <- b[, c(1,202)]
s2 <- sum(b1$sum)
b1$sum <- b1$sum/s2
c <- read.table("t98_h3k4me3_b1.b50.score.width.tss.pergene.txt", sep = "\t", stringsAsFactor = FALSE, header = TRUE)
c$sum <- rowSums(c[,2:201])
c1 <- c[, c(1,202)]
s3 <- sum(c1$sum)
c1$sum <- c1$sum/s3
names(a1)[2] <- paste("EZH2")
names(b1)[2] <- paste("H3K27me3")
names(c1)[2] <- paste("H3K4me3")
ab <- merge(a1, b1, by = "NAME")
abc <- merge(ab, c1, by = "NAME")
m <- read.table("../DESeq2_EZH2KO_WT_enrich_padj005_protcod_norm_counts_protcod_logfold_above_below_05_plus1_ratios.txt", header = TRUE, sep = "\t", stringsAsFactor = FALSE)
m1 <- m[, c(1,2)]
names(m1)[1] <- paste("NAME")
m1$X.1 <- 1:nrow(m1)
names(m1)[2] <- paste("V5")
m1$NAME <- paste(m1$V5, m1$NAME, sep = "_")
mabc <- merge(m1, abc, by = "NAME")
write.csv(mabc, "EZH2_reg_K27_k4_ez_vertical.csv")
m2 <- ddply(mabc, c("NAME"),summarise, V5=mean(V5), EZH2=mean(EZH2), K27=mean(H3K27me3), k4m3=mean(H3K4me3))
m3 <- m2[with(m2, order(V5)), ]
k1 <- m3[, c(2,3)]
k2 <- m3[, c(2,4)]
k3 <- m3[, c(2,5)]
k1$gp <- paste("EZH2")
k2$gp <- paste("H3K27me3")
k3$gp <- paste("H3K4me3")
names(k1)[2] <- paste("Peak")
names(k2)[2] <- paste("Peak")
names(k3)[2] <- paste("Peak")
k1$roll <- rollmean(k1$Peak, 50,fill = list(NA, NULL, NA))
k2$roll <- rollmean(k2$Peak, 50,fill = list(NA, NULL, NA))
k3$roll <- rollmean(k3$Peak, 50,fill = list(NA, NULL, NA))
kkk <- rbind(k1, k2, k3)
kkk[is.na(kkk)] <- 0
kkk$log <- log(kkk$roll, 2)
pdf("EZH2_reg_K27_k4_ez_vertical.pdf")
p<-ggplot(data=kkk, aes(x=V5, y=log, colour = gp)) + geom_line() + scale_color_manual(values=c("blue", "red", "green"))
p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
dev.off()


