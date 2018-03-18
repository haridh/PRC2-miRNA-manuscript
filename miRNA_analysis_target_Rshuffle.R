args <- commandArgs(trailingOnly = TRUE)
n=as.numeric(args[1])

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
zzz <- as.matrix(a1)
heatmap.2(zzz, trace="none", col = hmcol, dendrogram= "none", margin=c(10, 10))
dev.off()

test <- merge(s1, m1, by = "Target")
tt <- test$Target
tt <- data.frame(tt)
tt <- unique(tt)
print(nrow(tt))

File <- 1:n
Interactions <- 1:n
Genes_interacting <- 1:n
df <- data.frame(File, Interactions, Genes_interacting, stringsAsFactors = FALSE)
for(i in 1:n){
a <- read.csv("DIANA_miRNA_targets_pairs.csv", stringsAsFactor = FALSE)
a1 <- a$miR
a1 <- data.frame(a1)
a1 <- unique(a1)
b <- read.table("gencode.v27.annotation.protcod.bed", sep = "\t", stringsAsFactor = FALSE)
b1 <- b$V4
b1 <- data.frame(b1)
b1 <- unique(b1)
b2 <- b1[sample(nrow(b1)),]
b2 <- data.frame(b2)
names(b2)[1] <- paste("Target")
s1 <- read.table("Final_gene_set.bed", sep = "\t", stringsAsFactor = FALSE)
s1 <- s1[, c(4)]
s1 <- data.frame(s1)
names(s1)[1] <- paste("Target")
library(dplyr)
bs <- setdiff(b2, s1)
b3 <- bs[c(1:116),]
b3 <- data.frame(b3)
names(b3)[1] <- paste("Target")
aa <- merge(b3, a, by = "Target")
aa <- unique(aa)
k1 <- data.frame(aa)
k1 <- k1[, -c(2)]
m1 <- read.csv("EZH2_activated_miRs_DIANA_targets.csv", stringsAsFactor = FALSE)
m1 <- m1[, -c(1)]
names(m1)[1] <- paste("miR")
names(m1)[2] <- paste("Target")
m1 <- m1$miR
m1 <- unique(m1)
m1 <- data.frame(m1)
names(m1)[1] <- paste("miR")
mk <- merge(m1, k1, by = "miR")
mk <- unique(mk)
ct <- mk$Target
ct <- data.frame(ct)
ct <- unique(ct)
Interactions <- nrow(ct)
df$File[i] <- i
df$Interactions[i] <- nrow(mk)
df$Genes_interacting[i] <- nrow(ct)
mk$count <- 1
b4 <- merge(b3, mk, by = "Target", all.x = TRUE)
b4$miR <- ifelse(is.na(b4$miR), 'No_miR', as.character(b4$miR))
b4$count[is.na(b4$count)] <- 0
p1 <- reshape(b4, idvar = "miR", timevar = "Target", direction = "wide")
p1[is.na(p1)] <- 0
print(ncol(p1))
library(RColorBrewer)
library(gplots)
hmcol = colorRampPalette(brewer.pal(9, "Reds"))(100)
file <- paste("Target_shuffle", i, ".pdf", sep = "")
pdf(file)
write.table(p1, "test.txt", sep = "\t", row.names = FALSE, quote = FALSE)
a <- read.table("test.txt", sep = "\t", stringsAsFactor = FALSE, header = TRUE)
a1 <- data.frame(a, row.names = 1)
zzz <- as.matrix(a1)
heatmap.2(zzz, trace="none", col = hmcol, dendrogram= "none", margin=c(10, 10))
dev.off()
}
write.csv(df, "Target_shuffle_intersections_count.csv")


