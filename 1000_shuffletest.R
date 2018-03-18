args <- commandArgs(trailingOnly = TRUE)
n=as.numeric(args[1])
File <- 1:n
Interactions <- 1:n
Genes_interacting <- 1:n
df <- data.frame(File, Interactions, Genes_interacting, stringsAsFactors = FALSE)
for(i in 1:n){
a <- read.csv("DIANA_miRNA_targets_pairs.csv", stringsAsFactor = FALSE)
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
}
write.csv(df, "Target_1000_shuffle_intersections_count.csv")
