file.names <- dir("./", pattern = ".txt")
for(i in 1:length(file.names)){
p <- read.table(file.names[i],header=TRUE,sep = "\t", stringsAsFactor = FALSE)
q <- p[, c(1,7)]
if(i==1){
out <- q} else{
out <- merge(out, q, by = "Geneid")
}}
file.names
names(out)[2] <- paste("KO_NO_C1")
names(out)[3] <- paste("KO_NO_C3")
names(out)[4] <- paste("KO_VP_C1")
names(out)[5] <- paste("KO_VP_C3")
names(out)[6] <- paste("WT_NO_C11")
names(out)[7] <- paste("WT_NO_C12")
names(out)[8] <- paste("WT_VP_C11")
names(out)[9] <- paste("WT_VP_C12")
out <- unique(out)
gen <- read.table("gencode.v27.annotation.protcod.bed", sep = "\t", stringsAsFactor = FALSE)
gen1 <- gen$V4
gen1 <- unique(gen1)
gen1 <- data.frame(gen1)
names(gen1)[1] <- paste("Geneid")
out1 <- merge(out, gen1, by = "Geneid")
write.table(out1, "VP_protcod_Raw_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)
cts <- data.frame(out1, row.names = 1)
cts <- round(cts)
suppressMessages(library(DESeq2))
files = c("KO_NO_C1", "KO_NO_C3", "KO_VP_C1", "KO_VP_C3", "WT_NO_C11", "WT_NO_C12", "WT_VP_C11", "WT_VP_C12")
condition = c("KO_NO_C", "KO_NO_C", "KO_VP_C", "KO_VP_C", "WT_NO_C1", "WT_NO_C1", "WT_VP_C1", "WT_VP_C1")
sampleTable <- data.frame(sampleName = files, condition = condition)
dds <- DESeqDataSetFromMatrix(countData = cts, colData = sampleTable, design = ~ condition)
dds <- DESeq(dds)
norm <- counts(dds, normalized = TRUE)
write.table(norm, "Normalized_counts_VP55.txt", sep = "\t", quote = FALSE)

ka <- norm[, c(1:4)]
ka <- round(ka)
file1 <- c("KO_NO_C1", "KO_NO_C3", "KO_VP_C1", "KO_VP_C3")
condition1 <- c("KO_NO_C", "KO_NO_C", "KO_VP_C", "KO_VP_C")
sampleTable <- data.frame(sampleName = file1, condition = condition1)
dds <- DESeqDataSetFromMatrix(countData = ka, colData = sampleTable, design = ~ condition)
sizeFactors(dds) = 1
dds$condition <- relevel(dds$condition, "KO_NO_C")
dds <- DESeq(dds)
res <- results(dds)
write.csv(res, "DESeq2_KO_VP_vs_NO.csv")
k <- read.csv("DESeq2_KO_VP_vs_NO.csv", stringsAsFactor = FALSE, row.names = 1)
k1 <- subset(k, k$padj < 0.05)
kn <- merge(k1, norm, by = "row.names")
kn <- unique(kn)
write.csv(kn, "DESeq2_KO_VP_vs_NO_padj005_protcod_norm_counts.csv")

wa <- norm[, c(5:8)]
wa <- round(wa)
file2 <- c("WT_NO_C11", "WT_NO_C12", "WT_VP_C11", "WT_VP_C12")
condition2 <- c("WT_NO_C1", "WT_NO_C1", "WT_VP_C1", "WT_VP_C1")
sampleTable <- data.frame(sampleName = file2, condition = condition2)
dds <- DESeqDataSetFromMatrix(countData = wa, colData = sampleTable, design = ~ condition)
sizeFactors(dds) = 1
dds$condition <- relevel(dds$condition, "WT_NO_C1")
dds <- DESeq(dds)
res <- results(dds)
write.csv(res, "DESeq2_WT_VP_vs_NO.csv")
k <- read.csv("DESeq2_WT_VP_vs_NO.csv", stringsAsFactor = FALSE, row.names = 1)
k1 <- subset(k, k$padj < 0.05)
kn <- merge(k1, norm, by = "row.names")
kn <- unique(kn)
write.csv(kn, "DESeq2_WT_VP_vs_NO_padj005_protcod_norm_counts.csv")


wk <- norm[, c(3,4,7,8)]
wk <- round(wk)
file3 <- c("KO_VP_C1", "KO_VP_C3", "WT_VP_C11", "WT_VP_C12")
condition3 <- c("KO_VP_C", "KO_VP_C", "WT_VP_C1", "WT_VP_C1")
sampleTable <- data.frame(sampleName = file3, condition = condition3)
dds <- DESeqDataSetFromMatrix(countData = wk, colData = sampleTable, design = ~ condition)
sizeFactors(dds) = 1
dds$condition <- relevel(dds$condition, "WT_VP_C1")
dds <- DESeq(dds)
res <- results(dds)
write.csv(res, "DESeq2_KOVP_vs_WTVP.csv")
k <- read.csv("DESeq2_KOVP_vs_WTVP.csv", stringsAsFactor = FALSE, row.names = 1)
k1 <- subset(k, k$padj < 0.05)
kn <- merge(k1, norm, by = "row.names")
kn <- unique(kn)
write.csv(kn, "DESeq2_KOVP_vs_WTVP_padj005_protcod_norm_counts.csv")
