file.names <- dir("./", pattern = ".txt")
for(i in 1:length(file.names)){
p <- read.table(file.names[i],header=TRUE,sep = "\t", stringsAsFactor = FALSE)
q <- p[, c(1,7)]
if(i==1){
out <- q} else{
out <- merge(out, q, by = "Geneid")
}}
file.names
names(out)[2] <- paste("KO_AGO1")
names(out)[3] <- paste("KO_AGO2")
names(out)[4] <- paste("KO_input1")
names(out)[5] <- paste("KO_input2")
names(out)[6] <- paste("SUZKO_input1")
names(out)[7] <- paste("SUZKO_input2")
names(out)[8] <- paste("WT_AGO1")
names(out)[9] <- paste("WT_AGO2")
names(out)[10] <- paste("WT_input1")
names(out)[11] <- paste("WT_input2")
out <- unique(out)
gen <- read.table("gencode.v27.annotation.protcod.bed", sep = "\t", stringsAsFactor = FALSE)
gen1 <- gen$V4
gen1 <- unique(gen1)
gen1 <- data.frame(gen1)
names(gen1)[1] <- paste("Geneid")
out1 <- merge(out, gen1, by = "Geneid")
write.table(out1, "Raw_counts_protcod.txt", sep = "\t", row.names = FALSE, quote = FALSE)
cts <- data.frame(out1, row.names = 1)
cts <- round(cts)
suppressMessages(library(DESeq2))
files = c("KO_AGO1", "KO_AGO2", "KO_input1", "KO_input2", "SUZKO_input1", "SUZKO_input2", "WT_AGO1", "WT_AGO2", "WT_input1", "WT_input2")
condition = c("KO_AGO", "KO_AGO", "KO_input", "KO_input", "SUZKO_input", "SUZKO_input", "WT_AGO", "WT_AGO", "WT_input", "WT_input")
sampleTable <- data.frame(sampleName = files, condition = condition)
dds <- DESeqDataSetFromMatrix(countData = cts, colData = sampleTable, design = ~ condition)
dds <- DESeq(dds)
norm <- counts(dds, normalized = TRUE)
write.table(norm, "Normalized_counts_protcod.txt", sep = "\t", quote = FALSE)

ka <- norm[, c(1:4)]
ka <- round(ka)
file1 <- c("KO_AGO1", "KO_AGO2", "KO_input1", "KO_input2")
condition1 <- c("KO_AGO", "KO_AGO", "KO_input", "KO_input")
sampleTable <- data.frame(sampleName = file1, condition = condition1)
dds <- DESeqDataSetFromMatrix(countData = ka, colData = sampleTable, design = ~ condition)
sizeFactors(dds) = 1
dds$condition <- relevel(dds$condition, "KO_input")
dds <- DESeq(dds)
res <- results(dds)
write.csv(res, "DESeq2_AGO2_KO_enrich.csv")
k <- read.csv("DESeq2_AGO2_KO_enrich.csv", stringsAsFactor = FALSE, row.names = 1)
k1 <- subset(k, k$padj < 0.05)
kn <- merge(k1, norm, by = "row.names")
kn <- unique(kn)
write.csv(kn, "DESeq2_AGO2_KO_enrich_padj005_protcod_norm_counts.csv")


wa <- norm[, c(7:10)]
wa <- round(wa)
file2 <- c("WT_AGO1", "WT_AGO2", "WT_input1", "WT_input2")
condition2 <- c("WT_AGO", "WT_AGO", "WT_input", "WT_input")
sampleTable <- data.frame(sampleName = file2, condition = condition2)
dds <- DESeqDataSetFromMatrix(countData = wa, colData = sampleTable, design = ~ condition)
sizeFactors(dds) = 1
dds$condition <- relevel(dds$condition, "WT_input")
dds <- DESeq(dds)
res <- results(dds)
write.csv(res, "DESeq2_AGO2_WT_enrich.csv")
k <- read.csv("DESeq2_AGO2_WT_enrich.csv", stringsAsFactor = FALSE, row.names = 1)
k1 <- subset(k, k$padj < 0.05)
kn <- merge(k1, norm, by = "row.names")
kn <- unique(kn)
write.csv(kn, "DESeq2_AGO2_WT_enrich_padj005_protcod_norm_counts.csv")

su <- norm[, c(5,6,9,10)]
su <- round(su)
file3 <- c("SUZKO_input1", "SUZKO_input2", "WT_input1", "WT_input2")
condition3 <- c("SUZKO_input", "SUZKO_input", "WT_input", "WT_input")
sampleTable <- data.frame(sampleName = file3, condition = condition3)
dds <- DESeqDataSetFromMatrix(countData = su, colData = sampleTable, design = ~ condition)
sizeFactors(dds) = 1
dds$condition <- relevel(dds$condition, "WT_input")
dds <- DESeq(dds)
res <- results(dds)
write.csv(res, "DESeq2_SUZKO_WT_enrich.csv")
k <- read.csv("DESeq2_SUZKO_WT_enrich.csv", stringsAsFactor = FALSE, row.names = 1)
k1 <- subset(k, k$padj < 0.05)
kn <- merge(k1, norm, by = "row.names")
kn <- unique(kn)
write.csv(kn, "DESeq2_SUZKO_WT_enrich_padj005_protcod_norm_counts.csv")


ez <- norm[, c(3,4,9,10)]
ez <- round(ez)
file4 <- c("KO_input1", "KO_input2", "WT_input1", "WT_input2")
condition4 <- c("KO_input", "KO_input", "WT_input", "WT_input")
sampleTable <- data.frame(sampleName = file4, condition = condition4)
dds <- DESeqDataSetFromMatrix(countData = ez, colData = sampleTable, design = ~ condition)
sizeFactors(dds) = 1
dds$condition <- relevel(dds$condition, "WT_input")
dds <- DESeq(dds)
res <- results(dds)
write.csv(res, "DESeq2_EZH2KO_WT_enrich.csv")
k <- read.csv("DESeq2_EZH2KO_WT_enrich.csv", stringsAsFactor = FALSE, row.names = 1)
k1 <- subset(k, k$padj < 0.05)
kn <- merge(k1, norm, by = "row.names")
kn <- unique(kn)
write.csv(kn, "DESeq2_EZH2KO_WT_enrich_padj005_protcod_norm_counts.csv")


