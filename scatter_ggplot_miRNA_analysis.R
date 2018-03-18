library(ggplot2)
library(plyr)
library(ggrepel)
a <- read.table("AGO2_iCLIP_BR1_miR_counts.bed", sep = "\t", stringsAsFactor = FALSE)
a1 <- a[, c(4,7)]
a <- read.table("AGO2_iCLIP_BR2_miR_counts.bed", sep = "\t", stringsAsFactor = FALSE)
b1 <- a[, c(4,7)]
k <- read.table("miRNA_input_counts.bed", sep = "\t", stringsAsFactor = FALSE)
k1 <- k[, c(4,7)]
ab <- merge(a1, b1, by = "V4")
ab <- subset(ab, ab$V7.x > 0)
ab <- subset(ab, ab$V7.y > 0)
ab$clip_counts <- (ab$V7.x + ab$V7.y)/2
ab <- ab[, c(1,4)]
abk <- merge(ab, k1, by = "V4")
names(abk)[3] <- paste("Input")
abk1 <- subset(abk, abk$Input > 0)
abk1 <- unique(abk1)
b5 <- ddply(abk1, "V4", summarize, clip_counts=sum(clip_counts), Input = sum(Input))
b5$log_clip <- log(b5$clip_counts, 2)
b5$log_input <- log(b5$Input, 2)
library(splitstackshape)
p <- cSplit(b5, "V4", sep = "-")
p <- data.frame(p)
p$V4 <- paste(p$V4_2, p$V4_3, p$V4_4, sep = "-")
b5 <- p
pdf("scatter_miRNA_clip_input.pdf")
ggplot(b5, aes(x=log_clip, y=log_input)) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_point(color = "black", size = 1) + geom_label_repel(aes(label=ifelse(b5$V4 == "miR-27a-3p" | b5$V4 == "miR-21-5p" | b5$V4 == "miR-29a-3p" | b5$V4 == "miR-125b-5p" | b5$V4 == "miR-23a-3p" | b5$V4 == "miR-27b-3p" | b5$V4 == "miR-29b-3p" | b5$V4 == "miR-21-3p" | b5$V4 == "miR-193b-3p" | b5$V4 == "miR-222-3p" | b5$V4 == "miR-23b-3p" | b5$V4 == "let-7a-5p" | b5$V4 == "miR-29c-3p" | b5$V4 == "let-7b-5p" | b5$V4 == "miR-16-5p" | b5$V4 == "miR-193a-3p" | b5$V4 == "miR-484-NA" | b5$V4 == "miR-224-5p" | b5$V4 == "miR-22-3p" | b5$V4 == "miR-125a-5p", as.character(b5$V4),'')), size = 2, box.padding = 0.35, point.padding = 0.5, segment.color = 'blue')
dev.off()
