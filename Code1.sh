## Script for doing differential expression analysis using DESeq2 package in R
## Author: Dr. Rahul Kumar

rm(list=ls())
options(scipen=999)
library(DESeq2)
library(ggplot2)
library(reshape2)
library(ggfortify)
library(ComplexHeatmap)

## Data processing
data <- read.table("Counts.txt", sep="\t", header=TRUE, check.names=FALSE)
row.names(data) <- make.names(data$Gene, unique = TRUE, allow_ = FALSE)
tpm <- read.table("TPM.txt", sep="\t", header=TRUE, check.names=FALSE)

## Add median to the count matrix 
data$Untreated_median <- apply(data[,c(2:4)], 1, median)
data$Treated_median <- apply(data[,c(5:7)], 1, median)

## DESeq2 ##
## Design matrix
design <- read.table(
        file="designMatrix.txt",
        header=TRUE,
        sep="\t",
        row.names=1
        )

## Differential expression analysis (Untreated vs Treated)
this_analysis <- "Untreated_vs_Treated"
this_data <- data[,c(2:6)]
this_design <- design[c(1:5),]
	dds = DESeqDataSetFromMatrix(
		countData=this_data,
		colData=this_design,
		design= ~ Condition
		)
	dds <- dds[rowSums(counts(dds)) > 100,]	## Remove the low read count genes
	dds <- DESeq(dds)
	res <- results(dds,contrast=c("Condition", "Treated", "Untreated"))
	
	pdf(paste(this_analysis, ".pdf", sep=""), width=7, height=5)
	hist(res$pvalue, breaks=100, main="", col="skyblue", border="slateblue")
	plotDispEsts( dds )
	plotMA(res)
	dev.off()
	
## Writing outputs		
res <- as.data.frame(res)
res$Gene <- row.names(res)	
res$padj.BH <- p.adjust(res$pvalue, method="BH") 
res <- merge(res, data, by="Gene")
res <- merge(res, tpm, by="Gene")
res <- res[order(res$padj.BH),]
write.table(res, file = paste(this_analysis, "_deseq2_results_all_genes.txt", sep=""), sep = "\t", row.names=F, quote=F)
res_sig <- subset(res, res$padj.BH < 0.05 & (res$log2FoldChange > 0.5 | res$log2FoldChange < -0.5))
write.table(res_sig, file = paste(this_analysis, "_deseq2_results_significant_genes.txt", sep=""), sep = "\t", row.names=F, quote=F)
write.table(res_sig[,c(1,3)], file = paste(this_analysis, "_deseq2_results_significant_genes.rnk", sep=""), sep = "\t", row.names=F, quote=F, col.names=F)
write.table(res_sig[,c(1,17:22)], file = paste(this_analysis, "_input_for_heatmap.txt", sep=""), sep = "\t", row.names=F, quote=F, col.names=T)

## Volcano plot with color labelling of TGF-pathway genes ##
pdf(paste(this_analysis, "_volcano.pdf", sep=""), width=8, height=8)
with(res, plot(log2FoldChange, -log10(padj.BH), pch=20, ylim=c(0,50), cex=0.5, col="red", xlab="Log2 fold change", ylab="-log10 Padj value", cex.lab=1.5, cex.axis=1.5, cex.names=1.5))
with(subset(res, padj.BH >= 0.05 | (res$log2FoldChange < 0.5 & res$log2FoldChange > -0.5)), points(log2FoldChange, -log10(padj.BH), pch=20, col="gray", cex=0.5))
with(subset(res, Gene == "FOXP3"), points(log2FoldChange, -log10(padj.BH), pch=20, col="green", cex=1.5))

tgfb_genes <- c("FURIN", "PMEPA1", "SMAD7", "UBE2M", "SKI", "FKBP1A", "STRAP", "TGFB1")
with(subset(res, Gene %in% tgfb_genes), points(log2FoldChange, -log10(padj.BH), pch=20, col="blue", cex=1.5))
legend("topright", legend=c("TGF-beta signaling genes", "FOXP3", "Significant DEGs"), pch=20, col=c("blue", "green", "red"), cex=0.9)

abline(h=1.30103, col = "blue",lty=2)
abline(v=0.5, col = "blue",lty=2)
abline(v=-0.5, col = "blue",lty=2)
dev.off()
