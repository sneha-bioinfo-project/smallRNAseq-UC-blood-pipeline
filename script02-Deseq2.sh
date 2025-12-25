#!/usr/bin/env Rscript

library(DESeq2)  # standard for small RNA differential expression. [web:607][web:610]

counts_file <- "../results/miRNA_counts.tsv"
out_all     <- "../results/DESeq2_miRNA_results.csv"
out_sig     <- "../results/significant_miRNAs_padj_lt_0.1.csv"
out_top20   <- "../results/top20_miRNAs_DESeq2.csv"
plot_dir    <- "../results/plots"

dir.create("../results", showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)

cts <- read.table(counts_file, header = TRUE, row.names = 1,
                  sep = "\t", check.names = FALSE)

## EDIT this to match your samples
## Example: 2 controls, 2 UC
condition <- factor(c("control","control","UC","UC"))
coldata <- data.frame(row.names = colnames(cts), condition = condition)

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData   = coldata,
                              design    = ~ condition)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition","UC","control"))
res <- lfcShrink(dds, coef = "condition_UC_vs_control", res = res,
                 type = "apeglm")   # LFC shrinkage is recommended. [web:607][web:591]

res_df <- as.data.frame(res[order(res$padj), ])

write.csv(res_df, out_all)

sig <- subset(res_df, !is.na(padj) & padj < 0.1)
write.csv(sig, out_sig)

top20 <- head(res_df, 20)
write.csv(top20, out_top20)

png(file.path(plot_dir, "MAplot_miRNA.png"), width = 1200, height = 1000, res = 150)
plotMA(res, main = "DESeq2 MA-plot", ylim = c(-5, 5))
dev.off()

png(file.path(plot_dir, "volcano_miRNA.png"), width = 1200, height = 1000, res = 150)
with(res_df, {
  plot(log2FoldChange, -log10(padj),
       pch = 20, cex = 0.6,
       xlab = "log2 fold change (UC vs control)",
       ylab = "-log10 adjusted p-value",
       main = "Volcano plot")
})
dev.off()

