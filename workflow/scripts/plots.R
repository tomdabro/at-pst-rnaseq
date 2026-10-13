log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages({library(ggplot2); library(pheatmap)})
source(file.path(dirname(snakemake@scriptdir), "scripts", "common.R"))

coldata <- load_coldata(snakemake@params$samples, snakemake@params$colonisation)
alpha <- as.numeric(snakemake@params$alpha)
lfc_thr <- as.numeric(snakemake@params$lfc)

vst <- read.delim(snakemake@input$vst, check.names = FALSE, row.names = 1)
vst <- as.matrix(vst[, rownames(coldata)])

# PCA on the 500 most variable genes.
top <- head(order(apply(vst, 1, var), decreasing = TRUE), 500)
pca <- prcomp(t(vst[top, ]))
pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], coldata)
p <- ggplot(df, aes(PC1, PC2, colour = infection, shape = dpi)) +
  geom_point(size = 3) +
  labs(x = sprintf("PC1 (%.1f%%)", pct[1]), y = sprintf("PC2 (%.1f%%)", pct[2]),
       shape = "dpi") +
  theme_bw()
if (length(unique(coldata$colonisation)) > 1) p <- p + facet_wrap(~colonisation)
ggsave(snakemake@output$pca, p, width = 6, height = 5)

# Sample-to-sample Euclidean distances on vst values.
d <- as.matrix(dist(t(vst)))
ann <- coldata[, c("colonisation", "infection", "dpi")]
pdf(snakemake@output$distances, width = 8, height = 7)
pheatmap(d, annotation_col = ann, annotation_row = ann, show_colnames = FALSE)
dev.off()

classify <- function(res) {
  ifelse(!is.na(res$padj) & res$padj < alpha & abs(res$log2FoldChange) >= lfc_thr,
         ifelse(res$log2FoldChange > 0, "up", "down"), "n.s.")
}
cols <- c(up = "firebrick", down = "steelblue", n.s. = "grey70")

for (path in unlist(snakemake@input$results)) {
  ct <- parse_contrast(path)
  suffix <- paste0(ct$colonisation, "_d", ct$dpi)
  res <- read.delim(path)
  res$status <- classify(res)
  title <- sprintf("Pst vs mock, %s, %s dpi", ct$colonisation, ct$dpi)

  ma <- ggplot(res[res$baseMean > 0, ], aes(baseMean, log2FoldChange, colour = status)) +
    geom_point(size = 0.6) + scale_x_log10() + scale_colour_manual(values = cols) +
    geom_hline(yintercept = 0) + labs(title = title, y = "shrunken log2 fold change") +
    theme_bw()
  ggsave(grep(paste0("ma_", suffix, ".pdf"), unlist(snakemake@output$ma), value = TRUE, fixed = TRUE),
         ma, width = 6, height = 5)

  vol <- ggplot(res[!is.na(res$padj), ], aes(log2FoldChange, -log10(padj), colour = status)) +
    geom_point(size = 0.6) + scale_colour_manual(values = cols) +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed") +
    labs(title = title, x = "shrunken log2 fold change") + theme_bw()
  ggsave(grep(paste0("volcano_", suffix, ".pdf"), unlist(snakemake@output$volcano), value = TRUE, fixed = TRUE),
         vol, width = 6, height = 5)
}
