log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages(library(DESeq2))
source(file.path(dirname(snakemake@scriptdir), "scripts", "common.R"))

coldata <- load_coldata(snakemake@params$samples, snakemake@params$colonisation)

# featureCounts columns are named after the BAM paths: results/star/<sample>/Aligned...bam
fc <- read.delim(snakemake@input$counts, comment.char = "#", check.names = FALSE)
counts <- as.matrix(fc[, -(1:6)])
rownames(counts) <- fc$Geneid
colnames(counts) <- basename(dirname(colnames(counts)))
counts <- counts[, rownames(coldata)]
stopifnot(!anyNA(counts))

dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~group)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

write_table <- function(df, path) write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)

write_table(data.frame(gene = rownames(dds), counts(dds, normalized = TRUE), check.names = FALSE),
            snakemake@output$normalised)
write_table(data.frame(gene = rownames(dds), assay(vst(dds, blind = TRUE)), check.names = FALSE),
            snakemake@output$vst)

alpha <- as.numeric(snakemake@params$alpha)
for (out in unlist(snakemake@output$results)) {
  ct <- parse_contrast(out)
  contrast <- c("group",
                paste(ct$colonisation, "Pst", paste0("d", ct$dpi), sep = "__"),
                paste(ct$colonisation, "mock", paste0("d", ct$dpi), sep = "__"))
  res <- results(dds, contrast = contrast, alpha = alpha)
  shrunk <- lfcShrink(dds, contrast = contrast, res = res, type = "ashr")
  write_table(data.frame(gene = rownames(shrunk),
                         baseMean = shrunk$baseMean,
                         log2FoldChange = shrunk$log2FoldChange,  # shrunken (ashr)
                         lfcSE = shrunk$lfcSE,
                         log2FoldChange_unshrunk = res$log2FoldChange,
                         pvalue = shrunk$pvalue,
                         padj = shrunk$padj),
              out)
}
