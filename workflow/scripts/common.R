# Shared helpers for deseq2.R and plots.R.

load_coldata <- function(path, colonisation) {
  cd <- read.delim(path, stringsAsFactors = FALSE)
  cd <- cd[cd$colonisation %in% colonisation, ]
  rownames(cd) <- cd$sample
  cd$infection <- factor(cd$infection, levels = c("mock", "Pst"))
  cd$dpi <- factor(cd$dpi, levels = sort(unique(as.numeric(cd$dpi))))
  cd$group <- factor(paste(cd$colonisation, cd$infection, paste0("d", cd$dpi), sep = "__"))
  cd
}

# "Pst_vs_mock_control_d2.tsv" / "ma_control_d2.pdf" -> list(colonisation = "control", dpi = "2")
parse_contrast <- function(file) {
  m <- regmatches(basename(file), regexec("^(?:Pst_vs_mock|ma|volcano)_(.+)_d([0-9]+)\\.(?:tsv|pdf)$", basename(file), perl = TRUE))[[1]]
  if (length(m) != 3) stop("Cannot parse contrast from ", file)
  list(colonisation = m[2], dpi = m[3])
}
