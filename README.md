# at-pst-rnaseq

Snakemake RNA-seq pipeline for differential expression in *Arabidopsis thaliana* infected with *Pseudomonas syringae* DC3000 vs mock, from raw reads to DESeq2 results.

This repository reanalyses public RNA-seq data from Vogel et al. 2016 (ENA PRJEB13938, ArrayExpress E-MTAB-4683), comparing *Arabidopsis thaliana* (Col-0) leaves infected with *Pseudomonas syringae* pv. *tomato* DC3000 against mock-treated controls.

The Snakemake workflow covers read QC and trimming (fastp), alignment (STAR), count generation (featureCounts), and differential expression with DESeq2 in R. Reads are aligned to the TAIR10 assembly with Araport11 gene models (Ensembl Plants release 63).

Outputs include normalised counts, DESeq2 result tables with shrunken log2 fold changes, and QC and diagnostic plots (PCA, sample distances, MA plot, volcano plot).

## Data

The study has 48 single-end libraries: 3 colonisation conditions × mock / Pst DC3000 × 2 time points (2 and 7 days post infection) × 4 replicates. `config/samples.tsv` lists every run with its ENA FASTQ URL and MD5 (built from the E-MTAB-4683 SDRF and the ENA file report).

By default only the plain comparison is analysed (`colonisation: [control]`, plants in 10 mM MgCl2 without a commensal): 16 libraries, giving two contrasts, Pst vs mock at 2 dpi and at 7 dpi. Add `M_extorquens` and/or `S_melonis` to `colonisation` in `config/config.yaml` to include the commensal-colonised plants; each condition gets its own Pst-vs-mock contrasts.

## Usage

Requires [conda](https://docs.conda.io) and Snakemake ≥ 8 (tested with 9.26). All tools are installed per rule from `workflow/envs/`; `profiles/default/config.yaml` is picked up automatically and enables conda deployment and the greedy scheduler.

```bash
snakemake --cores 8
```

A full run of the default 16 libraries downloads the FASTQ files (about 1 GB each, deleted after trimming), builds the STAR index, and needs roughly 4 GB RAM for STAR.

## Outputs

| Path | Content |
| --- | --- |
| `results/qc/fastp/` | per-sample fastp HTML / JSON reports |
| `results/star/<sample>/` | sorted BAM and STAR `Log.final.out` |
| `results/counts/counts.tsv` | featureCounts gene × sample raw count matrix (+ `.summary`) |
| `results/deseq2/normalised_counts.tsv` | median-of-ratios normalised counts |
| `results/deseq2/vst_counts.tsv` | variance-stabilised counts (used for PCA / distances) |
| `results/deseq2/Pst_vs_mock_<colonisation>_d<dpi>.tsv` | DESeq2 results; `log2FoldChange` is ashr-shrunken, `log2FoldChange_unshrunk` the raw estimate |
| `results/plots/` | `pca.pdf`, `sample_distances.pdf`, `ma_*.pdf`, `volcano_*.pdf` |

## Validation

Run on all 48 libraries (macOS, Apple Silicon) and, with 1M-read truncated inputs, on Linux (aarch64 container). Cross-check against the independent count table for the same project on Zenodo (record 3627021: HISAT2 + featureCounts on Araport11), re-running the same per-timepoint `~ infection` DESeq2 model on both count matrices for the control plants:

| Contrast | log2FC Pearson r | DE genes (ours / Zenodo) | Overlap (Jaccard) |
| --- | --- | --- | --- |
| Pst vs mock, 2 dpi | 0.951 | 381 / 379 | 0.919 |
| Pst vs mock, 7 dpi | 0.974 | 1765 / 1757 | 0.894 |

Per-sample gene-count Spearman correlation between the two matrices is 0.94. The Zenodo record contains counts and a design table only, not a published differential gene list, so no such list is compared.

Pipeline DE counts depend on which libraries are analysed because the shared `~ group` model pools dispersion estimates (e.g. control 2 dpi: 305 DE genes with the 16 control libraries, 278 with all 48).

## Caveats

- **macOS:** the conda/upstream STAR 2.7.11b maps 0 reads on macOS, because it aliases its read buffers with `std::stringbuf::pubsetbuf()` and libc++ ignores that call (alexdobin/STAR#2663 shows the related spawn error). On macOS the workflow therefore builds a patched STAR (`workflow/patches/star-libcxx-streambuf.patch`) into `resources/bin/STAR` and uses it for indexing and alignment. This needs the Xcode command line tools and `brew install libomp`. The build omits STAR's SIMD adapter-clipping code (`--clip3pAdapterSeq`), which the workflow does not use. On Linux the conda STAR is used unchanged.
- The SDRF records library strandedness as "not applicable". `featurecounts.strandedness` defaults to 0 (unstranded); check the assigned fraction in `results/counts/counts.tsv.summary` and rerun with 1 or 2 if it is low.
- The greedy scheduler is set in the default profile because Snakemake's default ILP scheduler crashes on Apple Silicon (`Bad CPU type in executable` from the bundled x86 CBC solver).
- Shrinkage uses `ashr` (not `apeglm`) because the design is `~ group` with contrasts between group levels, which `apeglm` does not support.
