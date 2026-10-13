import os
import sys

import pandas as pd

configfile: "config/config.yaml"


samples = pd.read_csv(config["samples"], sep="\t", dtype=str)
samples = samples[samples["colonisation"].isin(config["colonisation"])].set_index(
    "sample", drop=False
)
if samples.empty:
    raise ValueError(f"No samples match colonisation={config['colonisation']}")

SAMPLES = list(samples.index)

# One Pst-vs-mock contrast per (colonisation, dpi) present in the selected samples.
CONTRASTS = sorted(
    {f"{r.colonisation}_d{r.dpi}" for r in samples.itertuples()},
    key=lambda c: (c.rsplit("_d", 1)[0], int(c.rsplit("_d", 1)[1])),
)


# On macOS the conda STAR maps 0 reads (libc++ ignores stringbuf::pubsetbuf); build a patched
# binary instead (see workflow/rules/align.smk). On Linux the conda STAR is used as is.
if sys.platform == "darwin":
    STAR = "resources/bin/STAR"
    STAR_DEPS = ["resources/bin/STAR"]
else:
    STAR = "STAR"
    STAR_DEPS = []


include: "workflow/rules/reference.smk"
include: "workflow/rules/reads.smk"
include: "workflow/rules/align.smk"
include: "workflow/rules/deseq2.smk"


rule all:
    input:
        "results/deseq2/normalised_counts.tsv",
        expand("results/deseq2/Pst_vs_mock_{contrast}.tsv", contrast=CONTRASTS),
        "results/plots/pca.pdf",
        "results/plots/sample_distances.pdf",
        expand("results/plots/ma_{contrast}.pdf", contrast=CONTRASTS),
        expand("results/plots/volcano_{contrast}.pdf", contrast=CONTRASTS),
