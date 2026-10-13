rule deseq2:
    input:
        counts="results/counts/counts.tsv",
    output:
        normalised="results/deseq2/normalised_counts.tsv",
        vst="results/deseq2/vst_counts.tsv",
        results=expand("results/deseq2/Pst_vs_mock_{contrast}.tsv", contrast=CONTRASTS),
    params:
        samples=config["samples"],
        colonisation=config["colonisation"],
        alpha=config["deseq2"]["alpha"],
    conda:
        "../envs/r.yaml"
    log:
        "logs/deseq2.log",
    script:
        "../scripts/deseq2.R"


rule plots:
    input:
        vst="results/deseq2/vst_counts.tsv",
        results=expand("results/deseq2/Pst_vs_mock_{contrast}.tsv", contrast=CONTRASTS),
    output:
        pca="results/plots/pca.pdf",
        distances="results/plots/sample_distances.pdf",
        ma=expand("results/plots/ma_{contrast}.pdf", contrast=CONTRASTS),
        volcano=expand("results/plots/volcano_{contrast}.pdf", contrast=CONTRASTS),
    params:
        samples=config["samples"],
        colonisation=config["colonisation"],
        alpha=config["deseq2"]["alpha"],
        lfc=config["deseq2"]["lfc_threshold"],
    conda:
        "../envs/r.yaml"
    log:
        "logs/plots.log",
    script:
        "../scripts/plots.R"
