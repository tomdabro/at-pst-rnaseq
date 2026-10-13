# macOS only: build STAR with the libc++ stream-buffer fix (workflow/patches/).
rule build_star_macos:
    output:
        "resources/bin/STAR",
    log:
        "logs/build_star_macos.log",
    shell:
        "workflow/scripts/build_star_macos.sh {output} > {log} 2>&1"


rule star_align:
    input:
        fastq="results/trimmed/{sample}.fastq.gz",
        index="resources/star_index",
        star=STAR_DEPS,
    output:
        bam="results/star/{sample}/Aligned.sortedByCoord.out.bam",
        log_final="results/star/{sample}/Log.final.out",
    params:
        star=STAR,
        prefix=lambda wc, output: os.path.dirname(output.bam) + "/",
    threads: 8
    conda:
        "../envs/star.yaml"
    log:
        "logs/star/{sample}.log",
    shell:
        """
        {params.star} --runThreadN {threads} --genomeDir {input.index} \
            --readFilesIn {input.fastq} --readFilesCommand gunzip -c \
            --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {params.prefix} \
            --outFilterMultimapNmax 20 > {log} 2>&1
        """


# One featureCounts call over all BAMs gives a single, consistently ordered count matrix.
rule featurecounts:
    input:
        bams=expand("results/star/{sample}/Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        gtf="resources/annotation.gtf",
    output:
        counts="results/counts/counts.tsv",
        summary="results/counts/counts.tsv.summary",
    params:
        strand=config["featurecounts"]["strandedness"],
    threads: 8
    conda:
        "../envs/subread.yaml"
    log:
        "logs/featurecounts.log",
    shell:
        """
        featureCounts -T {threads} -s {params.strand} -t exon -g gene_id \
            -a {input.gtf} -o {output.counts} {input.bams} > {log} 2>&1
        """
