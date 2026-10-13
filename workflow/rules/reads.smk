import hashlib


rule download_fastq:
    output:
        temp("resources/fastq/{sample}.fastq.gz"),
    params:
        # ENA serves the same path over HTTPS; parallel FTP connections get refused.
        url=lambda wc: samples.loc[wc.sample, "fastq_url"].replace("ftp://", "https://"),
        md5=lambda wc: samples.loc[wc.sample, "md5"],
    log:
        "logs/download/{sample}.log",
    resources:
        downloads=1,
    retries: 2
    run:
        shell("curl -sSfL --retry 5 --retry-delay 10 --retry-all-errors -o {output} {params.url} 2> {log}")
        h = hashlib.md5()
        with open(output[0], "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        if h.hexdigest() != params.md5:
            raise ValueError(f"MD5 mismatch for {wildcards.sample}: {h.hexdigest()} != {params.md5}")


# fastp does QC and trimming in one pass (single-end reads; adapters auto-detected).
rule fastp:
    input:
        "resources/fastq/{sample}.fastq.gz",
    output:
        fastq=temp("results/trimmed/{sample}.fastq.gz"),
        html="results/qc/fastp/{sample}.html",
        json="results/qc/fastp/{sample}.json",
    threads: 4
    conda:
        "../envs/fastp.yaml"
    log:
        "logs/fastp/{sample}.log",
    shell:
        """
        fastp -i {input} -o {output.fastq} --html {output.html} --json {output.json} \
            --thread {threads} --qualified_quality_phred 20 --length_required 36 > {log} 2>&1
        """
