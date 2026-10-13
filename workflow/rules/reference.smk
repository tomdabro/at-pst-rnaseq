rule get_reference:
    output:
        fasta="resources/genome.fa",
        gtf="resources/annotation.gtf",
    params:
        fasta_url=config["reference"]["fasta"],
        gtf_url=config["reference"]["gtf"],
    log:
        "logs/get_reference.log",
    shell:
        """
        (curl -sSfL {params.fasta_url} | gunzip -c > {output.fasta}
         curl -sSfL {params.gtf_url} | gunzip -c > {output.gtf}) 2> {log}
        """


rule star_index:
    input:
        fasta="resources/genome.fa",
        gtf="resources/annotation.gtf",
        star=STAR_DEPS,
    output:
        directory("resources/star_index"),
    params:
        star=STAR,
        overhang=config["star"]["sjdb_overhang"],
        sa_nbases=config["star"]["genome_sa_index_nbases"],
    threads: 8
    conda:
        "../envs/star.yaml"
    log:
        "logs/star_index.log",
    shell:
        """
        mkdir -p {output}
        {params.star} --runMode genomeGenerate --runThreadN {threads} \
            --genomeDir {output} --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} --sjdbOverhang {params.overhang} \
            --genomeSAindexNbases {params.sa_nbases} \
            --outFileNamePrefix logs/star_index_ > {log} 2>&1
        """
