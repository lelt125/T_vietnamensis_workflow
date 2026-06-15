from pathlib import Path
import pandas as pd

RUN_NAME = config["run_name"]
RESULTS_BASE = config["results_base"]
RESULTS_DIR = f"{RESULTS_BASE}/{RUN_NAME}"

ASSEMBLIES = pd.read_csv(config["assemblies_csv"], comment="#")

SAMPLES = ASSEMBLIES["sample_id"].tolist()
FASTA_BY_SAMPLE = dict(zip(ASSEMBLIES["sample_id"], ASSEMBLIES["fasta"]))

missing_fastas = [
    fasta for fasta in FASTA_BY_SAMPLE.values()
    if not Path(fasta).exists()
]

if missing_fastas:
    raise ValueError(
        "Some assembly FASTA files do not exist:\n" +
        "\n".join(missing_fastas)
    )

HIFI_READS_DIR = Path(config["hifi_reads_dir"]).resolve()

HIFI_READS = sorted(
    HIFI_READS_DIR.glob("*.hifi_reads.fastq.gz")
)

if len(HIFI_READS) == 0:
    raise ValueError(f"No HiFi read files found in {HIFI_READS_DIR}")

HIFI_READ_BY_CHUNK = {
    read.name.replace(".hifi_reads.fastq.gz", ""): str(read.resolve())
    for read in HIFI_READS
}

HIFI_CHUNKS = list(HIFI_READ_BY_CHUNK.keys())

BUSCO_LINEAGE = config["busco"]["lineage"]


rule all:
    input:
        expand(
            f"{RESULTS_DIR}/genome_qc/seqkit/{{sample}}.seqkit_stats.tsv",
            sample=SAMPLES
        ),
        expand(
            f"{RESULTS_DIR}/genome_qc/seqkit/{{sample}}.seqkit_fx2tab.tsv",
            sample=SAMPLES
        ),
        expand(
            f"{RESULTS_DIR}/genome_qc/busco/{{sample}}_busco",
            sample=SAMPLES
        ),
        expand(
            f"{RESULTS_DIR}/genome_qc/merqury/{{sample}}/{{sample}}.qv",
            sample=SAMPLES
        ),
        expand(
            f"{RESULTS_DIR}/genome_qc/merqury/{{sample}}/{{sample}}.completeness.stats",
            sample=SAMPLES
        )


rule seqkit_assembly_stats:
    input:
        fasta=lambda wildcards: FASTA_BY_SAMPLE[wildcards.sample]
    output:
        stats=f"{RESULTS_DIR}/genome_qc/seqkit/{{sample}}.seqkit_stats.tsv",
        fx2tab=f"{RESULTS_DIR}/genome_qc/seqkit/{{sample}}.seqkit_fx2tab.tsv"
    threads:
        config["seqkit"]["threads"]
    shell:
        r"""
        seqkit stats -j {threads} -a -T {input.fasta} > {output.stats}

        seqkit fx2tab -j {threads} \
            --name \
            --length \
            --gc \
            {input.fasta} > {output.fx2tab}
        """


rule prepare_busco_fasta:
    input:
        fasta=lambda wildcards: FASTA_BY_SAMPLE[wildcards.sample]
    output:
        fasta=f"{RESULTS_DIR}/genome_qc/busco/prepared_fasta/{{sample}}.busco_input.fa"
    threads:
        8
    shell:
        r"""
        seqkit seq -j {threads} -w 80 {input.fasta} > {output.fasta}
        """


rule busco_genome:
    input:
        fasta=f"{RESULTS_DIR}/genome_qc/busco/prepared_fasta/{{sample}}.busco_input.fa"
    output:
        directory(f"{RESULTS_DIR}/genome_qc/busco/{{sample}}_busco")
    threads:
        config["busco"]["threads"]
    resources:
        mem_gb=config["busco"]["mem_gb"]
    params:
        out_parent=f"{RESULTS_DIR}/genome_qc/busco",
        lineage=BUSCO_LINEAGE,
        download_path=lambda wildcards: str(Path(config["busco"]["download_path"]).resolve()),
        java_heap=config["busco"]["java_heap"]
    shell:
        r"""
        export _JAVA_OPTIONS="-Xmx{params.java_heap}"

        busco \
            -i {input.fasta} \
            -o {wildcards.sample}_busco \
            --out_path {params.out_parent} \
            -m genome \
            -l {params.lineage} \
            -c {threads} \
            --download_path {params.download_path} \
            --force
        """


rule meryl_count_hifi_chunk:
    input:
        hifi=lambda wildcards: HIFI_READ_BY_CHUNK[wildcards.chunk]
    output:
        meryl=directory(
            f"{RESULTS_DIR}/genome_qc/meryl/chunks/{{chunk}}.k{config['meryl']['kmer']}.meryl"
        )
    threads:
        config["meryl"]["count_threads"]
    resources:
        mem_gb=config["meryl"]["count_mem_gb"]
    params:
        kmer=config["meryl"]["kmer"]
    log:
        f"{RESULTS_DIR}/logs/meryl/chunks/{{chunk}}.meryl_count.log"
    shell:
        r"""
        mkdir -p $(dirname {log})

        meryl k={params.kmer} count \
            threads={threads} \
            memory={resources.mem_gb}g \
            {input.hifi} \
            output {output.meryl} \
            > {log} 2>&1
        """


rule meryl_merge_hifi_chunks:
    input:
        expand(
            f"{RESULTS_DIR}/genome_qc/meryl/chunks/{{chunk}}.k{config['meryl']['kmer']}.meryl",
            chunk=HIFI_CHUNKS
        )
    output:
        meryl=directory(
            f"{RESULTS_DIR}/genome_qc/meryl/hifi.combined.k{config['meryl']['kmer']}.meryl"
        )
    threads:
        config["meryl"]["merge_threads"]
    resources:
        mem_gb=config["meryl"]["merge_mem_gb"]
    log:
        f"{RESULTS_DIR}/logs/meryl/hifi.combined.meryl_merge.log"
    shell:
        r"""
        mkdir -p $(dirname {log})

        meryl union-sum \
            threads={threads} \
            memory={resources.mem_gb}g \
            {input} \
            output {output.meryl} \
            > {log} 2>&1
        """

rule merqury_assembly_qc:
    input:
        reads_meryl=f"{RESULTS_DIR}/genome_qc/meryl/hifi.combined.k{config['meryl']['kmer']}.meryl",
        assembly=lambda wildcards: FASTA_BY_SAMPLE[wildcards.sample]
    output:
        qv=f"{RESULTS_DIR}/genome_qc/merqury/{{sample}}/{{sample}}.qv",
        completeness=f"{RESULTS_DIR}/genome_qc/merqury/{{sample}}/{{sample}}.completeness.stats"
    threads:
        config["merqury"]["threads"]
    resources:
        mem_gb=config["merqury"]["mem_gb"]
    params:
        outdir=lambda wildcards: f"{RESULTS_DIR}/genome_qc/merqury/{wildcards.sample}",
        prefix=lambda wildcards: wildcards.sample,
        java_heap=config["merqury"].get("java_heap", "64g")
    log:
        f"{RESULTS_DIR}/logs/merqury/{{sample}}.merqury.log"
    shell:
        r"""
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})

        READS_MERYL=$(realpath {input.reads_meryl})
        ASSEMBLY=$(realpath {input.assembly})
        LOG=$(realpath {log})

        export _JAVA_OPTIONS="-Xmx{params.java_heap}"

        cd {params.outdir}

        mkdir -p logs

        merqury.sh \
            "$READS_MERYL" \
            "$ASSEMBLY" \
            {params.prefix} \
            > "$LOG" 2>&1
        """