# filter_contaminant_reads.smk
include: "rules/common_config_files.smk"


RUN_NAME, RESULTS_BASE, RESULTS_DIR = get_results_dir(config)
HIFI_READS_DIR, HIFI_READS, HIFI_READ_BY_CHUNK, HIFI_CHUNKS = load_hifi_reads(config)

RAW_ASSEMBLY = config["contaminant_filter"]["raw_assembly_fasta"]
CONTAMINANT_CONTIGS_LIST = config["contaminant_filter"]["contaminant_contigs_list"]

FILTERED_HIFI_OUTPUTS = expand(
    f"{RESULTS_DIR}/contaminant_filter/clean_reads/{{chunk}}.clean.hifi_reads.fastq.gz",
    chunk=HIFI_CHUNKS,
)


rule all:
    input:
        FILTERED_HIFI_OUTPUTS,
        f"{RESULTS_DIR}/contaminant_filter/summary/contaminant_read_filtering.summary.tsv",


rule map_hifi_chunk_to_raw_assembly:
    input:
        assembly=RAW_ASSEMBLY,
        hifi=lambda wildcards: HIFI_READ_BY_CHUNK[wildcards.chunk],
    output:
        bam=temp(f"{RESULTS_DIR}/contaminant_filter/hifi_bam/{{chunk}}.sorted.bam"),
    log:
        minimap2=f"{RESULTS_DIR}/logs/contaminant_filter/minimap2/{{chunk}}.log",
        sort=f"{RESULTS_DIR}/logs/contaminant_filter/samtools_sort/{{chunk}}.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/contaminant_filter/hifi_bam/{{chunk}}.tsv"
    threads: config["minimap2"]["threads"] + config["minimap2"].get("sort_threads", 4)
    resources:
        mem_gb=config["minimap2"]["mem_gb"],
    params:
        minimap2_threads=lambda wildcards, threads: max(
            1, threads - config["minimap2"].get("sort_threads", 4)
        ),
        sort_threads=config["minimap2"].get("sort_threads", 4),
        sort_mem=config["minimap2"].get("sort_mem", "2G"),
        index_size=config["minimap2"].get("index_size", "4G"),
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        mkdir -p $(dirname {log.minimap2})
        mkdir -p $(dirname {log.sort})

        minimap2 \
            -t {params.minimap2_threads} \
            -I {params.index_size} \
            -ax map-hifi \
            {input.assembly} \
            {input.hifi} \
            2>{log.minimap2} \
            | samtools sort \
                -@ {params.sort_threads} \
                -m {params.sort_mem} \
                -o {output.bam} \
                - \
                2>{log.sort}

        test -s {output.bam}
        samtools quickcheck {output.bam}
        """


rule identify_contaminant_read_ids:
    input:
        bam=f"{RESULTS_DIR}/contaminant_filter/hifi_bam/{{chunk}}.sorted.bam",
        contaminant_list=CONTAMINANT_CONTIGS_LIST,
    output:
        read_ids=f"{RESULTS_DIR}/contaminant_filter/read_ids/{{chunk}}.contaminant_read_ids.txt",
    log:
        f"{RESULTS_DIR}/logs/contaminant_filter/identify_reads/{{chunk}}.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/contaminant_filter/identify_reads/{{chunk}}.tsv"
    threads: 4
    shell:
        r"""
        mkdir -p $(dirname {output.read_ids})
        mkdir -p $(dirname {log})

        # Primary alignments only (exclude unmapped/secondary/supplementary: 0x904),
        # keep read name (col 1) and reference contig (col 3),
        # retain reads whose primary alignment lands on a contaminant contig.
        samtools view -@ {threads} -F 0x904 {input.bam} \
            2>{log} \
            | awk 'BEGIN{{
                while ((getline line < "{input.contaminant_list}") > 0) contam[line]=1
              }}
              {{ if ($3 in contam) print $1 }}' \
            | sort -u \
                >{output.read_ids}

        touch {output.read_ids}
        """


rule filter_contaminant_reads_from_chunk:
    input:
        hifi=lambda wildcards: HIFI_READ_BY_CHUNK[wildcards.chunk],
        read_ids=f"{RESULTS_DIR}/contaminant_filter/read_ids/{{chunk}}.contaminant_read_ids.txt",
    output:
        clean=f"{RESULTS_DIR}/contaminant_filter/clean_reads/{{chunk}}.clean.hifi_reads.fastq.gz",
    log:
        f"{RESULTS_DIR}/logs/contaminant_filter/seqkit_grep/{{chunk}}.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/contaminant_filter/seqkit_grep/{{chunk}}.tsv"
    threads: config["seqkit"]["threads"]
    shell:
        r"""
        mkdir -p $(dirname {output.clean})
        mkdir -p $(dirname {log})

        if [ -s {input.read_ids} ]; then
            seqkit grep -j {threads} -v -f {input.read_ids} {input.hifi} \
                -o {output.clean} \
                2>{log}
        else
            cp {input.hifi} {output.clean}
            echo "No contaminant read IDs found; copied input unchanged." >{log}
        fi

        test -s {output.clean}
        """


rule summarize_contaminant_filtering:
    input:
        read_ids=expand(
            f"{RESULTS_DIR}/contaminant_filter/read_ids/{{chunk}}.contaminant_read_ids.txt",
            chunk=HIFI_CHUNKS,
        ),
    output:
        summary=f"{RESULTS_DIR}/contaminant_filter/summary/contaminant_read_filtering.summary.tsv",
    log:
        f"{RESULTS_DIR}/logs/contaminant_filter/summary.log",
    shell:
        r"""
        mkdir -p $(dirname {output.summary})
        mkdir -p $(dirname {log})

        echo -e "chunk\tn_contaminant_reads" >{output.summary}
        for f in {input.read_ids}; do
            chunk=$(basename "$f" .contaminant_read_ids.txt)
            n=$(wc -l <"$f")
            echo -e "${{chunk}}\t${{n}}" >>{output.summary}
        done
        """
