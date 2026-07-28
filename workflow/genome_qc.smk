include: "rules/common_config_files.smk"

RUN_NAME, RESULTS_BASE, RESULTS_DIR = get_results_dir(config)
ASSEMBLIES, SAMPLES, FASTA_BY_SAMPLE = load_assemblies(config)
HIFI_READS_DIR, HIFI_READS, HIFI_READ_BY_CHUNK, HIFI_CHUNKS = load_hifi_reads(config)

VALID_BUSCO_LINEAGES, DEFAULT_BUSCO_LINEAGE = load_valid_busco_lineages(config)
BUSCO_LINEAGE_BY_SAMPLE = build_busco_lineage_by_sample(
    ASSEMBLIES,
    VALID_BUSCO_LINEAGES,
    DEFAULT_BUSCO_LINEAGE
)

SEQKIT_SAMPLES = get_qc_samples(ASSEMBLIES, "run_seqkit", default=True)
BUSCO_SAMPLES = get_qc_samples(ASSEMBLIES, "run_busco", default=False)
MERQURY_SAMPLES = get_qc_samples(ASSEMBLIES, "run_merqury", default=False)
MOSDEPTH_SAMPLES = get_qc_samples(ASSEMBLIES, "run_mosdepth", default=False)

SEQKIT_STATS_OUTPUTS = expand(
    f"{RESULTS_DIR}/genome_qc/seqkit/{{sample}}.seqkit_stats.tsv",
    sample=SEQKIT_SAMPLES
)

SEQKIT_FX2TAB_OUTPUTS = expand(
    f"{RESULTS_DIR}/genome_qc/seqkit/{{sample}}.seqkit_fx2tab.tsv",
    sample=SEQKIT_SAMPLES
)

BUSCO_OUTPUTS = [
    f"{RESULTS_DIR}/genome_qc/busco/{sample}/{BUSCO_LINEAGE_BY_SAMPLE[sample]}"
    for sample in BUSCO_SAMPLES
]

MERQURY_QV_OUTPUTS = expand(
    f"{RESULTS_DIR}/genome_qc/merqury/{{sample}}/{{sample}}.qv",
    sample=MERQURY_SAMPLES
)

MERQURY_COMPLETENESS_OUTPUTS = expand(
    f"{RESULTS_DIR}/genome_qc/merqury/{{sample}}/{{sample}}.completeness.stats",
    sample=MERQURY_SAMPLES
)

MOSDEPTH_SUMMARY_OUTPUTS = expand(
    f"{RESULTS_DIR}/genome_qc/mosdepth/{{sample}}/{{sample}}.hifi.mosdepth.summary.txt",
    sample=MOSDEPTH_SAMPLES
)

MOSDEPTH_GLOBAL_DIST_OUTPUTS = expand(
    f"{RESULTS_DIR}/genome_qc/mosdepth/{{sample}}/{{sample}}.hifi.mosdepth.global.dist.txt",
    sample=MOSDEPTH_SAMPLES
)

rule all:
    input:
        SEQKIT_STATS_OUTPUTS,
        SEQKIT_FX2TAB_OUTPUTS,
        BUSCO_OUTPUTS,
        MERQURY_QV_OUTPUTS,
        MERQURY_COMPLETENESS_OUTPUTS,
        MOSDEPTH_SUMMARY_OUTPUTS,
        MOSDEPTH_GLOBAL_DIST_OUTPUTS


rule seqkit_assembly_stats:
    input:
        fasta=lambda wildcards: FASTA_BY_SAMPLE[wildcards.sample]
    output:
        stats=f"{RESULTS_DIR}/genome_qc/seqkit/{{sample}}.seqkit_stats.tsv",
        fx2tab=f"{RESULTS_DIR}/genome_qc/seqkit/{{sample}}.seqkit_fx2tab.tsv"
    threads:
        config["seqkit"]["threads"]
    log:
        stats=f"{RESULTS_DIR}/logs/genome_qc/seqkit/{{sample}}.seqkit_stats.log",
        fx2tab=f"{RESULTS_DIR}/logs/genome_qc/seqkit/{{sample}}.seqkit_fx2tab.log"
    benchmark:
        f"{RESULTS_DIR}/benchmarks/genome_qc/seqkit/{{sample}}.tsv"
    shell:
        r"""
        mkdir -p $(dirname {output.stats})
        mkdir -p $(dirname {log.stats})

        seqkit stats -j {threads} -a -T {input.fasta} \
            > {output.stats} \
            2> {log.stats}

        seqkit fx2tab -j {threads} \
            --name \
            --length \
            --gc \
            {input.fasta} \
            > {output.fx2tab} \
            2> {log.fx2tab}
        """


rule prepare_busco_fasta:
    input:
        fasta=lambda wildcards: FASTA_BY_SAMPLE[wildcards.sample]
    output:
        fasta=f"{RESULTS_DIR}/genome_qc/busco_input/{{sample}}.busco_input.fa"
    threads:
        8
    log:
        f"{RESULTS_DIR}/logs/genome_qc/busco/prepare_fasta/{{sample}}.log"
    benchmark:
        f"{RESULTS_DIR}/benchmarks/genome_qc/busco/prepare_fasta/{{sample}}.tsv"
    shell:
        r"""
        mkdir -p $(dirname {output.fasta})
        mkdir -p $(dirname {log})

        seqkit seq -j {threads} -w 80 {input.fasta} \
            > {output.fasta} \
            2> {log}
        """


rule busco_genome:
    input:
        fasta=f"{RESULTS_DIR}/genome_qc/busco_input/{{sample}}.busco_input.fa"
    output:
        directory(f"{RESULTS_DIR}/genome_qc/busco/{{sample}}/{{lineage}}")
    threads:
        config["busco"]["threads"]
    resources:
        mem_gb=config["busco"]["mem_gb"]
    params:
        out_parent=lambda wildcards: f"{RESULTS_DIR}/genome_qc/busco/{wildcards.sample}",
        download_path=lambda wildcards: str(Path(config["busco"]["download_path"]).resolve()),
        java_heap=config["busco"]["java_heap"],
        extra=config["busco"].get("extra", "")
    log:
        f"{RESULTS_DIR}/logs/genome_qc/busco/{{sample}}/{{lineage}}.log"
    benchmark:
        f"{RESULTS_DIR}/benchmarks/genome_qc/busco/{{sample}}/{{lineage}}.tsv"
    shell:
        r"""
        mkdir -p {params.out_parent}
        mkdir -p $(dirname {log})

        export _JAVA_OPTIONS="-Xmx{params.java_heap}"

        busco \
            -i {input.fasta} \
            -o {wildcards.lineage} \
            -l {wildcards.lineage} \
            --out_path {params.out_parent} \
            -m genome \
            -c {threads} \
            --download_path {params.download_path} \
            --force \
            {params.extra} \
            > {log} 2>&1
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
        f"{RESULTS_DIR}/logs/genome_qc/meryl/chunks/{{chunk}}.meryl_count.log"
    benchmark:
        f"{RESULTS_DIR}/benchmarks/genome_qc/meryl/chunks/{{chunk}}.tsv"
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
        f"{RESULTS_DIR}/logs/genome_qc/meryl/hifi.combined.meryl_merge.log"
    benchmark:
        f"{RESULTS_DIR}/benchmarks/genome_qc/meryl/hifi.combined.merge.tsv"
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
        f"{RESULTS_DIR}/logs/genome_qc/merqury/{{sample}}.log"
    benchmark:
        f"{RESULTS_DIR}/benchmarks/genome_qc/merqury/{{sample}}.tsv"
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

rule map_hifi_chunk_to_assembly_bam:
    input:
        assembly=lambda wildcards: FASTA_BY_SAMPLE[wildcards.sample],
        hifi=lambda wildcards: HIFI_READ_BY_CHUNK[wildcards.chunk]
    output:
        bam=temp(f"{RESULTS_DIR}/read_mapping/hifi_bam/{{sample}}/chunks/{{chunk}}.sorted.bam")
    threads:
        config["minimap2"]["threads"] + config["minimap2"].get("sort_threads", 4)
    resources:
        mem_gb=config["minimap2"]["mem_gb"]
    params:
        minimap2_threads=lambda wildcards, threads: max(
            1,
            threads - config["minimap2"].get("sort_threads", 4)
        ),
        sort_threads=config["minimap2"].get("sort_threads", 4),
        sort_mem=config["minimap2"].get("sort_mem", "2G"),
        index_size=config["minimap2"].get("index_size", "4G"),
        split_prefix=lambda wildcards: f"{RESULTS_DIR}/read_mapping/hifi_bam/{wildcards.sample}/chunks/{wildcards.chunk}.split_idx"
    log:
        minimap2=f"{RESULTS_DIR}/logs/read_mapping/minimap2/hifi_bam/{{sample}}/{{chunk}}.minimap2.log",
        sort=f"{RESULTS_DIR}/logs/read_mapping/samtools_sort/hifi_bam/{{sample}}/{{chunk}}.sort.log"
    benchmark:
        f"{RESULTS_DIR}/benchmarks/read_mapping/hifi_bam/{{sample}}/{{chunk}}.tsv"
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        mkdir -p $(dirname {log.minimap2})
        mkdir -p $(dirname {log.sort})

        minimap2 \
            -t {params.minimap2_threads} \
            -I {params.index_size} \
            --split-prefix {params.split_prefix} \
            -ax map-hifi \
            {input.assembly} \
            {input.hifi} \
            2> {log.minimap2} \
        | samtools sort \
            -@ {params.sort_threads} \
            -m {params.sort_mem} \
            -o {output.bam} \
            - \
            2> {log.sort}

        test -s {output.bam}
        samtools quickcheck {output.bam}
        """


rule merge_hifi_to_assembly_bam:
    input:
        bams=lambda wildcards: expand(
            f"{RESULTS_DIR}/read_mapping/hifi_bam/{wildcards.sample}/chunks/{{chunk}}.sorted.bam",
            chunk=HIFI_CHUNKS
        )
    output:
        bam=f"{RESULTS_DIR}/read_mapping/hifi_bam/{{sample}}/{{sample}}.hifi_to_assembly.sorted.bam",
        bai=f"{RESULTS_DIR}/read_mapping/hifi_bam/{{sample}}/{{sample}}.hifi_to_assembly.sorted.bam.bai"
    threads:
        config["minimap2"].get("merge_threads", 8)
    resources:
        mem_gb=config["minimap2"].get("merge_mem_gb", 32)
    params:
        bam_list=f"{RESULTS_DIR}/read_mapping/hifi_bam/{{sample}}/{{sample}}.bam_list.txt"
    log:
        merge=f"{RESULTS_DIR}/logs/read_mapping/samtools_merge/hifi_bam/{{sample}}.merge.log",
        index=f"{RESULTS_DIR}/logs/read_mapping/samtools_index/hifi_bam/{{sample}}.index.log"
    benchmark:
        f"{RESULTS_DIR}/benchmarks/read_mapping/hifi_bam/{{sample}}.merge.tsv"
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        mkdir -p $(dirname {log.merge})
        mkdir -p $(dirname {log.index})

        printf "%s\n" {input.bams} > {params.bam_list}

        samtools merge \
            -@ {threads} \
            -f \
            -b {params.bam_list} \
            {output.bam} \
            > {log.merge} 2>&1

        samtools index \
            -@ {threads} \
            {output.bam} \
            > {log.index} 2>&1

        test -s {output.bam}
        test -s {output.bai}
        samtools quickcheck {output.bam}
        """


rule mosdepth_hifi_coverage:
    input:
        bam=f"{RESULTS_DIR}/read_mapping/hifi_bam/{{sample}}/{{sample}}.hifi_to_assembly.sorted.bam",
        bai=f"{RESULTS_DIR}/read_mapping/hifi_bam/{{sample}}/{{sample}}.hifi_to_assembly.sorted.bam.bai"
    output:
        summary=f"{RESULTS_DIR}/genome_qc/mosdepth/{{sample}}/{{sample}}.hifi.mosdepth.summary.txt",
        global_dist=f"{RESULTS_DIR}/genome_qc/mosdepth/{{sample}}/{{sample}}.hifi.mosdepth.global.dist.txt"
    threads:
        config["mosdepth"]["threads"]
    resources:
        mem_gb=config["mosdepth"]["mem_gb"]
    params:
        outdir=lambda wildcards: f"{RESULTS_DIR}/genome_qc/mosdepth/{wildcards.sample}",
        prefix=lambda wildcards: f"{RESULTS_DIR}/genome_qc/mosdepth/{wildcards.sample}/{wildcards.sample}.hifi",
        extra=config["mosdepth"].get("extra", "--fast-mode --no-per-base")
    log:
        f"{RESULTS_DIR}/logs/genome_qc/mosdepth/{{sample}}.log"
    benchmark:
        f"{RESULTS_DIR}/benchmarks/genome_qc/mosdepth/{{sample}}.tsv"
    shell:
        r"""
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})

        mosdepth \
            -t {threads} \
            {params.extra} \
            {params.prefix} \
            {input.bam} \
            > {log} 2>&1

        test -s {output.summary}
        test -s {output.global_dist}
        """