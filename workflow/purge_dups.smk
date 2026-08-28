include: "rules/common_config_files.smk"


import re

RUN_NAME, RESULTS_BASE, RESULTS_DIR = get_results_dir(config)
ASSEMBLIES, SAMPLES, FASTA_BY_SAMPLE = load_assemblies(config)
HIFI_READS_DIR, HIFI_READS, HIFI_READ_BY_CHUNK, HIFI_CHUNKS = load_hifi_reads(config)

PURGE_TESTS = list(config["purge_dups_tests"].keys())
PRESETS = config["purge_dups"]["self_align_presets"]
ALIGNER = config["purge_dups"].get("self_align_aligner", "minimap2")

WRAPPER = "v9.15.0"


wildcard_constraints:
    sample="|".join(re.escape(s) for s in SAMPLES),
    preset="|".join(re.escape(p) for p in PRESETS),
    test_id="|".join(re.escape(t) for t in PURGE_TESTS),
    chunk="|".join(re.escape(c) for c in HIFI_CHUNKS),


PURGE_DUP_OUTPUTS = expand(
    f"{RESULTS_DIR}/purge_dups/{{sample}}/{{preset}}/{{test_id}}/purged.fa",
    sample=SAMPLES,
    preset=PRESETS,
    test_id=PURGE_TESTS,
)

PURGE_HAP_OUTPUTS = expand(
    f"{RESULTS_DIR}/purge_dups/{{sample}}/{{preset}}/{{test_id}}/hap.fa",
    sample=SAMPLES,
    preset=PRESETS,
    test_id=PURGE_TESTS,
)

PURGE_AUTO_CUTOFFS = expand(
    f"{RESULTS_DIR}/purge_dups/{{sample}}/coverage/cutoffs.auto",
    sample=SAMPLES,
)


rule all:
    input:
        PURGE_DUP_OUTPUTS,
        PURGE_HAP_OUTPUTS,
        PURGE_AUTO_CUTOFFS,


rule map_hifi_chunk_to_assembly_paf:
    input:
        assembly=lambda wildcards: FASTA_BY_SAMPLE[wildcards.sample],
        hifi=lambda wildcards: HIFI_READ_BY_CHUNK[wildcards.chunk],
    output:
        paf=f"{RESULTS_DIR}/read_mapping/hifi_paf/{{sample}}/chunks/{{chunk}}.paf.gz",
    log:
        f"{RESULTS_DIR}/logs/read_mapping/minimap2/hifi_paf/{{sample}}/{{chunk}}.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/read_mapping/hifi_paf/{{sample}}/{{chunk}}.tsv"
    conda:
        "envs/minimap2.yaml"
    threads: config["minimap2"]["threads"]
    resources:
        mem_gb=config["minimap2"]["mem_gb"],
    params:
        pigz_threads=config["minimap2"].get("pigz_threads", 4),
        index_size=config["minimap2"].get("index_size", "4G"),
        split_prefix=lambda wildcards: (
            f"{RESULTS_DIR}/read_mapping/hifi_paf/{wildcards.sample}/chunks/"
            f"{wildcards.chunk}.split_idx"
        ),
    shell:
        r"""
        mkdir -p $(dirname {output.paf})
        mkdir -p $(dirname {log})

        minimap2 \
            -t {threads} \
            -I {params.index_size} \
            --split-prefix {params.split_prefix} \
            -x map-hifi \
            {input.assembly} \
            {input.hifi} \
            2>{log} \
            | pigz -p {params.pigz_threads} -c >{output.paf}
        """


rule merge_hifi_to_assembly_paf:
    input:
        lambda wildcards: expand(
            f"{RESULTS_DIR}/read_mapping/hifi_paf/{wildcards.sample}/chunks/{{chunk}}.paf.gz",
            chunk=HIFI_CHUNKS,
        ),
    output:
        paf=f"{RESULTS_DIR}/read_mapping/hifi_paf/{{sample}}/{{sample}}.hifi_to_assembly.paf.gz",
    log:
        f"{RESULTS_DIR}/logs/read_mapping/minimap2/hifi_paf/{{sample}}/merge.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/read_mapping/hifi_paf/{{sample}}.merge.tsv"
    shell:
        r"""
        mkdir -p $(dirname {output.paf})
        mkdir -p $(dirname {log})

        cat {input} >{output.paf} 2>{log}
        """


rule purge_dups_pbcstat:
    input:
        paf=f"{RESULTS_DIR}/read_mapping/hifi_paf/{{sample}}/{{sample}}.hifi_to_assembly.paf.gz",
    output:
        stat=f"{RESULTS_DIR}/purge_dups/{{sample}}/coverage/PB.stat",
        cov=f"{RESULTS_DIR}/purge_dups/{{sample}}/coverage/PB.base.cov",
    log:
        f"{RESULTS_DIR}/logs/purge_dups/pbcstat/{{sample}}.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/purge_dups/pbcstat/{{sample}}.tsv"
    conda:
        "envs/purge_dups.yaml"
    resources:
        mem_gb=config["purge_dups"]["pbcstat_mem_gb"],
        tmpdir=config["purge_dups"].get("tmpdir", "/home/lehoanglinh/tmp"),
    params:
        extra=config["purge_dups"].get("pbcstat_extra", ""),
    wrapper:
        f"{WRAPPER}/bio/purge_dups/pbcstat"


rule purge_dups_split_assembly_fasta:
    input:
        fasta=lambda wildcards: FASTA_BY_SAMPLE[wildcards.sample],
    output:
        split=f"{RESULTS_DIR}/purge_dups/{{sample}}/self_alignment/{{sample}}.split.fa",
    log:
        f"{RESULTS_DIR}/logs/purge_dups/split_fa/{{sample}}.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/purge_dups/split_fa/{{sample}}.tsv"
    conda:
        "envs/purge_dups.yaml"
    threads: 1
    params:
        extra=config["purge_dups"].get("split_fa_extra", ""),
    wrapper:
        f"{WRAPPER}/bio/purge_dups/split_fa"


rule purge_dups_self_align_assembly:
    input:
        split=f"{RESULTS_DIR}/purge_dups/{{sample}}/self_alignment/{{sample}}.split.fa",
    output:
        paf=f"{RESULTS_DIR}/purge_dups/{{sample}}/self_alignment/{{sample}}.{ALIGNER}.{{preset}}.self.paf.gz",
    log:
        f"{RESULTS_DIR}/logs/purge_dups/{ALIGNER}/self_alignment/{{sample}}.{{preset}}.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/purge_dups/self_alignment/{{sample}}.{{preset}}.tsv"
    conda:
        "envs/mm2plus.yaml"
    threads: config["purge_dups"]["self_align_threads"]
    resources:
        mem_gb=config["purge_dups"]["self_align_mem_gb"],
    params:
        aligner=ALIGNER,
        extra=lambda wildcards: config["purge_dups"].get("self_align_extra", ""),
        pigz_threads=config["purge_dups"].get("pigz_threads", 4),
        index_size=config["purge_dups"].get("self_align_index_size", "4G"),
        split_prefix=lambda wildcards: (
            f"{RESULTS_DIR}/purge_dups/{wildcards.sample}/self_alignment/"
            f"{wildcards.sample}.{wildcards.preset}.split_idx"
        ),
    shell:
        r"""
        mkdir -p $(dirname {output.paf})
        mkdir -p $(dirname {log})

        {params.aligner} \
            -t {threads} \
            -I {params.index_size} \
            --split-prefix {params.split_prefix} \
            -x {wildcards.preset} \
            -DP \
            {params.extra} \
            {input.split} \
            {input.split} \
            2>{log} \
            | pigz -p {params.pigz_threads} -c >{output.paf}

        test -s {output.paf}
        """


rule purge_dups_auto_cutoffs:
    input:
        stat=f"{RESULTS_DIR}/purge_dups/{{sample}}/coverage/PB.stat",
    output:
        cutoffs=f"{RESULTS_DIR}/purge_dups/{{sample}}/coverage/cutoffs.auto",
    log:
        f"{RESULTS_DIR}/logs/purge_dups/calcuts/{{sample}}.log",
    conda:
        "envs/purge_dups.yaml"
    threads: 1
    params:
        extra=config["purge_dups"].get("calcuts_extra", ""),
    wrapper:
        f"{WRAPPER}/bio/purge_dups/calcuts"


rule purge_dups_write_manual_cutoffs:
    input:
        stat=f"{RESULTS_DIR}/purge_dups/{{sample}}/coverage/PB.stat",
    output:
        cutoffs=f"{RESULTS_DIR}/purge_dups/{{sample}}/{{preset}}/{{test_id}}/cutoffs",
    log:
        f"{RESULTS_DIR}/logs/purge_dups/calcuts/{{sample}}/{{preset}}/{{test_id}}.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/purge_dups/calcuts/{{sample}}/{{preset}}/{{test_id}}.tsv"
    conda:
        "envs/purge_dups.yaml"
    threads: 1
    params:
        extra=lambda wildcards: build_calcuts_extra(config, wildcards.test_id),
    wrapper:
        f"{WRAPPER}/bio/purge_dups/calcuts"


rule purge_dups_call_dups:
    input:
        self_paf=f"{RESULTS_DIR}/purge_dups/{{sample}}/self_alignment/{{sample}}.{ALIGNER}.{{preset}}.self.paf.gz",
        cov=f"{RESULTS_DIR}/purge_dups/{{sample}}/coverage/PB.base.cov",
        cutoffs=f"{RESULTS_DIR}/purge_dups/{{sample}}/{{preset}}/{{test_id}}/cutoffs",
    output:
        bed=f"{RESULTS_DIR}/purge_dups/{{sample}}/{{preset}}/{{test_id}}/dups.bed",
    log:
        f"{RESULTS_DIR}/logs/purge_dups/purge_dups/{{sample}}/{{preset}}/{{test_id}}.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/purge_dups/purge_dups/{{sample}}/{{preset}}/{{test_id}}.tsv"
    conda:
        "envs/purge_dups.yaml"
    threads: config["purge_dups"].get("purge_threads", 1)
    resources:
        mem_gb=config["purge_dups"]["purge_mem_gb"],
    params:
        extra=config["purge_dups"].get("purge_extra", "-2"),
    shell:
        r"""
        mkdir -p $(dirname {output.bed})
        mkdir -p $(dirname {log})

        awk 'NR==1 && NF==6 {{ok=1}} END {{exit !ok}}' {input.cutoffs} || {{
                                                                            echo "[ERROR] malformed cutoffs (expected 6 tab-separated fields): {input.cutoffs}" >&2
                                                                            cat {input.cutoffs} >&2
                                                                            exit 1
                                                                        }}

        purge_dups \
            {params.extra} \
            -T {input.cutoffs} \
            -c {input.cov} \
            {input.self_paf} \
            >{output.bed} 2>{log}
        """


rule purge_dups_get_seqs:
    input:
        bed=f"{RESULTS_DIR}/purge_dups/{{sample}}/{{preset}}/{{test_id}}/dups.bed",
        fasta=lambda wildcards: FASTA_BY_SAMPLE[wildcards.sample],
    output:
        purged=f"{RESULTS_DIR}/purge_dups/{{sample}}/{{preset}}/{{test_id}}/purged.fa",
        hap=f"{RESULTS_DIR}/purge_dups/{{sample}}/{{preset}}/{{test_id}}/hap.fa",
    log:
        f"{RESULTS_DIR}/logs/purge_dups/get_seqs/{{sample}}/{{preset}}/{{test_id}}.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/purge_dups/get_seqs/{{sample}}/{{preset}}/{{test_id}}.tsv"
    conda:
        "envs/purge_dups.yaml"
    params:
        mode=lambda wildcards: config["purge_dups_tests"][wildcards.test_id].get(
            "get_seqs_mode", config["purge_dups"].get("get_seqs_mode", "-e")
        ),
    shell:
        r"""
        mkdir -p $(dirname {output.purged})
        mkdir -p $(dirname {log})

        WORKDIR=$(realpath $(dirname {output.purged}))
        FASTA=$(realpath {input.fasta})
        BED=$(realpath {input.bed})
        LOG=$(realpath {log})

        cd "$WORKDIR"

        rm -f purged.fa hap.fa

        get_seqs \
            {params.mode} \
            "$BED" \
            "$FASTA" \
            >"$LOG" 2>&1

        test -s purged.fa
        [ -e hap.fa ] || touch hap.fa
        """
