include: "rules/common_config_files.smk"

RUN_NAME, RESULTS_BASE, RESULTS_DIR = get_results_dir(config)
ASSEMBLIES, SAMPLES, FASTA_BY_SAMPLE = load_assemblies(config)
HIFI_READS_DIR, HIFI_READS, HIFI_READ_BY_CHUNK, HIFI_CHUNKS = load_hifi_reads(config)

PURGE_TESTS = list(config["purge_dups_tests"].keys())


PURGE_DUP_OUTPUTS = expand(
    f"{RESULTS_DIR}/purge_dups/{{sample}}/{{test_id}}/purged.fa",
    sample=SAMPLES,
    test_id=PURGE_TESTS
)

PURGE_HAP_OUTPUTS = expand(
    f"{RESULTS_DIR}/purge_dups/{{sample}}/{{test_id}}/hap.fa",
    sample=SAMPLES,
    test_id=PURGE_TESTS
)


rule all:
    input:
        PURGE_DUP_OUTPUTS,
        PURGE_HAP_OUTPUTS


rule map_hifi_chunk_to_assembly_paf:
    input:
        assembly=lambda wildcards: FASTA_BY_SAMPLE[wildcards.sample],
        hifi=lambda wildcards: HIFI_READ_BY_CHUNK[wildcards.chunk]
    output:
        paf=f"{RESULTS_DIR}/read_mapping/hifi_paf/{{sample}}/chunks/{{chunk}}.paf.gz"
    threads:
        config["minimap2_read_mapping"]["map_threads"]
    resources:
        mem_gb=config["minimap2_read_mapping"]["map_mem_gb"]
    params:
        pigz_threads=config["minimap2_read_mapping"].get("pigz_threads", 4)
    log:
        f"{RESULTS_DIR}/logs/read_mapping/minimap2/hifi_paf/{{sample}}/{{chunk}}.log"
    shell:
        r"""
        mkdir -p $(dirname {output.paf})
        mkdir -p $(dirname {log})

        minimap2 \
            -t {threads} \
            -x map-hifi \
            {input.assembly} \
            {input.hifi} \
            2> {log} \
        | pigz -p {params.pigz_threads} -c > {output.paf}
        """


rule merge_hifi_to_assembly_paf:
    input:
        lambda wildcards: expand(
            f"{RESULTS_DIR}/read_mapping/hifi_paf/{wildcards.sample}/chunks/{{chunk}}.paf.gz",
            chunk=HIFI_CHUNKS
        )
    output:
        paf=f"{RESULTS_DIR}/read_mapping/hifi_paf/{{sample}}/{{sample}}.hifi_to_assembly.paf.gz"
    log:
        f"{RESULTS_DIR}/logs/read_mapping/minimap2/hifi_paf/{{sample}}/merge.log"
    shell:
        r"""
        mkdir -p $(dirname {output.paf})
        mkdir -p $(dirname {log})

        cat {input} > {output.paf} 2> {log}
        """


rule purge_dups_pbcstat:
    input:
        paf=f"{RESULTS_DIR}/read_mapping/hifi_paf/{{sample}}/{{sample}}.hifi_to_assembly.paf.gz"
    output:
        stat=f"{RESULTS_DIR}/purge_dups/{{sample}}/coverage/PB.stat",
        cov=f"{RESULTS_DIR}/purge_dups/{{sample}}/coverage/PB.base.cov"
    resources:
        mem_gb=config["purge_dups"]["pbcstat_mem_gb"]
    params:
        outdir=lambda wildcards: f"{RESULTS_DIR}/purge_dups/{wildcards.sample}/coverage",
        extra=config["purge_dups"].get("pbcstat_extra", "")
    log:
        f"{RESULTS_DIR}/logs/purge_dups/pbcstat/{{sample}}.log"
    shell:
        r"""
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})

        pbcstat \
            -O {params.outdir} \
            {params.extra} \
            {input.paf} \
            > {log} 2>&1

        test -s {output.stat}
        test -s {output.cov}
        """


rule purge_dups_auto_cutoffs:
    input:
        stat=f"{RESULTS_DIR}/purge_dups/{{sample}}/coverage/PB.stat"
    output:
        cutoffs=f"{RESULTS_DIR}/purge_dups/{{sample}}/coverage/cutoffs.auto"
    log:
        f"{RESULTS_DIR}/logs/purge_dups/calcuts/{{sample}}.log"
    shell:
        r"""
        mkdir -p $(dirname {output.cutoffs})
        mkdir -p $(dirname {log})

        calcuts {input.stat} > {output.cutoffs} 2> {log}
        """


rule purge_dups_write_manual_cutoffs:
    output:
        cutoffs=f"{RESULTS_DIR}/purge_dups/{{sample}}/{{test_id}}/cutoffs"
    params:
        low=lambda wildcards: config["purge_dups_tests"][wildcards.test_id]["low"],
        mid=lambda wildcards: config["purge_dups_tests"][wildcards.test_id]["mid"],
        high=lambda wildcards: config["purge_dups_tests"][wildcards.test_id]["high"]
    shell:
        r"""
        mkdir -p $(dirname {output.cutoffs})
        echo "{params.low} {params.mid} {params.high}" > {output.cutoffs}
        """


rule purge_dups_split_assembly_fasta:
    input:
        fasta=lambda wildcards: FASTA_BY_SAMPLE[wildcards.sample]
    output:
        split=f"{RESULTS_DIR}/purge_dups/{{sample}}/self_alignment/{{sample}}.split.fa"
    log:
        f"{RESULTS_DIR}/logs/purge_dups/split_fa/{{sample}}.log"
    shell:
        r"""
        mkdir -p $(dirname {output.split})
        mkdir -p $(dirname {log})

        split_fa {input.fasta} > {output.split} 2> {log}
        """


rule purge_dups_self_align_assembly:
    input:
        split=f"{RESULTS_DIR}/purge_dups/{{sample}}/self_alignment/{{sample}}.split.fa"
    output:
        paf=f"{RESULTS_DIR}/purge_dups/{{sample}}/self_alignment/{{sample}}.{config['purge_dups'].get('self_align_aligner', 'minimap2')}.{config['purge_dups'].get('self_align_preset', 'asm5')}.self.paf.gz"
    threads:
        config["purge_dups"]["self_align_threads"]
    resources:
        mem_gb=config["purge_dups"]["self_align_mem_gb"]
    params:
        aligner=config["purge_dups"].get("self_align_aligner", "minimap2"),
        preset=config["purge_dups"].get("self_align_preset", "asm5"),
        extra=config["purge_dups"].get("self_align_extra", ""),
        pigz_threads=config["purge_dups"].get("pigz_threads", 4)
    log:
        f"{RESULTS_DIR}/logs/purge_dups/{config['purge_dups'].get('self_align_aligner', 'minimap2')}/self_alignment/{{sample}}.{config['purge_dups'].get('self_align_preset', 'asm5')}.log"
    conda:
        "envs/mm2plus.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.paf})
        mkdir -p $(dirname {log})

        {params.aligner} \
            -t {threads} \
            -x {params.preset} \
            -DP \
            {params.extra} \
            {input.split} \
            {input.split} \
            2> {log} \
        | pigz -p {params.pigz_threads} -c > {output.paf}

        test -s {output.paf}
        """


rule purge_dups_call_dups:
    input:
        self_paf=f"{RESULTS_DIR}/purge_dups/{{sample}}/self_alignment/{{sample}}.{config['purge_dups'].get('self_align_aligner', 'minimap2')}.{config['purge_dups'].get('self_align_preset', 'asm5')}.self.paf.gz",
        cov=f"{RESULTS_DIR}/purge_dups/{{sample}}/coverage/PB.base.cov",
        cutoffs=f"{RESULTS_DIR}/purge_dups/{{sample}}/{{test_id}}/cutoffs"
    output:
        bed=f"{RESULTS_DIR}/purge_dups/{{sample}}/{{test_id}}/dups.bed"
    threads:
        config["purge_dups"].get("purge_threads", 1)
    resources:
        mem_gb=config["purge_dups"]["purge_mem_gb"]
    params:
        extra=config["purge_dups"].get("purge_extra", "-2")
    log:
        f"{RESULTS_DIR}/logs/purge_dups/purge_dups/{{sample}}/{{test_id}}.log"
    benchmark:
        f"{RESULTS_DIR}/benchmarks/purge_dups/purge_dups/{{sample}}/{{test_id}}.tsv"
    shell:
        r"""
        mkdir -p $(dirname {output.bed})
        mkdir -p $(dirname {log})

        purge_dups \
            {params.extra} \
            -T {input.cutoffs} \
            -c {input.cov} \
            {input.self_paf} \
            > {output.bed} 2> {log}
        """


rule purge_dups_get_seqs:
    input:
        bed=f"{RESULTS_DIR}/purge_dups/{{sample}}/{{test_id}}/dups.bed",
        fasta=lambda wildcards: FASTA_BY_SAMPLE[wildcards.sample]
    output:
        purged=f"{RESULTS_DIR}/purge_dups/{{sample}}/{{test_id}}/purged.fa",
        hap=f"{RESULTS_DIR}/purge_dups/{{sample}}/{{test_id}}/hap.fa"
    params:
        mode=lambda wildcards: config["purge_dups_tests"][wildcards.test_id].get(
            "get_seqs_mode",
            config["purge_dups"].get("get_seqs_mode", "-e")
        )
    log:
        f"{RESULTS_DIR}/logs/purge_dups/get_seqs/{{sample}}/{{test_id}}.log"
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
            > "$LOG" 2>&1

        test -s purged.fa
        [ -e hap.fa ] || touch hap.fa
        """