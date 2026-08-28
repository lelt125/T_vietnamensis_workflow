from pathlib import Path
from snakemake.shell import shell


include: "rules/common_config_files.smk"


RUN_NAME, RESULTS_BASE, RESULTS_DIR = get_results_dir(config)
ASSEMBLIES, SAMPLES, FASTA_BY_SAMPLE = load_assemblies(config)
HIFI_READS_DIR, HIFI_READS, HIFI_READ_BY_CHUNK, HIFI_CHUNKS = load_hifi_reads(config)

VALID_BUSCO_LINEAGES, DEFAULT_BUSCO_LINEAGE = load_valid_busco_lineages(config)
BUSCO_LINEAGE_BY_SAMPLE = build_busco_lineage_by_sample(
    ASSEMBLIES, VALID_BUSCO_LINEAGES, DEFAULT_BUSCO_LINEAGE
)

SEQKIT_SAMPLES = get_qc_samples(ASSEMBLIES, "run_seqkit", default=True)
BUSCO_SAMPLES = get_qc_samples(ASSEMBLIES, "run_busco", default=False)
MERQURY_SAMPLES = get_qc_samples(ASSEMBLIES, "run_merqury", default=False)
MOSDEPTH_SAMPLES = get_qc_samples(ASSEMBLIES, "run_mosdepth", default=False)

SEQKIT_STATS_OUTPUTS = expand(
    f"{RESULTS_DIR}/genome_qc/seqkit/{{sample}}.seqkit_stats.tsv", sample=SEQKIT_SAMPLES
)

SEQKIT_FX2TAB_OUTPUTS = expand(
    f"{RESULTS_DIR}/genome_qc/seqkit/{{sample}}.seqkit_fx2tab.tsv",
    sample=SEQKIT_SAMPLES,
)

BUSCO_OUTPUTS = [
    f"{RESULTS_DIR}/genome_qc/busco/{sample}/{BUSCO_LINEAGE_BY_SAMPLE[sample]}"
    for sample in BUSCO_SAMPLES
]

MERQURY_QV_OUTPUTS = expand(
    f"{RESULTS_DIR}/genome_qc/merqury/{{sample}}/{{sample}}.qv", sample=MERQURY_SAMPLES
)

MERQURY_COMPLETENESS_OUTPUTS = expand(
    f"{RESULTS_DIR}/genome_qc/merqury/{{sample}}/{{sample}}.completeness.stats",
    sample=MERQURY_SAMPLES,
)

MERQURY_SPECTRA_CN_OUTPUTS = expand(
    f"{RESULTS_DIR}/genome_qc/merqury/{{sample}}/{{sample}}.asm.spectra-cn.ln.png",
    sample=MERQURY_SAMPLES,
)

MOSDEPTH_SUMMARY_OUTPUTS = expand(
    f"{RESULTS_DIR}/genome_qc/mosdepth/{{sample}}/{{sample}}.hifi.mosdepth.summary.txt",
    sample=MOSDEPTH_SAMPLES,
)

MOSDEPTH_GLOBAL_DIST_OUTPUTS = expand(
    f"{RESULTS_DIR}/genome_qc/mosdepth/{{sample}}/{{sample}}.hifi.mosdepth.global.dist.txt",
    sample=MOSDEPTH_SAMPLES,
)

MERYL_K = config["meryl"]["kmer"]

MERYL_CHUNK_DB = f"{RESULTS_DIR}/genome_qc/meryl/chunks/{{chunk}}.k{MERYL_K}.meryl"
MERYL_CHUNK_STAT = (
    f"{RESULTS_DIR}/genome_qc/meryl/chunks/{{chunk}}.k{MERYL_K}.statistics"
)

MERYL_READS_DB = f"{RESULTS_DIR}/genome_qc/meryl/hifi.combined.k{MERYL_K}.meryl"
MERYL_READS_STAT = f"{RESULTS_DIR}/genome_qc/meryl/hifi.combined.k{MERYL_K}.statistics"
MERYL_READS_VALIDATED = (
    f"{RESULTS_DIR}/genome_qc/meryl/hifi.combined.k{MERYL_K}.validated"
)

MERYL_ASM_DB = f"{RESULTS_DIR}/genome_qc/meryl/assemblies/{{sample}}.k{MERYL_K}.meryl"
MERYL_ASM_STAT = (
    f"{RESULTS_DIR}/genome_qc/meryl/assemblies/{{sample}}.k{MERYL_K}.statistics"
)
MERYL_ASM_VALIDATED = (
    f"{RESULTS_DIR}/genome_qc/meryl/assemblies/{{sample}}.k{MERYL_K}.validated"
)

PROJECT_ROOT = Path(workflow.basedir).parent
ASSERT_MERYL_SCRIPT = str(PROJECT_ROOT / "scripts" / "assert_meryl_kmers.py")


wildcard_constraints:
    sample=r"[^/]+",
    chunk=r"[^/]+",


rule all:
    input:
        SEQKIT_STATS_OUTPUTS,
        SEQKIT_FX2TAB_OUTPUTS,
        BUSCO_OUTPUTS,
        MERQURY_QV_OUTPUTS,
        MERQURY_COMPLETENESS_OUTPUTS,
        MERQURY_SPECTRA_CN_OUTPUTS,
        MOSDEPTH_SUMMARY_OUTPUTS,
        MOSDEPTH_GLOBAL_DIST_OUTPUTS,


rule seqkit_assembly_stats:
    input:
        fasta=lambda wildcards: FASTA_BY_SAMPLE[wildcards.sample],
    output:
        stats=f"{RESULTS_DIR}/genome_qc/seqkit/{{sample}}.seqkit_stats.tsv",
        fx2tab=f"{RESULTS_DIR}/genome_qc/seqkit/{{sample}}.seqkit_fx2tab.tsv",
    log:
        stats=f"{RESULTS_DIR}/logs/genome_qc/seqkit/{{sample}}.seqkit_stats.log",
        fx2tab=f"{RESULTS_DIR}/logs/genome_qc/seqkit/{{sample}}.seqkit_fx2tab.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/genome_qc/seqkit/{{sample}}.tsv"
    threads: config["seqkit"]["threads"]
    shell:
        r"""
        seqkit stats -j {threads} -a -T {input.fasta} \
            >{output.stats} \
            2>{log.stats}

        seqkit fx2tab -j {threads} \
            --name \
            --length \
            --gc \
            {input.fasta} \
            >{output.fx2tab} \
            2>{log.fx2tab}
        """


rule prepare_busco_fasta:
    input:
        fasta=lambda wildcards: FASTA_BY_SAMPLE[wildcards.sample],
    output:
        fasta=f"{RESULTS_DIR}/genome_qc/busco_input/{{sample}}.busco_input.fa",
    log:
        f"{RESULTS_DIR}/logs/genome_qc/busco/prepare_fasta/{{sample}}.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/genome_qc/busco/prepare_fasta/{{sample}}.tsv"
    threads: 8
    shell:
        r"""
        seqkit seq -j {threads} -w 80 {input.fasta} \
            >{output.fasta} \
            2>{log}
        """


rule busco_genome:
    input:
        fasta=f"{RESULTS_DIR}/genome_qc/busco_input/{{sample}}.busco_input.fa",
    output:
        directory(f"{RESULTS_DIR}/genome_qc/busco/{{sample}}/{{lineage}}"),
    log:
        f"{RESULTS_DIR}/logs/genome_qc/busco/{{sample}}/{{lineage}}.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/genome_qc/busco/{{sample}}/{{lineage}}.tsv"
    threads: config["busco"]["threads"]
    resources:
        mem_gb=config["busco"]["mem_gb"],
    params:
        out_parent=lambda wildcards: f"{RESULTS_DIR}/genome_qc/busco/{wildcards.sample}",
        download_path=lambda wildcards: str(
            Path(config["busco"]["download_path"]).resolve()
        ),
        extra=config["busco"].get("extra", ""),
    shell:
        r"""
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
            >{log} 2>&1
        """


# ---------------------------------------------------------------------------
# Meryl / Merqury
#
# Counting rules only count and emit `meryl statistics` output. Validation
# (comparing `present` against an independently-derived expected count) is
# a separate rule per database, producing a `.validated` marker. Downstream
# rules depend on the marker, not just the database, so Snakemake cannot
# schedule merqury on an unvalidated database.
# ---------------------------------------------------------------------------


rule meryl_count_hifi_chunk:
    input:
        hifi=lambda wildcards: HIFI_READ_BY_CHUNK[wildcards.chunk],
    output:
        meryl=directory(MERYL_CHUNK_DB),
        stats=MERYL_CHUNK_STAT,
    log:
        meryl_count=f"{RESULTS_DIR}/logs/genome_qc/meryl/chunks/{{chunk}}.meryl_count.log",
        stats=f"{RESULTS_DIR}/logs/genome_qc/meryl/chunks/{{chunk}}.meryl_statistics.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/genome_qc/meryl/chunks/{{chunk}}.tsv"
    threads: config["meryl"]["count_threads"]
    resources:
        mem_gb=config["meryl"]["count_mem_gb"],
    params:
        kmer=MERYL_K,
    shell:
        r"""
        rm -rf {output.meryl}

        meryl k={params.kmer} count \
            threads={threads} \
            memory={resources.mem_gb}g \
            {input.hifi} \
            output {output.meryl} \
            >{log.meryl_count} 2>&1

        meryl statistics {output.meryl} \
            >{output.stats} \
            2>{log.stats}
        """


rule meryl_merge_hifi_chunks:
    input:
        meryl=expand(MERYL_CHUNK_DB, chunk=HIFI_CHUNKS),
    output:
        meryl=directory(MERYL_READS_DB),
        stats=MERYL_READS_STAT,
    log:
        merge=f"{RESULTS_DIR}/logs/genome_qc/meryl/hifi.combined.meryl_merge.log",
        stats=f"{RESULTS_DIR}/logs/genome_qc/meryl/hifi.combined.meryl_statistics.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/genome_qc/meryl/hifi.combined.merge.tsv"
    threads: config["meryl"]["merge_threads"]
    resources:
        mem_gb=config["meryl"]["merge_mem_gb"],
    shell:
        r"""
        rm -rf {output.meryl}

        meryl union-sum \
            threads={threads} \
            memory={resources.mem_gb}g \
            {input.meryl} \
            output {output.meryl} \
            >{log.merge} 2>&1

        meryl statistics {output.meryl} \
            >{output.stats} \
            2>{log.stats}
        """


rule assert_meryl_kmers_combined:
    input:
        stats=MERYL_READS_STAT,
        parts=expand(MERYL_CHUNK_STAT, chunk=HIFI_CHUNKS),
    output:
        MERYL_READS_VALIDATED,
    log:
        f"{RESULTS_DIR}/logs/genome_qc/meryl/hifi.combined.assert_kmers.log",
    params:
        mode="sum",
        tolerance=0.0,
        label="hifi.combined",
    script:
        ASSERT_MERYL_SCRIPT


rule prepare_merqury_fasta:
    input:
        fasta=lambda wildcards: FASTA_BY_SAMPLE[wildcards.sample],
    output:
        fasta=temp(f"{RESULTS_DIR}/genome_qc/merqury_input/{{sample}}.fasta"),
    log:
        f"{RESULTS_DIR}/logs/genome_qc/merqury/prepare_fasta/{{sample}}.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/genome_qc/merqury/prepare_fasta/{{sample}}.tsv"
    threads: 8
    shell:
        r"""
        seqkit seq -j {threads} {input.fasta} \
            >{output.fasta} \
            2>{log}

        test -s {output.fasta}
        grep -q "^>" {output.fasta}
        """


rule meryl_count_assembly:
    input:
        fasta=f"{RESULTS_DIR}/genome_qc/merqury_input/{{sample}}.fasta",
    output:
        meryl=directory(MERYL_ASM_DB),
        stats=MERYL_ASM_STAT,
    log:
        meryl_count=f"{RESULTS_DIR}/logs/genome_qc/meryl/assemblies/{{sample}}.meryl_count.log",
        stats=f"{RESULTS_DIR}/logs/genome_qc/meryl/assemblies/{{sample}}.meryl_statistics.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/genome_qc/meryl/assemblies/{{sample}}.tsv"
    threads: config["meryl"]["asm_count_threads"]
    resources:
        mem_gb=config["meryl"]["asm_count_mem_gb"],
    params:
        kmer=MERYL_K,
    shell:
        r"""
        rm -rf {output.meryl}

        meryl k={params.kmer} count \
            threads={threads} \
            memory={resources.mem_gb}g \
            {input.fasta} \
            output {output.meryl} \
            >{log.meryl_count} 2>&1

        meryl statistics {output.meryl} \
            >{output.stats} \
            2>{log.stats}
        """


rule assert_meryl_kmers_assembly:
    input:
        stats=MERYL_ASM_STAT,
        seqkit=f"{RESULTS_DIR}/genome_qc/seqkit/{{sample}}.seqkit_stats.tsv",
    output:
        MERYL_ASM_VALIDATED,
    log:
        f"{RESULTS_DIR}/logs/genome_qc/meryl/assemblies/{{sample}}.assert_kmers.log",
    params:
        mode="seqkit",
        kmer=MERYL_K,
        tolerance=config["meryl"].get("asm_tolerance", 0.005),
        label=lambda wildcards: wildcards.sample,
    script:
        ASSERT_MERYL_SCRIPT


rule merqury_assembly_qc:
    input:
        fasta=f"{RESULTS_DIR}/genome_qc/merqury_input/{{sample}}.fasta",
        reads_meryl=MERYL_READS_DB,
        reads_validated=MERYL_READS_VALIDATED,
        asm_meryl=MERYL_ASM_DB,
        asm_validated=MERYL_ASM_VALIDATED,
    output:
        qv=f"{RESULTS_DIR}/genome_qc/merqury/{{sample}}/{{sample}}.qv",
        completeness=f"{RESULTS_DIR}/genome_qc/merqury/{{sample}}/{{sample}}.completeness.stats",
        asm_qv=f"{RESULTS_DIR}/genome_qc/merqury/{{sample}}/{{sample}}.asm.qv",
        spectra_cn_hist=f"{RESULTS_DIR}/genome_qc/merqury/{{sample}}/{{sample}}.asm.spectra-cn.hist",
        spectra_cn_png=f"{RESULTS_DIR}/genome_qc/merqury/{{sample}}/{{sample}}.asm.spectra-cn.ln.png",
        spectra_asm_hist=f"{RESULTS_DIR}/genome_qc/merqury/{{sample}}/{{sample}}.spectra-asm.hist",
        spectra_asm_png=f"{RESULTS_DIR}/genome_qc/merqury/{{sample}}/{{sample}}.spectra-asm.ln.png",
    log:
        f"{RESULTS_DIR}/logs/genome_qc/merqury/{{sample}}.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/genome_qc/merqury/{{sample}}.tsv"
    threads: config["merqury"]["threads"]
    resources:
        mem_gb=config["merqury"]["mem_gb"],
    params:
        outdir=lambda wildcards: f"{RESULTS_DIR}/genome_qc/merqury/{wildcards.sample}",
        prefix=lambda wildcards: wildcards.sample,
        java_heap=config["merqury"].get("java_heap", "64g"),
        kmer=MERYL_K,
    shell:
        r"""
        mkdir -p {params.outdir}

        READS_MERYL=$(realpath {input.reads_meryl})
        ASM_MERYL=$(realpath {input.asm_meryl})
        FASTA=$(realpath {input.fasta})
        LOG=$(realpath {log})

        export _JAVA_OPTIONS="-Xmx{params.java_heap}"

        cd {params.outdir}
        mkdir -p logs

        # Remove intermediates from any earlier run: derivatives of a
        # previously broken read database would otherwise be silently reused.
        rm -rf read.k{params.kmer}.asm.*.meryl asm.0.meryl asm.solid.meryl asm.meryl
        rm -f asm.fasta *.filt *.hist.ploidy

        # merqury names the assembly database from the FASTA basename and
        # skips counting when it already exists - this is what makes it
        # reuse our pre-built, validated assembly database instead of
        # recounting (and re-risking the silent batch bug) itself.
        ln -s "$FASTA" asm.fasta
        ln -s "$ASM_MERYL" asm.meryl

        merqury.sh "$READS_MERYL" asm.fasta {params.prefix} >"$LOG" 2>&1

        if grep -qE "Counting [0-9]+ \(estimated\) (billion|million)" "$LOG"; then
            echo "FATAL: merqury re-counted the assembly instead of reusing the pre-built, validated database. Inspect $LOG." >&2
            exit 1
        fi

        for f in \
            {params.prefix}.qv \
            {params.prefix}.completeness.stats \
            {params.prefix}.asm.qv \
            {params.prefix}.asm.spectra-cn.hist \
            {params.prefix}.asm.spectra-cn.ln.png \
            {params.prefix}.spectra-asm.hist \
            {params.prefix}.spectra-asm.ln.png; do
            test -s "$f" || {{ echo "FATAL: missing or empty merqury output: $f" >&2; exit 1; }}
        done
        """


rule map_hifi_chunk_to_assembly_bam:
    input:
        assembly=lambda wildcards: FASTA_BY_SAMPLE[wildcards.sample],
        hifi=lambda wildcards: HIFI_READ_BY_CHUNK[wildcards.chunk],
    output:
        bam=temp(
            f"{RESULTS_DIR}/read_mapping/hifi_bam/{{sample}}/chunks/{{chunk}}.sorted.bam"
        ),
    log:
        minimap2=f"{RESULTS_DIR}/logs/read_mapping/minimap2/hifi_bam/{{sample}}/{{chunk}}.minimap2.log",
        sort=f"{RESULTS_DIR}/logs/read_mapping/samtools_sort/hifi_bam/{{sample}}/{{chunk}}.sort.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/read_mapping/hifi_bam/{{sample}}/{{chunk}}.tsv"
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
        split_prefix=lambda wildcards: f"{RESULTS_DIR}/read_mapping/hifi_bam/{wildcards.sample}/chunks/{wildcards.chunk}.split_idx",
    shell:
        r"""
        minimap2 \
            -t {params.minimap2_threads} \
            -I {params.index_size} \
            --split-prefix {params.split_prefix} \
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


rule merge_hifi_to_assembly_bam:
    input:
        bams=lambda wildcards: expand(
            f"{RESULTS_DIR}/read_mapping/hifi_bam/{wildcards.sample}/chunks/{{chunk}}.sorted.bam",
            chunk=HIFI_CHUNKS,
        ),
    output:
        bam=f"{RESULTS_DIR}/read_mapping/hifi_bam/{{sample}}/{{sample}}.hifi_to_assembly.sorted.bam",
        bai=f"{RESULTS_DIR}/read_mapping/hifi_bam/{{sample}}/{{sample}}.hifi_to_assembly.sorted.bam.bai",
    log:
        merge=f"{RESULTS_DIR}/logs/read_mapping/samtools_merge/hifi_bam/{{sample}}.merge.log",
        index=f"{RESULTS_DIR}/logs/read_mapping/samtools_index/hifi_bam/{{sample}}.index.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/read_mapping/hifi_bam/{{sample}}.merge.tsv"
    threads: config["minimap2"].get("merge_threads", 8)
    resources:
        mem_gb=config["minimap2"].get("merge_mem_gb", 32),
    params:
        bam_list=f"{RESULTS_DIR}/read_mapping/hifi_bam/{{sample}}/{{sample}}.bam_list.txt",
    shell:
        r"""
        printf "%s\n" {input.bams} >{params.bam_list}

        samtools merge \
            -@ {threads} \
            -f \
            -b {params.bam_list} \
            {output.bam} \
            >{log.merge} 2>&1

        samtools index \
            -@ {threads} \
            {output.bam} \
            >{log.index} 2>&1

        test -s {output.bam}
        test -s {output.bai}
        samtools quickcheck {output.bam}
        """


rule mosdepth_hifi_coverage:
    input:
        bam=f"{RESULTS_DIR}/read_mapping/hifi_bam/{{sample}}/{{sample}}.hifi_to_assembly.sorted.bam",
        bai=f"{RESULTS_DIR}/read_mapping/hifi_bam/{{sample}}/{{sample}}.hifi_to_assembly.sorted.bam.bai",
    output:
        summary=f"{RESULTS_DIR}/genome_qc/mosdepth/{{sample}}/{{sample}}.hifi.mosdepth.summary.txt",
        global_dist=f"{RESULTS_DIR}/genome_qc/mosdepth/{{sample}}/{{sample}}.hifi.mosdepth.global.dist.txt",
    log:
        f"{RESULTS_DIR}/logs/genome_qc/mosdepth/{{sample}}.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/genome_qc/mosdepth/{{sample}}.tsv"
    threads: config["mosdepth"]["threads"]
    resources:
        mem_gb=config["mosdepth"]["mem_gb"],
    params:
        outdir=lambda wildcards: f"{RESULTS_DIR}/genome_qc/mosdepth/{wildcards.sample}",
        prefix=lambda wildcards: f"{RESULTS_DIR}/genome_qc/mosdepth/{wildcards.sample}/{wildcards.sample}.hifi",
        extra=config["mosdepth"].get("extra", "--fast-mode --no-per-base"),
    shell:
        r"""
        mosdepth \
            -t {threads} \
            {params.extra} \
            {params.prefix} \
            {input.bam} \
            >{log} 2>&1

        test -s {output.summary}
        test -s {output.global_dist}
        """
