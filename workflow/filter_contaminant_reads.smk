include: "rules/common_config_files.smk"


RUN_NAME, RESULTS_BASE, RESULTS_DIR = get_results_dir(config)
ASSEMBLIES, SAMPLES, FASTA_BY_SAMPLE = load_assemblies(config)
HIFI_READS_DIR, HIFI_READS, HIFI_READ_BY_CHUNK, HIFI_CHUNKS = load_hifi_reads(config)

CONTAM = config["contaminant_filter"]

# Assembly the contaminant contig IDs were called against. The IDs must match
# the headers of this FASTA exactly, so it has to be the pre-purge, pre-scaffold
# assembly that FCS-GX screened.
CONTAM_SAMPLE = CONTAM["reference_sample"]
CONTAM_ASSEMBLY = FASTA_BY_SAMPLE[CONTAM_SAMPLE]
CONTAM_IDS = CONTAM["contaminant_ids"]

FILTERED_READS = expand(
    f"{RESULTS_DIR}/contaminant_filter/filtered_reads/{{chunk}}.filtered.fastq.gz",
    chunk=HIFI_CHUNKS,
)

FILTER_REPORT = f"{RESULTS_DIR}/contaminant_filter/contaminant_filter_report.tsv"


rule all:
    input:
        FILTERED_READS,
        FILTER_REPORT,


rule contaminant_reference_fasta:
    """
    Pull the flagged contigs out of the screened assembly. This is the stage-1
    bait: a small reference assembled from the same library as the contaminant
    reads, so true contaminant reads align to it at near-100% identity.
    """
    input:
        assembly=CONTAM_ASSEMBLY,
        ids=CONTAM_IDS,
    output:
        fasta=f"{RESULTS_DIR}/contaminant_filter/reference/contaminants.fa",
    log:
        f"{RESULTS_DIR}/logs/contaminant_filter/reference/contaminants.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/contaminant_filter/reference/contaminants.tsv"
    threads: CONTAM.get("seqkit_threads", 8)
    shell:
        r"""
        mkdir -p $(dirname {output.fasta})
        mkdir -p $(dirname {log})

        seqkit grep -j {threads} -f {input.ids} {input.assembly} \
            >{output.fasta} \
            2>{log}

        test -s {output.fasta}

        N_REQUESTED=$(grep -vc '^#' {input.ids})
        N_FOUND=$(grep -c '^>' {output.fasta})

        echo "requested=$N_REQUESTED found=$N_FOUND" >>{log}

        if [ "$N_REQUESTED" -ne "$N_FOUND" ]; then
            echo "ERROR: $N_REQUESTED IDs requested but $N_FOUND contigs found." >&2
            echo "Contig IDs do not match the assembly headers." >&2
            exit 1
        fi
        """


rule map_hifi_chunk_to_contaminant_ref:
    """
    Stage 1. Map every read against the contaminant contigs alone. The
    reference is ~1 Gb, so this is cheap, and the overwhelming majority of
    reads produce no alignment at all and drop out here.
    """
    input:
        ref=f"{RESULTS_DIR}/contaminant_filter/reference/contaminants.fa",
        hifi=lambda wildcards: HIFI_READ_BY_CHUNK[wildcards.chunk],
    output:
        paf=temp(f"{RESULTS_DIR}/contaminant_filter/stage1_paf/{{chunk}}.paf.gz"),
    log:
        f"{RESULTS_DIR}/logs/contaminant_filter/stage1_minimap2/{{chunk}}.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/contaminant_filter/stage1_minimap2/{{chunk}}.tsv"
    threads: CONTAM.get("stage1_threads", config["minimap2"]["threads"])
    resources:
        mem_gb=CONTAM.get("stage1_mem_gb", 64),
    params:
        index_size=CONTAM.get("stage1_index_size", "8G"),
        preset=CONTAM.get("preset", "map-hifi"),
        pigz_threads=config["minimap2"].get("pigz_threads", 4),
    shell:
        r"""
        mkdir -p $(dirname {output.paf})
        mkdir -p $(dirname {log})

        minimap2 \
            -t {threads} \
            -I {params.index_size} \
            -x {params.preset} \
            --secondary=no \
            {input.ref} \
            {input.hifi} \
            2>{log} \
            | pigz -p {params.pigz_threads} -c >{output.paf}
        """


rule candidate_contaminant_read_ids:
    """
    Stage 1 selection, deliberately permissive: any read with an alignment
    block over min_align_len becomes a candidate. Being loose here costs
    nothing, since stage 2 makes the actual call. Being strict here would
    silently drop contaminant reads before they were ever assessed.
    """
    input:
        paf=f"{RESULTS_DIR}/contaminant_filter/stage1_paf/{{chunk}}.paf.gz",
    output:
        ids=f"{RESULTS_DIR}/contaminant_filter/candidates/{{chunk}}.candidate_ids.txt",
    log:
        f"{RESULTS_DIR}/logs/contaminant_filter/candidates/{{chunk}}.log",
    params:
        min_align_len=CONTAM.get("stage1_min_align_len", 1000),
    shell:
        r"""
        mkdir -p $(dirname {output.ids})
        mkdir -p $(dirname {log})

        zcat {input.paf} \
            | awk -v minlen={params.min_align_len} '$11 >= minlen {{ print $1 }}' \
            | sort -u \
                >{output.ids} \
                2>{log}

        echo "candidates=$(wc -l <{output.ids})" >>{log}
        """


rule extract_candidate_reads:
    input:
        hifi=lambda wildcards: HIFI_READ_BY_CHUNK[wildcards.chunk],
        ids=f"{RESULTS_DIR}/contaminant_filter/candidates/{{chunk}}.candidate_ids.txt",
    output:
        fastq=temp(
            f"{RESULTS_DIR}/contaminant_filter/candidates/{{chunk}}.candidates.fastq.gz"
        ),
    log:
        f"{RESULTS_DIR}/logs/contaminant_filter/extract_candidates/{{chunk}}.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/contaminant_filter/extract_candidates/{{chunk}}.tsv"
    threads: CONTAM.get("seqkit_threads", 8)
    shell:
        r"""
        mkdir -p $(dirname {output.fastq})
        mkdir -p $(dirname {log})

        if [ ! -s {input.ids} ]; then
            echo "No candidates; writing empty FASTQ." >{log}
            : | gzip -c >{output.fastq}
            exit 0
        fi

        seqkit grep -j {threads} -f {input.ids} {input.hifi} \
            | gzip -c >{output.fastq} \
                2>{log}
        """


rule map_candidate_reads_to_full_assembly:
    """
    Stage 2. Re-map only the candidates against the FULL screened assembly,
    contaminant contigs included. This is what makes the comparison
    competitive: a read from a conserved or repetitive region can align to
    both, and minimap2 assigns the primary alignment to whichever it fits
    better. Mapping against the contaminant contigs in isolation and deleting
    every hit is how conserved tetrapod sequence gets thrown away.
    """
    input:
        assembly=CONTAM_ASSEMBLY,
        fastq=f"{RESULTS_DIR}/contaminant_filter/candidates/{{chunk}}.candidates.fastq.gz",
    output:
        paf=temp(f"{RESULTS_DIR}/contaminant_filter/stage2_paf/{{chunk}}.paf.gz"),
    log:
        f"{RESULTS_DIR}/logs/contaminant_filter/stage2_minimap2/{{chunk}}.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/contaminant_filter/stage2_minimap2/{{chunk}}.tsv"
    threads: CONTAM.get("stage2_threads", config["minimap2"]["threads"])
    resources:
        mem_gb=CONTAM.get("stage2_mem_gb", config["minimap2"]["mem_gb"]),
    params:
        index_size=CONTAM.get(
            "stage2_index_size", config["minimap2"].get("index_size", "24G")
        ),
        preset=CONTAM.get("preset", "map-hifi"),
        pigz_threads=config["minimap2"].get("pigz_threads", 4),
        split_prefix=lambda wildcards: f"{RESULTS_DIR}/contaminant_filter/stage2_paf/{wildcards.chunk}.split_idx",
    shell:
        r"""
        mkdir -p $(dirname {output.paf})
        mkdir -p $(dirname {log})

        if [ "$(zcat {input.fastq} | head -c 1 | wc -c)" -eq 0 ]; then
            echo "Empty candidate FASTQ; writing empty PAF." >{log}
            : | pigz -p {params.pigz_threads} -c >{output.paf}
            exit 0
        fi

        minimap2 \
            -t {threads} \
            -I {params.index_size} \
            --split-prefix {params.split_prefix} \
            -x {params.preset} \
            --secondary=no \
            {input.assembly} \
            {input.fastq} \
            2>{log} \
            | pigz -p {params.pigz_threads} -c >{output.paf}
        """


rule confirm_contaminant_read_ids:
    """
    Stage 2 decision. A read is removed only if its PRIMARY alignment is to a
    contaminant contig AND clears the identity and query-coverage thresholds.

    Identity is nmatch/alen (PAF columns 10/11) and coverage is the aligned
    query span over read length, summed across primary blocks to handle split
    alignments. Reads whose primary alignment lands on salamander sequence are
    rescued here regardless of how well they also matched the bait.

    Anything ambiguous — no alignment, low identity, poor coverage — is kept.
    The asymmetry is deliberate: keeping a contaminant read costs a little
    assembly noise, dropping a real read costs sequence you cannot recover.
    """
    input:
        paf=f"{RESULTS_DIR}/contaminant_filter/stage2_paf/{{chunk}}.paf.gz",
        contam_ids=CONTAM_IDS,
    output:
        ids=f"{RESULTS_DIR}/contaminant_filter/confirmed/{{chunk}}.contaminant_read_ids.txt",
    log:
        f"{RESULTS_DIR}/logs/contaminant_filter/confirmed/{{chunk}}.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/contaminant_filter/confirmed/{{chunk}}.tsv"
    params:
        min_identity=CONTAM.get("min_identity", 0.99),
        min_coverage=CONTAM.get("min_query_coverage", 0.90),
        min_mapq=CONTAM.get("min_mapq", 0),
    shell:
        r"""
        mkdir -p $(dirname {output.ids})
        mkdir -p $(dirname {log})

        zcat {input.paf} \
            | awk \
                -v minid={params.min_identity} \
                -v mincov={params.min_coverage} \
                -v minmapq={params.min_mapq} \
                '
            FNR == NR {{
                if ($0 !~ /^#/ && NF > 0) contam[$1] = 1
                next
            }}
            {{
                tp = ""
                for (i = 13; i <= NF; i++) {{
                    if ($i ~ /^tp:A:/) tp = substr($i, 6, 1)
                }}
                if (tp != "P") next
                if ($12 < minmapq) next

                read = $1
                qlen[read] = $2
                aligned[read] += ($4 - $3)
                nmatch[read] += $10
                alen[read] += $11

                if ($10 > best_score[read]) {{
                    best_score[read] = $10
                    best_target[read] = $6
                }}
            }}
            END {{
                for (read in qlen) {{
                    if (!(best_target[read] in contam)) continue
                    if (alen[read] == 0 || qlen[read] == 0) continue
                    if (nmatch[read] / alen[read] < minid) continue
                    if (aligned[read] / qlen[read] < mincov) continue
                    print read
                }}
            }}
            ' {input.contam_ids} - \
            | sort -u \
                >{output.ids} \
                2>{log}

        echo "confirmed=$(wc -l <{output.ids})" >>{log}
        """


rule filter_contaminant_reads:
    input:
        hifi=lambda wildcards: HIFI_READ_BY_CHUNK[wildcards.chunk],
        ids=f"{RESULTS_DIR}/contaminant_filter/confirmed/{{chunk}}.contaminant_read_ids.txt",
    output:
        fastq=f"{RESULTS_DIR}/contaminant_filter/filtered_reads/{{chunk}}.filtered.fastq.gz",
    log:
        f"{RESULTS_DIR}/logs/contaminant_filter/filter_reads/{{chunk}}.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/contaminant_filter/filter_reads/{{chunk}}.tsv"
    threads: CONTAM.get("seqkit_threads", 8)
    params:
        pigz_threads=config["minimap2"].get("pigz_threads", 4),
    shell:
        r"""
        mkdir -p $(dirname {output.fastq})
        mkdir -p $(dirname {log})

        if [ ! -s {input.ids} ]; then
            echo "Nothing to remove; passing reads through." >{log}
            ln -sf $(realpath {input.hifi}) {output.fastq}
            exit 0
        fi

        seqkit grep -j {threads} -v -f {input.ids} {input.hifi} \
            | pigz -p {params.pigz_threads} -c >{output.fastq} \
                2>{log}

        test -s {output.fastq}
        """


rule contaminant_filter_report:
    """
    Per-chunk accounting: how many reads were baited at stage 1, how many
    survived the competitive check, and how many were rescued because their
    primary alignment preferred salamander sequence. The rescue count is the
    number worth reporting in the manuscript — it quantifies how much conserved
    sequence a naive bait-and-delete filter would have discarded.
    """
    input:
        candidates=expand(
            f"{RESULTS_DIR}/contaminant_filter/candidates/{{chunk}}.candidate_ids.txt",
            chunk=HIFI_CHUNKS,
        ),
        confirmed=expand(
            f"{RESULTS_DIR}/contaminant_filter/confirmed/{{chunk}}.contaminant_read_ids.txt",
            chunk=HIFI_CHUNKS,
        ),
    output:
        report=FILTER_REPORT,
    params:
        chunks=HIFI_CHUNKS,
        candidates_dir=f"{RESULTS_DIR}/contaminant_filter/candidates",
        confirmed_dir=f"{RESULTS_DIR}/contaminant_filter/confirmed",
    shell:
        r"""
        mkdir -p $(dirname {output.report})

        printf "chunk\tcandidates\tconfirmed\trescued\n" >{output.report}

        TOTAL_CAND=0
        TOTAL_CONF=0

        for CHUNK in {params.chunks}; do
            CAND=$(wc -l <{params.candidates_dir}/$CHUNK.candidate_ids.txt)
            CONF=$(wc -l <{params.confirmed_dir}/$CHUNK.contaminant_read_ids.txt)
            RESCUED=$((CAND - CONF))

            printf "%s\t%s\t%s\t%s\n" "$CHUNK" "$CAND" "$CONF" "$RESCUED" \
                >>{output.report}

            TOTAL_CAND=$((TOTAL_CAND + CAND))
            TOTAL_CONF=$((TOTAL_CONF + CONF))
        done

        printf "TOTAL\t%s\t%s\t%s\n" \
            "$TOTAL_CAND" "$TOTAL_CONF" "$((TOTAL_CAND - TOTAL_CONF))" \
            >>{output.report}

        cat {output.report}
        """
