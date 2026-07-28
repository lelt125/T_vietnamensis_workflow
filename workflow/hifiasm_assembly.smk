# rules/hifiasm_assembly.smk
#
# De novo assembly with hifiasm.
# Supports Hi-C phased mode (--h1/--h2) and HiFi-only mode, toggled via config.
#
# Hi-C mode is recommended for T. vietnamensis: the Hi-C phasing resolves the
# heterozygous allelic duplication that purge_dups (coverage-based) cannot,
# producing two fully-phased haplotype assemblies rather than a het-inflated primary.

include: "rules/common_config_files.smk"

RUN_NAME, RESULTS_BASE, RESULTS_DIR = get_results_dir(config)
ASSEMBLIES, SAMPLES, FASTA_BY_SAMPLE = load_assemblies(config)
HIFI_READS_DIR, HIFI_READS, HIFI_READ_BY_CHUNK, HIFI_CHUNKS = load_hifi_reads(config)

# --- Hi-C read loading (only required when hifiasm.mode == "hic") ---------
HIFIASM_MODE = config["hifiasm"].get("mode", "hic")   # "hic" or "hifi_only"

if HIFIASM_MODE == "hic":
    HIC_R1_LIST, HIC_R2_LIST = load_hic_reads(config)
    HIC_R1 = ",".join(HIC_R1_LIST)   # comma-joined for a single --h1 flag
    HIC_R2 = ",".join(HIC_R2_LIST)
else:
    HIC_R1_LIST, HIC_R2_LIST = [], []
    HIC_R1, HIC_R2 = "", ""

# The sample name for the assembly being produced (distinct from the input
# assembly samples in FASTA_BY_SAMPLE — this rule file *creates* an assembly).
ASM_NAME = config["hifiasm"]["asm_name"]

HIFI_INPUT = " ".join(str(r) for r in HIFI_READS)


# --- Expected outputs -----------------------------------------------------
# hifiasm output filenames depend on mode:
#   hic mode : {prefix}.hic.hap1.p_ctg.gfa, {prefix}.hic.hap2.p_ctg.gfa  (NO primary p_ctg)
#   bp  mode : {prefix}.bp.p_ctg.gfa, {prefix}.bp.hap1.p_ctg.gfa, {prefix}.bp.hap2.p_ctg.gfa

if HIFIASM_MODE == "hic":
    ASM_INFIX = "hic"
    PRIMARY_TAGS = ["hap1.p_ctg", "hap2.p_ctg"]
else:
    ASM_INFIX = "bp"
    PRIMARY_TAGS = ["p_ctg", "hap1.p_ctg", "hap2.p_ctg"]

HIFIASM_GFAS = expand(
    f"{RESULTS_DIR}/hifiasm/{ASM_NAME}/{ASM_NAME}.{ASM_INFIX}.{{tag}}.gfa",
    tag=PRIMARY_TAGS
)

HIFIASM_FASTAS = expand(
    f"{RESULTS_DIR}/hifiasm/{ASM_NAME}/{ASM_NAME}.{ASM_INFIX}.{{tag}}.fa",
    tag=PRIMARY_TAGS
)


rule all:
    input:
        HIFIASM_GFAS,
        HIFIASM_FASTAS


rule hifiasm_assemble:
    """Run hifiasm. Produces GFA graphs for each haplotype (+ primary in bp mode).

    All HiFi reads are passed at once; hifiasm handles its own parallelism.
    The prefix (-o) drives all output filenames.
    """
    input:
        hifi=HIFI_READS,
        # Hi-C inputs only tracked when in hic mode; empty list otherwise so
        # Snakemake doesn't demand files that aren't used.
        hic_r1=(config["hifiasm"].get("hic_r1", []) if HIFIASM_MODE == "hic" else []),
        hic_r2=(config["hifiasm"].get("hic_r2", []) if HIFIASM_MODE == "hic" else [])
    output:
        gfas=HIFIASM_GFAS,
        # hifiasm also writes a bunch of intermediates (.ec.bin, .ovlp.*.bin);
        # we don't track them as outputs but they land in the same dir.
        marker=f"{RESULTS_DIR}/hifiasm/{ASM_NAME}/{ASM_NAME}.hifiasm.done"
    threads:
        config["hifiasm"]["threads"]
    resources:
        mem_gb=config["hifiasm"]["mem_gb"]
    params:
        prefix=lambda wildcards: f"{RESULTS_DIR}/hifiasm/{ASM_NAME}/{ASM_NAME}",
        hifi_input=HIFI_INPUT,
        mode=HIFIASM_MODE,
        hic_r1=HIC_R1,
        hic_r2=HIC_R2,
        # -l0 disables hifiasm's OWN purge for hom samples; for a highly het
        # wild sample you generally want purge ON (default -l2/-l3). Leave to
        # config so you can experiment. Default: don't override (hifiasm auto).
        purge_level=config["hifiasm"].get("purge_level", ""),
        extra=config["hifiasm"].get("extra", "")
    log:
        f"{RESULTS_DIR}/logs/hifiasm/{ASM_NAME}.log"
    benchmark:
        f"{RESULTS_DIR}/benchmarks/hifiasm/{ASM_NAME}.tsv"
    shell:
        r"""
        mkdir -p $(dirname {params.prefix})
        mkdir -p $(dirname {log})

        # Build optional flag strings
        PURGE_FLAG=""
        if [ -n "{params.purge_level}" ]; then
            PURGE_FLAG="-l {params.purge_level}"
        fi

        HIC_FLAGS=""
        if [ "{params.mode}" = "hic" ]; then
            HIC_FLAGS="--h1 {params.hic_r1} --h2 {params.hic_r2}"
        fi

        hifiasm \
            -o {params.prefix} \
            -t {threads} \
            $PURGE_FLAG \
            $HIC_FLAGS \
            {params.extra} \
            {params.hifi_input} \
            > {log} 2>&1

        # Verify expected GFAs exist and are non-empty
        for gfa in {output.gfas}; do
            test -s "$gfa"
        done

        touch {output.marker}
        """


rule hifiasm_gfa_to_fasta:
    """Convert a hifiasm p_ctg GFA to FASTA.

    GFA 'S' lines carry the segment sequence in column 3; emit as FASTA.
    Wildcard {tag} matches e.g. 'hic.hap1.p_ctg'.
    """
    input:
        gfa=f"{RESULTS_DIR}/hifiasm/{ASM_NAME}/{ASM_NAME}.{{tag}}.gfa"
    output:
        fasta=f"{RESULTS_DIR}/hifiasm/{ASM_NAME}/{ASM_NAME}.{{tag}}.fa"
    log:
        f"{RESULTS_DIR}/logs/hifiasm/gfa_to_fasta/{ASM_NAME}.{{tag}}.log"
    benchmark:
        f"{RESULTS_DIR}/benchmarks/hifiasm/gfa_to_fasta/{ASM_NAME}.{{tag}}.tsv"
    shell:
        r"""
        mkdir -p $(dirname {output.fasta})
        mkdir -p $(dirname {log})

        awk '/^S/{{ print ">"$2"\n"$3 }}' {input.gfa} \
            > {output.fasta} \
            2> {log}

        test -s {output.fasta}
        """