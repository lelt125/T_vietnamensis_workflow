from pathlib import Path
import pandas as pd


def load_assemblies(config):
    assemblies = pd.read_csv(config["assemblies_csv"], comment="#")

    required_columns = {"sample_id", "fasta"}
    missing_columns = required_columns - set(assemblies.columns)

    if missing_columns:
        raise ValueError(
            "assemblies_csv is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    if assemblies["sample_id"].duplicated().any():
        dupes = assemblies.loc[
            assemblies["sample_id"].duplicated(),
            "sample_id"
        ].tolist()
        raise ValueError(f"Duplicate sample_id values in assemblies_csv: {dupes}")

    samples = assemblies["sample_id"].astype(str).tolist()
    fasta_by_sample = dict(zip(assemblies["sample_id"], assemblies["fasta"]))

    missing_fastas = [
        fasta for fasta in fasta_by_sample.values()
        if not Path(fasta).exists()
    ]

    if missing_fastas:
        raise ValueError(
            "Some assembly FASTA files do not exist:\n"
            + "\n".join(missing_fastas)
        )

    return assemblies, samples, fasta_by_sample


def load_hifi_reads(config):
    hifi_reads_dir = Path(config["hifi_reads_dir"]).resolve()

    hifi_reads = sorted(
        hifi_reads_dir.glob("*.hifi_reads.fastq.gz")
    )

    if len(hifi_reads) == 0:
        raise ValueError(f"No HiFi read files found in {hifi_reads_dir}")

    hifi_read_by_chunk = {
        read.name.replace(".hifi_reads.fastq.gz", ""): str(read.resolve())
        for read in hifi_reads
    }

    hifi_chunks = list(hifi_read_by_chunk.keys())

    return hifi_reads_dir, hifi_reads, hifi_read_by_chunk, hifi_chunks


def load_hic_reads(config):
    """Load Hi-C paired-end reads for hifiasm Hi-C phasing (--h1 / --h2).

    Supports two config layouts under config["hifiasm"]:

    1. Explicit file lists (preferred, unambiguous R1/R2 pairing):
           hifiasm:
             hic_r1: /path/lib1_R1.fastq.gz          # str or list of str
             hic_r2: /path/lib1_R2.fastq.gz

       Multiple libraries:
             hic_r1: [/path/lib1_R1.fastq.gz, /path/lib2_R1.fastq.gz]
             hic_r2: [/path/lib1_R2.fastq.gz, /path/lib2_R2.fastq.gz]

    2. Directory glob (like load_hifi_reads), pairing by _R1/_R2 in filenames:
           hifiasm:
             hic_reads_dir: /path/to/hic

    Returns:
        (hic_r1_list, hic_r2_list): sorted lists of resolved absolute-path
        strings. R1[i] pairs with R2[i]. Both lists have equal length.

    hifiasm accepts comma-separated files for a single --h1/--h2 flag, so the
    caller typically joins each list with ",".
    """
    hifiasm_cfg = config.get("hifiasm", {})

    def _resolve_existing(paths, label):
        resolved = []
        missing = []
        for p in paths:
            path = Path(p)
            if path.exists():
                resolved.append(str(path.resolve()))
            else:
                missing.append(str(p))
        if missing:
            raise ValueError(
                f"Hi-C {label} files do not exist:\n" + "\n".join(missing)
            )
        return resolved

    # --- Layout 1: explicit hic_r1 / hic_r2 ---
    if "hic_r1" in hifiasm_cfg or "hic_r2" in hifiasm_cfg:
        if "hic_r1" not in hifiasm_cfg or "hic_r2" not in hifiasm_cfg:
            raise ValueError(
                "Both config['hifiasm']['hic_r1'] and ['hic_r2'] must be set "
                "when using explicit Hi-C file paths"
            )

        r1_raw = hifiasm_cfg["hic_r1"]
        r2_raw = hifiasm_cfg["hic_r2"]

        # Normalize single string -> single-element list
        r1_list = [r1_raw] if isinstance(r1_raw, str) else list(r1_raw)
        r2_list = [r2_raw] if isinstance(r2_raw, str) else list(r2_raw)

        if len(r1_list) != len(r2_list):
            raise ValueError(
                f"Hi-C R1 and R2 must have equal file counts: "
                f"{len(r1_list)} R1 vs {len(r2_list)} R2"
            )

        if len(r1_list) == 0:
            raise ValueError("hifiasm.hic_r1 / hic_r2 are empty")

        hic_r1 = _resolve_existing(r1_list, "R1")
        hic_r2 = _resolve_existing(r2_list, "R2")

        return hic_r1, hic_r2

    # --- Layout 2: directory glob ---
    if "hic_reads_dir" in hifiasm_cfg:
        hic_dir = Path(hifiasm_cfg["hic_reads_dir"]).resolve()

        if not hic_dir.exists():
            raise ValueError(f"Hi-C reads dir does not exist: {hic_dir}")

        # Match common R1/R2 naming: *_R1*.f*q.gz and *_R2*.f*q.gz
        r1_files = sorted(
            f for pattern in ("*_R1*.fastq.gz", "*_R1*.fq.gz")
            for f in hic_dir.glob(pattern)
        )
        r2_files = sorted(
            f for pattern in ("*_R2*.fastq.gz", "*_R2*.fq.gz")
            for f in hic_dir.glob(pattern)
        )

        if len(r1_files) == 0 or len(r2_files) == 0:
            raise ValueError(
                f"No Hi-C R1/R2 files found in {hic_dir} "
                f"(expected *_R1*/*_R2* .fastq.gz or .fq.gz)"
            )

        if len(r1_files) != len(r2_files):
            raise ValueError(
                f"Unequal Hi-C R1/R2 file counts in {hic_dir}: "
                f"{len(r1_files)} R1 vs {len(r2_files)} R2. "
                f"Check filenames pair cleanly by _R1/_R2."
            )

        hic_r1 = [str(f.resolve()) for f in r1_files]
        hic_r2 = [str(f.resolve()) for f in r2_files]

        return hic_r1, hic_r2

    raise ValueError(
        "Hi-C mode requires either config['hifiasm']['hic_r1']+['hic_r2'] "
        "or config['hifiasm']['hic_reads_dir']"
    )


def as_bool(value, default=False):
    if pd.isna(value):
        return default

    value = str(value).strip().lower()

    if value == "":
        return default

    if value in {"true", "1", "yes", "y"}:
        return True

    if value in {"false", "0", "no", "n"}:
        return False

    raise ValueError(f"Invalid boolean value in CSV: {value}")


def get_qc_samples(assemblies, column_name, default=False):
    if column_name not in assemblies.columns:
        return assemblies["sample_id"].tolist() if default else []

    values = assemblies[column_name].apply(lambda x: as_bool(x, default=default))

    return assemblies.loc[values, "sample_id"].tolist()


def get_results_dir(config):
    run_name = config["run_name"]
    results_base = config["results_base"]
    results_dir = f"{results_base}/{run_name}"

    return run_name, results_base, results_dir


def load_valid_busco_lineages(config):
    valid_lineages_file = Path(config["busco"]["valid_lineages_file"])

    valid_lineages = {
        line.strip()
        for line in valid_lineages_file.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    }

    default_lineage = config["busco"]["lineage"]

    if default_lineage not in valid_lineages:
        raise ValueError(f"Invalid default BUSCO lineage: {default_lineage}")

    return valid_lineages, default_lineage


def build_busco_lineage_by_sample(assemblies, valid_lineages, default_lineage):
    lineage_by_sample = {}

    for _, row in assemblies.iterrows():
        sample = row["sample_id"]

        if "busco_lineage" not in assemblies.columns:
            lineage = default_lineage
        else:
            lineage = row["busco_lineage"]

            if pd.isna(lineage) or str(lineage).strip() == "":
                lineage = default_lineage
            else:
                lineage = str(lineage).strip()

        if lineage not in valid_lineages:
            raise ValueError(f"Invalid BUSCO lineage for {sample}: {lineage}")

        lineage_by_sample[sample] = lineage

    return lineage_by_sample