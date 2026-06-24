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