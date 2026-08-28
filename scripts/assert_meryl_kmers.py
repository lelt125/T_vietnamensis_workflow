"""Assert a meryl database holds the expected k-mer count.

Invoked via Snakemake's `script:` directive. Two modes, chosen by
`snakemake.params.mode`:

  seqkit  - expected count derived from a `seqkit stats -T` report
            (requires params.kmer)
  sum     - expected count is the sum of `present` across several
            upstream statistics files (snakemake.input.parts)

Guards against meryl's silent batch failure: a database that is
structurally valid but holds only a fraction of the expected k-mers.
"""

from pathlib import Path


def parse_present(stats_path: str) -> int:
    for line in Path(stats_path).read_text().splitlines():
        fields = line.split()
        if len(fields) >= 2 and fields[0] == "present":
            return int(fields[1].replace(",", ""))
    raise ValueError(f"no 'present' line found in {stats_path}")


def expected_from_seqkit(seqkit_path: str, kmer: int) -> int:
    rows = [
        line.split("\t")
        for line in Path(seqkit_path).read_text().splitlines()
        if line.strip()
    ]
    if len(rows) < 2:
        raise ValueError(f"no data rows in {seqkit_path}")
    header = rows[0]
    num_seqs_col = header.index("num_seqs")
    sum_len_col = header.index("sum_len")
    return sum(
        int(row[sum_len_col]) - int(row[num_seqs_col]) * (kmer - 1)
        for row in rows[1:]
    )


def main() -> None:
    stats = snakemake.input.stats
    label = snakemake.params.label
    tolerance = float(snakemake.params.get("tolerance", 0.0))
    mode = snakemake.params.mode

    present = parse_present(stats)

    if mode == "seqkit":
        kmer = int(snakemake.params.kmer)
        expected = expected_from_seqkit(snakemake.input.seqkit, kmer)
        basis = "sequence length"
    elif mode == "sum":
        expected = sum(parse_present(p) for p in snakemake.input.parts)
        basis = "sum of parts"
    else:
        raise ValueError(f"unknown mode: {mode!r}")

    with open(snakemake.log[0], "w") as log:
        if expected <= 0:
            msg = f"FATAL [{label}]: nonsensical expected count ({expected})"
            print(msg, file=log)
            raise SystemExit(msg)

        deviation = abs(present - expected) / expected
        print(
            f"[{label}] present={present} expected={expected} "
            f"({basis}) rel_diff={deviation:.6f}",
            file=log,
        )

        if deviation > tolerance:
            msg = (
                f"FATAL [{label}]: k-mer count deviates by {deviation:.4%} "
                f"(tolerance {tolerance:.4%}). A shortfall near 50% indicates "
                f"the silent meryl batch failure."
            )
            print(msg, file=log)
            raise SystemExit(msg)

    Path(snakemake.output[0]).write_text(
        f"present={present} expected={expected} rel_diff={deviation:.6f}\n"
    )


main()