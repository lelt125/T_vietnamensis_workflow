#!/usr/bin/env python3
"""
summarize_benchmarks.py

Aggregate Snakemake `benchmark:` TSV files into a single summary table.

Usage:
    python summarize_benchmarks.py <benchmarks_dir> [--out summary.tsv] [--sort-by max_rss]

Example:
    python summarize_benchmarks.py \
        data/2024_Tylototriton_vietnamensis_Wild_Genomes/Linh_Analysis/results/production_run/benchmarks \
        --out benchmark_summary.tsv

Assumes the standard Snakemake benchmark TSV columns:
    s  h:m:s  max_rss  max_vms  max_uss  max_pss  io_in  io_out  mean_load  cpu_time

The rule name and wildcards are inferred from the file's path relative to
the benchmarks directory, e.g.:
    benchmarks/read_mapping/hifi_bam/04_hifi_asm_.../XDOVE_....tsv
        -> rule   = read_mapping/hifi_bam
        -> sample = 04_hifi_asm_...
        -> chunk  = XDOVE_....tsv (stem)
This is a best-effort label, not a strict re-derivation of Snakemake wildcards.
"""

import argparse
import sys
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    sys.exit(
        "This script requires pandas. Install it with:\n"
        "    pip install pandas --break-system-packages\n"
        "or activate the conda env that has it."
    )


def infer_labels(tsv_path: Path, base_dir: Path):
    """Best-effort extraction of a rule label and wildcard-ish path parts."""
    rel = tsv_path.relative_to(base_dir)
    parts = rel.parts  # e.g. ("read_mapping", "hifi_bam", "sample", "chunk.tsv")

    if len(parts) == 1:
        rule = "unknown"
        wildcards = parts[0].removesuffix(".tsv")
    else:
        rule = "/".join(parts[:-1])
        wildcards = parts[-1].removesuffix(".tsv")

    return rule, wildcards


def load_benchmarks(base_dir: Path) -> pd.DataFrame:
    tsv_files = sorted(base_dir.rglob("*.tsv"))
    if not tsv_files:
        sys.exit(f"No .tsv files found under {base_dir}")

    rows = []
    for tsv_path in tsv_files:
        try:
            df = pd.read_csv(tsv_path, sep="\t")
        except Exception as e:
            print(f"WARNING: could not read {tsv_path}: {e}", file=sys.stderr)
            continue

        if df.empty:
            print(f"WARNING: empty benchmark file, skipping: {tsv_path}", file=sys.stderr)
            continue

        rule, wildcards = infer_labels(tsv_path, base_dir)
        row = df.iloc[0].to_dict()
        row["rule"] = rule
        row["wildcards"] = wildcards
        row["file"] = str(tsv_path)
        rows.append(row)

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description="Summarize Snakemake benchmark TSVs.")
    parser.add_argument(
        "benchmarks_dir",
        type=Path,
        help="Directory to search recursively for benchmark .tsv files.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Optional path to write the full summary table (TSV). "
             "If omitted, only prints to stdout.",
    )
    parser.add_argument(
        "--sort-by",
        default="max_rss",
        help="Column to sort by, descending (default: max_rss). "
             "Common choices: max_rss, s, cpu_time, mean_load.",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=20,
        help="Number of rows to print per-rule and overall (default: 20). "
             "Use -1 to print all rows.",
    )
    args = parser.parse_args()

    base_dir = args.benchmarks_dir
    if not base_dir.is_dir():
        sys.exit(f"Not a directory: {base_dir}")

    df = load_benchmarks(base_dir)

    if args.sort_by not in df.columns:
        sys.exit(
            f"--sort-by '{args.sort_by}' not found in benchmark columns: "
            f"{list(df.columns)}"
        )

    df = df.sort_values(args.sort_by, ascending=False).reset_index(drop=True)

    display_cols = [
        "rule", "wildcards", "s", "h:m:s", "max_rss", "max_vms",
        "mean_load", "cpu_time", "io_in", "io_out",
    ]
    display_cols = [c for c in display_cols if c in df.columns]

    n = None if args.top == -1 else args.top

    print(f"\n=== Top {args.top if n else 'ALL'} jobs by {args.sort_by} (descending) ===\n")
    print(df[display_cols].head(n).to_string(index=False))

    print("\n=== Per-rule peak resource usage ===\n")
    per_rule = (
        df.groupby("rule")
        .agg(
            n_jobs=("rule", "count"),
            max_rss_MB=("max_rss", "max"),
            mean_rss_MB=("max_rss", "mean"),
            max_runtime_s=("s", "max"),
            mean_runtime_s=("s", "mean"),
            max_cpu_time_s=("cpu_time", "max") if "cpu_time" in df.columns else ("s", "max"),
        )
        .sort_values("max_rss_MB", ascending=False)
    )
    print(per_rule.to_string())

    if args.out:
        df.to_csv(args.out, sep="\t", index=False)
        print(f"\nFull summary written to: {args.out}")


if __name__ == "__main__":
    main()