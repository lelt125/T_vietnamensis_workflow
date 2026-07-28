import json
from pathlib import Path
from typing import ClassVar

import pandas as pd


class BuscoSummarizer:
    METRIC_MAP: ClassVar[dict[str, tuple[str, str]]] = {
        "complete_pct": ("results", "Complete percentage"),
        "complete_n": ("results", "Complete BUSCOs"),
        "single_copy_pct": ("results", "Single copy percentage"),
        "single_copy_n": ("results", "Single copy BUSCOs"),
        "multi_copy_pct": ("results", "Multi copy percentage"),
        "multi_copy_n": ("results", "Multi copy BUSCOs"),
        "fragmented_pct": ("results", "Fragmented percentage"),
        "fragmented_n": ("results", "Fragmented BUSCOs"),
        "missing_pct": ("results", "Missing percentage"),
        "missing_n": ("results", "Missing BUSCOs"),
        "n_markers": ("results", "n_markers"),
        "avg_identity": ("results", "avg_identity"),
        "internal_stop_codon_count": ("results", "internal_stop_codon_count"),
        "internal_stop_codon_pct": ("results", "internal_stop_codon_percent"),
        "n_scaffolds": ("metrics", "Number of scaffolds"),
        "n_contigs": ("metrics", "Number of contigs"),
        "total_length_bp": ("metrics", "Total length"),
        "percent_gaps": ("metrics", "Percent gaps"),
        "scaffold_n50": ("metrics", "Scaffold N50"),
        "contig_n50": ("metrics", "Contigs N50"),
    }
    NUMERIC_COLS: ClassVar[list[str]] = [
        "n_scaffolds",
        "n_contigs",
        "total_length_bp",
        "scaffold_n50",
        "contig_n50",
    ]
    CONDENSED_COLS: ClassVar[list[str]] = [
        "sample",
        "complete_pct",
        "single_copy_pct",
        "multi_copy_pct",
        "fragmented_pct",
        "missing_pct",
        "n_contigs",
        "total_length_bp",
        "contig_n50",
    ]

    def __init__(self, busco_dir: Path, lineage: str, run_names: list[str]):
        self.busco_dir = Path(busco_dir)
        self.lineage = lineage
        self.run_names = run_names

    def _json_path(self, run_name: str) -> Path:
        path = (
            self.busco_dir
            / run_name
            / self.lineage
            / f"short_summary.specific.{self.lineage}.{run_name}_busco.json"
        )
        if not path.is_file():
            raise FileNotFoundError(f"BUSCO short_summary JSON not found: {path}")
        return path

    def _parse_one(self, run_name: str) -> dict:
        json_path = self._json_path(run_name)
        with open(json_path) as f:
            data = json.load(f)

        row = {
            "sample": run_name,
            "lineage": data.get("lineage_dataset", {}).get("name"),
        }
        row.update(
            {
                key: data.get(section, {}).get(field)
                for key, (section, field) in self.METRIC_MAP.items()
            }
        )
        row["one_line_summary"] = data.get("results", {}).get("one_line_summary")
        row["busco_version"] = data.get("versions", {}).get("busco")
        row["gene_predictor"] = data.get("parameters", {}).get("gene_predictor")
        row["json_path"] = str(json_path)
        return row

    def summary(self) -> pd.DataFrame:
        df = pd.DataFrame(self._parse_one(run_name) for run_name in self.run_names)
        df[self.NUMERIC_COLS] = df[self.NUMERIC_COLS].apply(
            pd.to_numeric, errors="coerce"
        )
        df["percent_gaps"] = df["percent_gaps"].str.rstrip("%").astype(float)
        return df

    def condensed_summary(self) -> pd.DataFrame:
        df = self.summary()[self.CONDENSED_COLS].rename(
            columns={
                "complete_pct": "complete_%",
                "single_copy_pct": "single_%",
                "multi_copy_pct": "duplicated_%",
                "fragmented_pct": "fragmented_%",
                "missing_pct": "missing_%",
                "total_length_bp": "total_length_gb",
                "contig_n50": "contig_n50_mb",
            }
        )
        df["total_length_gb"] = (df["total_length_gb"] / 1e9).round(2)
        df["contig_n50_mb"] = (df["contig_n50_mb"] / 1e6).round(2)
        return df
