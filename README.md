# Tylototriton vietnamensis genome assembly pipeline

Assembly pipeline for *de novo* genome assembly, scaffolding and quality control of *Tylototriton vietnamensis* (Vietnamese crocodile newt), an endangered, CITES-listed and conservation-priority salamander. The pipeline covers raw PacBio HiFi reads through assembly, quality checks, duplicate purging and scaffolding.

## Overview
- Species: *Tylototriton vietnamensis* (Vietnamese crocodile newt)
- Genome scale: ~24 Gb haploid, 12 chromosomes
- Sequencing: PacBio HiFi (~23x coverage)
- Workflow engine: Snakemake (9.x), config-driven, conda-managed environments

This repository documents the assembly QC and duplicate-purging snakemake workflows. Other stages referenced in the broader project (contaminant read filtering, assembly and scaffolding) are included as commands.

## Repository structure
```
.
├── config_production.yaml
├── config_test.yaml
├── notebooks                           # exploratory/diagnostic analyses, not part of the pipeline
│   ├── analyze_contaminants.qmd
│   ├── analyze_purge_dups_output.qmd
│   ├── investigate_duplicates.ipynb
│   └── sum_assembly_stats.qmd
├── pyproject.toml
├── _quarto.yml
├── README.md
├── scripts                             # standalone helper scripts used by Snakemake rules
│   ├── assert_meryl_kmers.py
│   ├── hist_plot.py
│   ├── __init__.py
│   ├── summarize_benchmarks.py
│   └── summarizer.py
├── uv.lock
└── workflow                            # Snakemake rules, configs, and conda envs
    ├── envs
    │   ├── busco.yaml
    │   ├── genome_assessment_minimal.yaml
    │   ├── genome_assessment.yaml
    │   ├── minimap2.yaml
    │   ├── mm2plus.yaml
    │   └── purge_dups.yaml
    ├── filter_contaminant_reads.smk
    ├── genome_qc.smk
    ├── hifiasm_assembly.smk
    ├── purge_dups.smk
    ├── rules
    │   └── common_config_files.smk
```

## Requirements
- Snakemake ≥ 9.x
- conda / mamba
- uv (Python environment management, optional)

Command-line tools invoked by the rules (installed via conda environments per rule, or available on $PATH):

| Stage | Tools |
|---|---|
| Assembly stats | `seqkit` |
| Completeness | `busco` |
| K-mer QC | `meryl`, `merqury.sh` |
| Read mapping | `minimap2` / `mm2plus`, `samtools`, `pigz` |
| Coverage | `mosdepth` |
| Duplicate purging | `purge_dups` (`pbcstat`, `calcuts`, `split_fa`, `purge_dups`, `get_seqs`) |

## Configuration

All paths and resources live in a single YAML config (see `config_production.yaml`).

### Output location

```yaml
run_name: production_run
results_base: data/results
```

Outputs go to `{results_base}/{run_name}`.

### Sample sheet

`assemblies_csv` points to a table with one row per assembly to QC:

| Column | Purpose |
|---|---|
| `sample_id`, `fasta` | required — assembly name and path |
| `run_seqkit`, `run_busco`, `run_merqury`, `run_mosdepth` | optional toggles, per QC tool |
| `busco_lineage` | optional override of the default lineage |

QC steps are toggled per row, not globally — so a raw hifiasm assembly and a purged version can run through different QC subsets.

### HiFi reads

`hifi_reads_dir` points to a directory of chunked FASTQ files (`*.hifi_reads.fastq.gz`). Each file is one "chunk," parallelized independently for mapping and k-mer counting.

### Resource blocks

Each tool gets its own block for `threads`, `mem_gb`, and (where relevant) `java_heap`:
`seqkit`, `busco`, `meryl`, `merqury`, `minimap2`, `purge_dups`, `mosdepth`.

### purge_dups cutoffs

`purge_dups_tests` defines named cutoff sets (`low` / `mid` / `high`, optional `get_seqs_mode`), letting you run several manual cutoff configurations in parallel. Manual cutoffs are used here because automatic valley-finding is unreliable on this genome's coverage histogram.

## Usage
For each section, adjust the number of processing cores and RAM memory according to your
own requirements. Remove the -n flag to begin the actual run.

### Run genome QC (seqkit, BUSCO, Merqury, mosdepth coverage)
```
snakemake -s workflow/genome_qc.smk --configfile config_production.yaml --sdm conda \
    --cores 40 --resources mem_gb=300 -p -n
```

### Run duplicate purging across all configured cutoff tests
```
snakemake -s workflow/purge_dups.smk --configfile config_production.yaml --sdm conda \
    --cores 96 --resources mem_gb=900 -p -n
```
## Citation
If you use this pipeline, please cite the associated Data Descriptor (citation to be added upon publication in Scientific Data).

## License
This repository is licensed under the [MIT License](LICENSE). You are free to use, modify, and distribute this pipeline with attribution.