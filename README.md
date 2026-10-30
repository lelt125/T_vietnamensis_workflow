# Tylototriton vietnamensis genome assembly pipeline

Assembly pipeline for *de novo* genome assembly, scaffolding and quality control of *Tylototriton vietnamensis* (Vietnamese crocodile newt), an endemic, CITES-listed and conservation-priority salamander. The pipeline covers raw PacBio HiFi reads through assembly, quality checks, duplicate purging and scaffolding.

## Overview
- Species: *Tylototriton vietnamensis* (Vietnamese crocodile newt)
- Genome scale: ~24 Gb haploid, 12 chromosomes
- Sequencing: PacBio HiFi (~23x coverage)
- Workflow engine: Snakemake (9.x), config-driven, conda-managed environments

Assembly QC and duplicate purging run as Snakemake workflows. Contig assembly (hifiasm) and Hi-C scaffolding (bwa-mem2 + pairtools + YaHS) are documented as

## Requirements
- conda / mamba (per-rule tool environments via `--sdm conda`)
- uv (Python environment management for the analysis/notebook stack)
- quarto (for running and viewing Quarto notebooks)

Command-line tools invoked by the rules (installed via conda environments per rule, or available on $PATH):

| Stage | Tools |
|---|---|
| Contig assembly | `hifiasm` |
| Contaminant screening | `FCS-GX`, `FCS-adaptor` |
| Duplicate purging | `purge_dups` (`pbcstat`, `calcuts`, `split_fa`, `purge_dups`, `get_seqs`) |
| Hi-C mapping | `bwa-mem2`, `pairtools`, `samtools` |
| Scaffolding | `yahs` |
| Assembly stats | `seqkit` |
| Completeness | `busco` |
| K-mer QC | `meryl`, `merqury.sh` |
| Read mapping | `minimap2` / `mm2plus`, `samtools`, `pigz` |
| Coverage | `mosdepth` |

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

### Resource blocks

Each tool gets its own block for `threads`, `mem_gb`, and (where relevant) `java_heap`:
`seqkit`, `busco`, `meryl`, `merqury`, `minimap2`, `purge_dups`, `mosdepth`.

### purge_dups cutoffs

`purge_dups_tests` defines named cutoff sets (`low` / `mid` / `high`, optional `get_seqs_mode`), letting you run several manual cutoff configurations in parallel. Manual cutoffs are used here because automatic valley-finding is unreliable on this genome's coverage histogram.

## Usage
### Assembly and scaffolding commands
```
# Contig-level genome assembly
hifiasm -o T_vietnamensis -t 32 ${HIFI_READS_DIR}/*.hifi_reads.fastq.gz
awk '/^S/{print ">"$2; print $3}' T_vietnamensis.bp.p_ctg.gfa > T_vietnamensis.p_ctg.fa

# Index the contigs after removing contaminants, purging duplicates with purge_dups.smk
# and filtering short contigs
samtools faidx ${CONTIGS}
cut -f1,2 ${CONTIGS}.fai > contigs.chrom.sizes      # needed by pairtools parse
bwa-mem2 index -p ${CONTIGS} ${CONTIGS} > bwa_mem2_index.log

# Map, parse to pairs, deduplicate
mkdir -p tmp_pairtools

bwa-mem2 mem -5SP -T0 -t ${THREADS} ${CONTIGS} ${HIC_R1} ${HIC_R2} \
    2> hic_bwa_mem2.log \
| pairtools parse \
    --min-mapq 10 \
    --walks-policy 5unique \
    --max-inter-align-gap 30 \
    --chroms-path contigs.chrom.sizes \
    --nproc-in 8 --nproc-out 8 \
| pairtools sort \
    --tmpdir tmp_pairtools \
    --nproc 16 \
| pairtools dedup \
    --output-stats hic_pairtools.dedup_stats.txt \
    --nproc-in 8 --nproc-out 8 \
| pairtools split \
    --output-pairs hic_to_contigs.pairs.gz \
    --output-sam - \
| samtools view -@ 16 -b -o hic_to_contigs.pairtools.bam -

# Scaffolding with YaHS (Omni-C: no restriction motif)
yahs -o yahs_out ${CONTIGS} hic_to_contigs.pairtools.bam > yahs.log 2>&1

# Use genome_qc.smk to assess genome assembly
```

### Snakemake pipeline
For each section, adjust the number of processing cores and RAM memory according to your
own requirements. Remove the -n flag to begin the actual run.

#### Run genome QC (seqkit, BUSCO, Merqury, mosdepth coverage)
```
snakemake -s workflow/genome_qc.smk --configfile config_production.yaml --sdm conda \
    --cores 40 --resources mem_gb=300 -p -n
```

#### Run duplicate purging across all configured cutoff tests
```
snakemake -s workflow/purge_dups.smk --configfile config_production.yaml --sdm conda \
    --cores 96 --resources mem_gb=900 -p -n
```
## Citation
If you use this pipeline, please cite the associated Data Descriptor (citation to be added upon publication in Scientific Data).

## License
This repository is licensed under the [MIT License](LICENSE). You are free to use, modify, and distribute this pipeline with attribution.