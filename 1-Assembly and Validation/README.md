# Assembly and Validation

[![Status: In Progress](https://img.shields.io/badge/status-active-success.svg)](#)
[![Pipeline: Genomics](https://img.shields.io/badge/pipeline-genome__assembly-blue.svg)](#)

This directory contains the scripts, cluster submission workflows, and custom utilities used to generate, scaffold, and validate the genome assembly for this publication. 

The workflow is divided into two major phases: **Assembly & Scaffolding** followed by **Assembly Validation**.

---

## Method 1: Assembly and Scaffolding

This phase focuses on merging initial contigs, scaffolding sequences using reference genomes or genetic maps, and managing long/short-read alignment workflows.

### Scripts & Workflows

| Script / File | Type | Description |
| :--- | :--- | :--- |
| `QuickmergeCommands.txt` | Reference / Notes | *Shell commands for structurally aligning and merging two Y chromosome assemblies using MUMmer and quickmerge. This workflow is used to improve overall genome assembly contiguity.* |
| `How2Winnow.txt` | Reference / Notes | *Command-line examples for sequence alignment using Winnowmap and Minimap2. Used for synteny analysis between repetitive Y chromosome assemblies* |
| `ragtagA3.slurm` | SLURM Script | *A Slurm batch script that scaffolds the YA3 assembly contigs using our merged ISO-1 assembly as a reference using RagTag.* |
| `ragtagA4.slurm` | SLURM Script | *A Slurm batch script that scaffolds the YA4 assembly contigs using our scaffolded A3 assembly as a reference using RagTag* |
| `ReadMapping.slurm` | SLURM Script | *A Slurm Script that maps long reads to our assembly using Winnowmap with repeat-aware parameters. Used for validating read support for chromosome assemblies* |


### Part 2: Method 2 (Assembly Validation)

---

## Method 2: Assembly Validation

Once the scaffolding is finalized, these scripts evaluate the structural integrity, continuity, base accuracy, and coverage distribution of the resulting genome sequence.

### Scripts & Workflows

| Script / File | Type | Description |
| :--- | :--- | :--- |
| `samtoolsdepth.slurm` | SLURM Script | *Calculates per-base coverage across all genomic coordinates using samtools. Generates the raw text input required for downstream depth profiling and gap analysis scripts.* |
| `ReadValidation.slurm` | SLURM Script | *Isolates primary read alignments and uses BEDtools to identify regions completely lacking coverage, subsequently masking them.* |
| `CoveragePlots.py` | Python Utility | *Calculates average read depth across scaffolds using 2kb sliding windows.* |
| `DotplotsCustom.py` | Python Utility | *Generates dotplots from PAF alignment files to compare Y chromosome assemblies against the ISO-1 reference.* |
| `GapFinder.py` | Python Utility | *Scans coverage depth files to identify continuous genomic regions with zero read support. Outputs a coordinate table (CSV) and a size distribution histogram (PDF) for the identified gaps.* |


---

## Cluster Environment Note
> **HPRC Requirements:** All `.slurm` scripts in this directory are configured for the cluster. Module names, version, and dependencies may vary between clusters. 