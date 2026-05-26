# Targeted Assembly of rDNA Structural Variants and Their Characterization

[![Status: Completed](https://img.shields.io/badge/status-complete-success.svg)](#)
[![Analysis: rDNA](https://img.shields.io/badge/analysis-rDNA__variants-red.svg)](#)

This directory contains the pipeline, specialized scripts, reference sequences, and workflows used for the targeted assembly, manual annotation, quality validation, and structural characterization of the ribosomal DNA (rDNA) array repeat units.

---

## Method 1: Targeted Assembly, Validation, and Annotation of rDNA Repeat Variants suing RIBOTIN

This phase utilizes long-read targeted assembly workflows to isolate specific rDNA units, track distinct structural morph variants, and annotate coding boundaries across chromosomes.

### Numbered Execution Pipeline & Notes

| Script / File | Type | Description |
| :--- | :--- | :--- |
| `'1)RMtoBED.py'` | Python Utility | *Converts standard RepeatMasker .out files into a clean 6-column BED format. This facilitates annotation in the next step and is required as the main input* |
| `'2)rDNA_ManualAnnotation.py'` | Python Utility | *Parses RepeatMasker BED files to annotate ribosomal DNA (rDNA) repeats via their RepeatMasker annotations. Identifies canonical and TE-isnerted units. One clean, non-inserted unit should be manually isolated to act as the reference in the ribotin-ref step* |
| `'3)rDNA_mapping.slurm'` | SLURM Script | *Maps long sequence reads to a specific rDNA locus using Winnowmap. It sequentially subsamples the mapped reads to targeted coverage depths (e.g., 5X, 15X, 30X) to compare variant discovery between different ribotin runs* |
| `'4)ribotin-ref.txt'` | Reference / Notes | *A quick reference guide containing the exact command-line parameters for executing ribotin-ref.* |
| `'5)MorphMapping.slurm'` | SLURM Script | *Validates the genomic placement of assembled rDNA morphs back onto the reference assembly using a dual-pass Minimap2 alignment strategy.* |
| `'6)ArrayEstimations.slurm'` | SLURM Script | *Maps the complete, non-subsampled long-read dataset against the resolved morph consensus libraries. Calculates the estimated overall copy number and total base-pair abundance for each specific rDNA variant within the massive array.* |

### Isolated Reference and Variant Assemblies (FASTA)

| Reference Dataset | Type | Target Description |
| :--- | :--- | :--- |
| `ISO1X_rDNA.fasta` | FASTA Reference | *rDNA units extracted from its respective assembly* |
| `ISO1Y_rDNA.fasta` | FASTA Reference | *rDNA units extracted from its respective assembly* |
| `ISO1Xmorphs.fasta` | FASTA Reference | *rDNA variant units assembled by RIBOTIN* |
| `ISO1YMorphs.fasta` | FASTA Reference | *rDNA variant units assembled by RIBOTIN* |
| `A3X_rDNA.fasta` | FASTA Assembly | *rDNA units extracted from its respective assembly* |
| `A3Y_rDNA.fasta` | FASTA Assembly | *rDNA units extracted from its respective assembly* |
| `A3Xmorphs.fasta` | FASTA Assembly | *rDNA variant units assembled by RIBOTIN* |
| `A3Ymorphs.fasta` | FASTA Assembly | *rDNA variant units assembled by RIBOTIN* |
| `A4X_rDNA.fasta` | FASTA Assembly | *rDNA units extracted from its respective assembly* |
| `A4Y_rDNA.fasta` | FASTA Assembly | *rDNA units extracted from its respective assembly* |
| `A4Xmorphs.fasta` | FASTA Assembly | *rDNA variant units assembled by RIBOTIN* |
| `A4Ymorphs.fasta` | FASTA Assembly | *rDNA variant units assembled by RIBOTIN* |

## Method 2: Homogenization Patterns of rDNA Arrays

Scripts focused on checking sequence homogeneity, calculating variant distribution depths, tracking gene conversion forces, and outputting visualization metrics.

### Scripts & Workflows

| Script / File | Type | Description |
| :--- | :--- | :--- |
| `rDNA_HiFiCoveragemapping.slurm` | SLURM Script | *Maps long HiFi reads to X and Y chromosome assemblies utilizing Winnowmap* |
| `rDNA_FFTalignment.slurm` | SLURM Script | *Aggregates clustered rDNA structural variants (morphs) from all strains and chromosomes into a unified dataset. Uses MAFFT to generate a master multiple sequence alignment required for downstream cross-strain homogenization analysis.* |
| `morphAnnotation.py` | Python Utility | *Analyzes assembled rDNA variants to generate structural summary plots of internal subunits like 18S, 28S, and internal transcribed spacers (ITS). Outputs IGV-ready color-coded BED tracks and abundance tables to map specific morph distributions across the genome.* |

### Figures and Array Visualizations

| Script / File | Type | Description |
| :--- | :--- | :--- |
| `rDNA_coveragefigures.py` | Python Utility | *Generates high-resolution architectural plots visualizing read depth across targeted rDNA arrays.* |
| `histograms.py` | Python Utility | *Generates log-scaled histograms to visualize sequence identity and structural homogenization within tandem arrays.* |
| `rDNA_heatmaps.py` | Python Utility | *Generates contiguous sequence identity heatmaps to visualize intra-array divergence among rDNA units.Explicitly color-codes structural variants to contrast canonical units against transposon-inserted elements across the X and Y chromosomes* |
| `rDNA_heatmaps_gapspenalized.py` | Python Utility | *Visualizes inter-array sequence homogenization by generating comparative heatmaps that strictly penalize alignment gaps. Highlights structural and sequence divergence when directly comparing specific rDNA populations across different sex chromosomes or genetic lineages.* |

---