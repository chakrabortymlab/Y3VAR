# Gene and Repeat Annotation

[![Status: In Progress](https://img.shields.io/badge/status-active-success.svg)](#)
[![Pipeline: Annotation](https://img.shields.io/badge/pipeline-gene__annotation-orange.svg)](#)

This directory contains the scripts, reference datasets, and workflows used to manually curate and annotate the 13 major Y-linked fertility factor genes and repetitive elements. 

The directory is split into three core methods: **Manual Curation of Y-linked Loci**, **Direct RNA sequencing (dRNA-seq) and Transcriptomic Validation**, and **Repeat Annotation**.

---

## Method 1: Manual Curation of Y-linked Loci

Given the high repetitive content and duplications and sequencing bias associated with these regions, specific loci required targeted assembly and manual curation using known fertility factors and specialized scripts.

### Scripts & Reference Files

| Script / File | Type | Description |
| :--- | :--- | :--- |
| `FertilityGenesCDS_Carvalho_et_al.fasta` | FASTA Reference | *Coding sequences of known Y-linked fertility factors and marker genes, obtained from Carvalho et al., 2026 . Serves as the baseline query reference for gene mapping, annotation, and reconstruction workflows.* |
| `Pp1Y1_targetedassembly.slurm` | SLURM Script | *Slurm pipeline that uses Minimap2 and Seqtk to extract long reads specific to the Pp1-Y1 locus. Runs Flye on the isolated sequence pool to generate a targeted, de novo assembly.* |
| `wholecontigscuration.py` | Python Utility | *Executes multi-step sequence integrations to patch missing gene fragments into the primary assembly. This script specifically integrates the Illumina-based targeted assemblies performed by the Carvalho et al., 2026 study* |
| `kl3andkl5curation.py` | Python Utility | *Performs targeted sequence integration on specific Y chromosome assembly scaffolds using Biopython. Directly patches known sequence gaps inside the massive kl-3 and kl-5 fertility genes by inserting missing exon fragments. Sequence was carried over from the Chang & Larracuente 2019 assembly* |
| `PprYcuration.py` | Python Utility | *Patches and corrects exon duplications within the Ppr-Y gene locus by parsing structural alignment variants from a SAM file and executing targeted insertions, deletions, or replacements.* |
| `PRYcuration.py` | Python Utility | *Integrates missing PRY coding sequence into specific scaffold coordinates.* |
| `WDYcuration.py` | Python Utility | *Resolves structural assemblies near the WDY gene locus by relocating exonic sequence to its syntenic position* |
| `removePRY.py` | Python Utility | *Removes a duplicated portion of the PRY gene resulting from a misassembly from a target scaffold and updates the corresponding annotation file.* |

## Method 2: Direct RNA Sequencing & Transcriptomic Validation

To validate the manually curated structures and generate final gene models, raw direct RNA-seq data was mapped back to the assembly. Coding sequences were reconstructed through a custom pipeline and exon coverages were calculated to confirm male-linkage.

### Scripts & Workflows

| Script / File | Type | Description |
| :--- | :--- | :--- |
| `dRNAseqMapping.slurm` | SLURM Script | *Uses Minimap2 to map direct long-read RNA-seq data to the curated ISO-1 assembly.* |
| `BLAT.slurm` | SLURM Script | *Slurm script that maps reference coding sequences (CDS) against multiple Y chromosome assemblies using BLAT. Optimized to successfully span massive Y-linked introns* |
| `ExonCoverages.py` | Python Utility | *Extracts coordinates from a reference GFF3 file and quantifies average read depth across multiple BAM datasets, including male and female samples. Generates a TSV support matrix evaluating individual exon validation thresholds.* |
| `CDSreconstruction.py` | Python Utility | *Reconstructs target gene coding sequences from genomic alignments using a reference sequence template. Outputs the reconstructed FASTA files alongside a detailed Excel spreadsheet auditing sequence completeness and documenting multi-mapping* |
| `GFFannotations.py` | Python Utility | *Converts raw PSL BLAT alignment files into GFF3 annotations. Color-codes genes for viewing on IGV* |
| `AnnotationProduction.py` | Python Utility | *Parses PSL alignment files from BLAT to stitch fragmented hits into comprehensive GFF3 gene structures.* |

## Method 3: Repeat Annotation

Handles the identification, classification, and comparative masking of transposable elements (TEs) and satellite repeats across assemblies.

### Scripts & Workflows

| Script / File | Type | Description |
| :--- | :--- | :--- |
| `ComparativeRepeats.py` | Python Utility | *Aggregates transposable element and satellite DNA abundances across multiple genome assemblies from RepeatMasker outputs. Organizes the compiled dataset into a multi-sheet Excel file highlighting repeat classes and highly repetitive contigs.* |

---