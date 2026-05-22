# Crystal Stellate Analysis

[![Status: Active](https://img.shields.io/badge/status-active-success.svg)](#)
[![Analysis: Crystal-Stellate](https://img.shields.io/badge/analysis-array__evolution-purple.svg)](#)

This directory contains the scripts, SLURM cluster workflows, references, and visualization utilities used to look at the structural organization, copy number dynamics, transcriptional activity, and evolutionary history of the X-linked *Stellate* and Y-linked *Su(Ste)* and *PCKR* tandem arrays.

---

## 🗺️ Method 1: Structural Annotation of PCKR/Su(Ste) Array

Primary identification, boundary definition, and target locus verification of the repetitive *PCKR* and *Su(Ste)* units across the assembled sequence.

### 📋 Scripts & Reference Files

| Script / File | Type | Description |
| :--- | :--- | :--- |
| `PCKR_SuSte_All_Units.fasta` | 🧬 FASTA Reference | *[All individually extracted PCKR and Su(Ste) repeat units spanning the target assemblies. Serves as the primary sequence input for downstream multiple sequence alignments and phylogenetic tree construction.]* |
| `RepSequences.fasta` | 🧬 FASTA Reference | *[Standard consensus sequences for Stellate family genes.]* |
| `FlybaseGenes.slurm` | 🚀 SLURM Script | *[Parses BLAST alignments for βNACtes1 and Ssl. It automatically re-formats these filtered hits into standard BED track format for viewing on IGV.]* |
| `FlybaseTranscripts.slurm` | 🚀 SLURM Script | *[BLAST alignments of the FlyBase transcript library (dmel-all-transcript-r6.67.fasta).]* |
| `AnchorFinding.slurm` | 🚀 SLURM Script | *[Script that locates structural anchors between assemblies. Isolates negative-strand (inverted) repeat coordinates from BED annotations and exports isolated homologous sequences into grouped multi-FASTA files for alignment.]* |
| `AnchorEvaluation.py` | 🐍 Python Utility | *[Parses the multi-FASTA alignments generated from the spatial groups and calculates pairwise sequence identities while safely ignoring double-gaps. If a trio maintains an average sequence identity of 95%, it is recorded as a highly conserved structural anchor.]* |

---

## 🎙️ Method 2: Male-Specific Transcription of PCKR/Su(Ste)

Workflows dedicated to assessing the sex-specific expression dynamics of the Y-linked suppressors.

### 📋 Scripts & Workflows

| Script / File | Type | Description |
| :--- | :--- | :--- |
| `MaleSpecificConsensusSequence.slurm` | 🚀 SLURM Script | *[Maps male and female long-read RNA-seq (dRNA-seq) data directly to isolated repeat consensus sequences. Generates a Python script to profile expression coverage and validate male-specific transcription patterns.]* |
| `StellateMapping.slurm` | 🚀 SLURM Script | *[Maps Su(Ste) and PCKR elements and merges their genomic footprints using BEDtools.]* |
| `BNACtes1BLAST.slurm` | 🚀 SLURM Script | *[BLAST alignment for the location of BNACtes1 promoter sequences. βNACtes1 (FBgn0030538) and Ssl (FBgn0015300) as queries.]* |

## 🌿 Method 3: Comparative and Phylogenetic Analysis

Scripts designed to assess the evolutionary relationships, sequence variance, and insertion/deletion dynamics between different individual repeat units across strains or assemblies.

### 📋 Scripts & Workflows

| Script / File | Type | Description |
| :--- | :--- | :--- |
| `phylo.slurm` | 🚀 SLURM Script | *[Constructs Maximum Likelihood phylogenetic trees from the aligned repeat sequences using IQ-TREE.]* |
| `PhyloTreeFigure.py` | 🐍 Python Utility | *[Visualizes the generated phylogenetic trees using Biopython and Matplotlib.]* |
| `PairwiseAlignment_ystellate.slurm` | 🚀 SLURM Script | *[Aligns compiled Stellate-related repeat units using MAFFT L-INS-i algorithm]* |
| `InDel_Variation.py` | 🐍 Python Utility | *[Parses multi-strain sequence alignments to build a global consensus backbone and calculates positional variation metrics, including deletions, insertions, and divergence. Generates resultant figures]* |

---

## ❌ Method 4: X-Linked Stellate Consensus

Analyses targeting the corresponding X-chromosomal *Stellate* loci to investigate gene homogenization patterns and sequence features across the array.

### 📋 Scripts & Workflows

| Script / File | Type | Description |
| :--- | :--- | :--- |
| `XlinkedStellateMapping.slurm` | 🚀 SLURM Script | *[Maps Euchromatic and Heterochromatic Stellate consensus sequences to X chromosome assemblies using BLASTn.]* |
| `XlinkedPairwiseAlignment.slurm` | 🚀 SLURM Script | *[Aligns all structurally resolved X-linked Stellate repeat units using MAFFT's L-INS-i algorithm.]* |
| `XlinkedHomogenization.py` | 🐍 Python Utility | *[Quantifies the degree of sequence homogenization within and across X-linked Stellate arrays by parsing alignment matrices.]* |

## 🗺️ Method 5: Hi-C Interactions

Chromosome Conformation Capture (Hi-C) pipelines used to determine physical interactions, spatial positioning, and 3D boundaries surrounding the repeat arrays.

### 📋 Scripts & Documentation

| Script / File | Type | Description |
| :--- | :--- | :--- |
| `HiCreadTrimming.slurm` | 🚀 SLURM Script | *[Quality trimming of raw paired-end Illumina Hi-C reads.]* |
| `HiCgenerate_site_positions.py` | 🐍 Python Utility | *[A core structural utility (aligned with the Juicer processing suite) that parses reference genome FASTA files to locate every single instance of a specific restriction enzyme recognition motif.]* |
| `Juicer.slurm` | 🚀 SLURM Script | *[Executes the complete Juicer Hi-C pipeline on a custom genome assembly.]* |
| `How2HIC.txt` | 📜 Reference / Notes | *[A command-line reference guide for downstream Hi-C analysis, detailing how to convert Juicer outputs to multi-resolution .mcool formats. Provides copy-paste examples for generating whole-genome contact maps, calling TADs, and plotting annotations using HiCExplorer.]* |
| `How2Juicer.txt` | 📜 Reference / Notes | *[Reference documentation for preprocessing raw paired-end Hi-C reads with fastp and running the core Juicer pipeline. Outlines the necessary steps to index the reference genome, generate restriction site profiles, and execute the contact matrix assembly.]* |

---

## 🔢 Method 6: Calculation of Similarity Between Array Units

Custom identity matrices and geographic heatmaps utilized to visualize structural clustering patterns and track sequence drift over spatial distance across the arrays.

### 📋 Scripts & Workflows

| Script / File | Type | Description |
| :--- | :--- | :--- |
| `AdjacentUnitsIdentityPlots.py` | 🐍 Python Utility | *[Tracks structural variation along the array by plotting sequence identity drops specifically between immediately adjacent repeat units.]* |
| `XlinkedAdjacentIdentityPlot.py` | 🐍 Python Utility | *[Tracks and visualizes structural sequence variation strictly between physically adjacent Stellate repeat units along the X chromosome]* |
| `Heatmaps.py` | 🐍 Python Utility | *[Computes all-versus-all nucleotide identity matrices from combined sequence sets.]* |
| `Xlinkedheatmaps.py` | 🐍 Python Utility | *[Generates high-contrast, triangular heatmaps to visualize the all-versus-all sequence identity of X-linked Stellate arrays]* |
| `group_spatial_fastas.py` | 🐍 Python Utility | *[Spatial grouping of repeat domains for utilization in the previous HiC scripts.]* |

## 📊 Method 7: Su(Ste) and PCKR Copy Number and Array Size

Workflows designed to calculate assembly collapse and true array size

### 📋 Scripts & Workflows

| Script / File | Type | Description |
| :--- | :--- | :--- |
| `HiFiandIllumina_ArrayCollapse.py` | 🐍 Python Utility | *[Generates robust coverage plots from HiFi and Illumina read depths while applying a RepeatMasker-derived transposable element (TE) mask to filter out noise.]* |
| `HiFiReadCoverageUnmasked.slurm` | 🚀 SLURM Script | *[A Slurm pipeline that maps PacBio HiFi reads against unmasked sequence arrays using Winnowmap2 and Meryl k-mer databases]* |
| `IlluminaReadCoverageMasked.slurm` | 🚀 SLURM Script | *[A Slurm batch script that maps Illumina short reads to hard-masked array references using BWA-MEM.]* |
| `HiFiReadDepth.slurm` | 🚀 SLURM Script | *[Calculates per-base and 10kb rolling window read depths from aligned HiFi BAM files using SAMtools. Computes comprehensive summary statistics—including mean, median, and breadth of coverage—to assess sequence assembly quality across repetitive regions.]* |
| `YlinkedvalidationofCG_contig.slurm` | 🚀 SLURM Script | *[Validates the true Y-linkage of a targeted sequence contig by strictly mapping both male and female HiFi reads using Minimap2. Computes per-base sequence coverage to confirm strict male-specific read pileups.]* |
| `MannWhitneyTest.py` | 🐍 Python Utility | *[Performs Mann-Whitney U statistical tests to compare sequence identity distributions between intra-strain and inter-strain array structures.]* |

---

## 🎨 Figure Generation Scripts

These scripts generate a few of our supplementary figures

| Script / File | Type | Description |
| :--- | :--- | :--- |
| `supplementaryfigure5.py` | 🐍 Python Utility | *[Generates structural architecture diagrams for the PCKR repeat unit, rendering integrated elements like HeT-A and Copia]* |
| `supplementaryfigure6.py` | 🐍 Python Utility | *[Visualizes normalized, log-scaled transcriptional read coverage across specific genomic contigs using both male and female depth data.]* |
| `'supplementaryfigures7&8.py'` | 🐍 Python Utility | *[Generates diagrams of the structural architecture and transcriptional activity of the PCKR and Su(Ste) repeat arrays]* |

---