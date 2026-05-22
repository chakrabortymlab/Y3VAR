# Maintenance of Functional Repeat Arrays on the *Drosophila* Y Chromosome

[![Static Badge](https://img.shields.io/badge/Release-v1.0.0-blue)](https://github.com/)
[![Static Badge](https://img.shields.io/badge/Data_Type-Genomics%20%26%20Transcriptomics-darkgreen)](https://github.com/)
[![Static Badge](https://img.shields.io/badge/Organism-Drosophila%20melanogaster-red)](https://github.com/)

This repository contains the custom scripts, SLURM cluster workflows, and documentation utilized for the study **"Maintenance of functional repeat arrays on the *Drosophila* Y chromosome"**.
The codebase is organized into modular directories corresponding directly to the major phases of our methods section. Each directory contains its own dedicated `README.md` detailing script executions, parameters, and cluster environment requirements.

---

## Repository Structure & Core Pipelines

Navigate to the individual directories below for specific step-by-step documentation and file listings:

### [1. Assembly and Validation]
* **Focus:** Assembly workflows, scaffolding scripts, and long-read alignment validations.
* **Key Data:** Evaluates quality and continuity across our three distinct Y chromosome assemblies.

### [2. Gene and Repeat Annotation]
* **Focus:** Specialized annotations and manual curation for repetitive and previously poorly characterized genomic landscapes.
* **Key Data:** Manual curation of the 13 major Y-linked fertility genes in our ISO-1 assembly, transcriptomic validation using A4 male and female Nanopore direct RNA sequencing (dRNA-seq) data, and annotation of repeats.

### [3. Crystal Stellate Analysis]
* **Focus:** Structural annotation, copy number estimations, and evolution of X and Y-linked Stellate tandem repeats
* **Key Data:** Scripts for the structural variation, spatial sequence identity, and 3D Hi-C interactions of the *PCKR* and *Su(Ste)* tandem suppressor arrays.

### [4. rDNA Assembly and Characterization]
* **Focus:** Annotation and targeted assembly of ribosomal RNA repeat units.
* **Key Data:** Scripts for the isolation, variant characterization, and array homogenization analysis of both X- and Y-linked rDNA tandem arrays.

---

## System & Environment Requirements

* **Platform:** The majority of computational steps (especially `.slurm` scripts) are optimized to run on High-Performance Research Computing (HPRC) Linux clusters using the SLURM workload manager.
* **Dependencies:** Environment requirements (such as Conda environments, Python scripts, or specific module names and versions) are detailed within their respective sub-directories.

---

## ✍Citation & Code Availability

If you use the scripts or resources in this repository for your research, please cite our publication:

> *Adolfo A. Delgado, Alejandra Samano, Mahul Chakraborty. (2026). Maintenance of functional repeat arrays on the Drosophila Y chromosome. [Journal Name, Volume, Pages, DOI placeholder]*

```bibtex
@article{YourCitation2026,
  author    = {Author, Anonymized and Author, Peer},
  title     = {Maintenance of functional repeat arrays on the Drosophila Y chromosome},
  journal   = {Journal Placeholder},
  year      = {2026},
  doi       = {10.XXXX/placeholder}
}