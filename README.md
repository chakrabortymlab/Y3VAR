
This repository contains the custom scripts, SLURM cluster workflows, and documentation utilized for the study **"Higher-order Architecture Shapes Concerted Evolution in a Y-linked repeat array"**.
The codebase is organized into modular directories corresponding directly to the major phases of our methods section. Each directory contains its own dedicated `README.md` detailing script executions, parameters, and cluster environment requirements.

---

## Repository Structure & Core Pipelines

Navigate to the individual directories below for specific step-by-step documentation and file listings:

### 1. Assembly and Validation
* **Focus:** Assembly workflows, scaffolding scripts, and long-read alignment validations.
* **Key Data:** Evaluates quality and continuity across our three distinct Y chromosome assemblies.

### 2. Gene and Repeat Annotation
* **Focus:** Specialized annotations and manual curation for repetitive and previously poorly characterized genomic landscapes.
* **Key Data:** Manual curation of the 13 major Y-linked fertility genes in our ISO-1 assembly, transcriptomic validation using A4 male and female Nanopore direct RNA sequencing (dRNA-seq) data, and annotation of repeats.

### 3. Crystal Stellate Analysis
* **Focus:** Structural annotation, copy number estimations, and evolution of X and Y-linked Stellate tandem repeats
* **Key Data:** Scripts for the structural variation, spatial sequence identity, and 3D Hi-C interactions of the *PCKR* and *Su(Ste)* tandem suppressor arrays.

### 4. rDNA Assembly and Characterization
* **Focus:** Annotation and targeted assembly of ribosomal RNA repeat units.
* **Key Data:** Scripts for the isolation, variant characterization, and array homogenization analysis of both X- and Y-linked rDNA tandem arrays.

---

## System & Environment Requirements

* **Platform:** The majority of computational steps (especially `.slurm` scripts) are optimized to run on High-Performance Research Computing (HPRC) Linux clusters using the SLURM workload manager.
* **Dependencies:** Environment requirements (such as Conda environments, Python scripts, or specific module names and versions) are detailed within their respective sub-directories. May vary for individuals. 

---

## Citation

Coming Soon
