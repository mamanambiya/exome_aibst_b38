# African Exome Analysis Workflow

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20639159.svg)](https://doi.org/10.5281/zenodo.20639159)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A comprehensive Nextflow DSL2 pipeline for analysing whole-exome sequencing
data from African populations, with a focus on pharmacogenomic (PGx) variants,
population structure, and highly differentiated variants (HDV).

This repository contains the analysis pipeline and scripts underlying the
manuscript *"Pharmacogenomic diversity across 12 African populations"* (Mbiyavanga
et al.). Whole-exome sequencing data are deposited in the European
Genome-phenome Archive (EGA) under accession **EGAS00001008456** (controlled
access).

## Overview

The workflow processes joint-called VCFs through quality control, annotation,
population-level analysis, and reporting:

- Variant annotation (SnpEff, dbSNP, ClinVar, COSMIC, CADD, AlphaMissense, dbNSFP)
- Population allele frequencies
- Pharmacogenomic (PGx) star-allele and phenotype calling (PyPGx)
- HLA typing (OptiType / arcasHLA)
- Highly differentiated variants (HDV) between populations
- F_ST analysis for population differentiation (Weir & Cockerham)
- Principal Component Analysis (PCA, EIGENSOFT smartpca)
- ADMIXTURE ancestry estimation
- Novel-variant discovery (absent from dbSNP, gnomAD, AGVP)

## Quick Start

### Prerequisites

- Nextflow (≥ 22.10)
- Singularity (or Docker)
- A SLURM HPC cluster (configs provided for ILIFU, HPC, CHPC)

### Run

```bash
nextflow run main_2025.nf \
    -profile singularity,slurm \
    -c ILIFU/exome_analysis_nextflow.config
```

The pipeline entry point is `main_2025.nf` (DSL2). Cluster-specific submission
scripts live under `ILIFU/`, `HPC/`, and `CHPC/`.

## Repository Layout

```text
.
├── main_2025.nf            # Pipeline entry point (DSL2)
├── module/                 # Process modules
│   ├── vcf_qc.nf           # Quality control
│   ├── annotate_vcf.nf     # SnpEff / dbSNP / ClinVar / COSMIC annotation
│   ├── pathogenicity_scores.nf  # CADD, AlphaMissense, dbNSFP
│   ├── freq_vcf.nf         # Allele-frequency calculation
│   ├── fst.nf              # F_ST analysis
│   ├── hdv_vcf.nf          # HDV detection
│   ├── hdv_fst_plots.nf    # HDV / F_ST plotting
│   ├── novel_vcf.nf        # Novel-variant discovery
│   ├── csq.nf              # Consequence analysis
│   ├── pgx_calling.nf      # PyPGx star-allele calling
│   ├── hla_typing.nf       # HLA typing
│   ├── structure.nf        # PCA (EIGENSOFT smartpca)
│   ├── structure_advanced.nf    # ADMIXTURE ancestry estimation
│   ├── chrx_ploidy.nf      # chrX sex inference + ploidy correction
│   └── subset_vcf.nf       # VCF subsetting
├── templates/              # Process scripts (Python / R)
├── bin/                    # Helper scripts
├── ILIFU/                  # ILIFU cluster configs + submission scripts
├── docs/                   # Workflow architecture and graph documentation
└── supp/                   # Supplementary tables published with the manuscript
```

## Configuration

Pipeline parameters (input datasets, reference paths, population groups,
F_ST cutoffs, chunk sizes) are set in the cluster config files under `ILIFU/`.
The DeepVariant WES +2 kb-padded GRCh38 track used for the manuscript is
configured in `ILIFU/dv_wes_padded_2kb.config`. Absolute paths in these configs
are specific to the authors' ILIFU environment and must be adapted for other
sites.

## Documentation

- `docs/WORKFLOW_GRAPH.md` — Mermaid graph of the full workflow
- `docs/WORKFLOW_ARCHITECTURE.md` — technical architecture
- `docs/ANNOTATION_PIPELINE_WORKFLOW.md` — annotation chain detail

## Software (via containers)

BCFtools, VCFtools, PLINK, SnpEff/SnpSift, EIGENSOFT, ADMIXTURE, PyPGx,
OptiType, arcasHLA, Python 3 (pandas, numpy, scipy), R (ggplot2, tidyverse).

## Citation

If you use this pipeline, please cite both the software and the manuscript:

> **Software:** Mbiyavanga M, et al. *African Exome Analysis Workflow.*
> Zenodo. <https://doi.org/10.5281/zenodo.20639159>
>
> **Manuscript:** Mbiyavanga M, et al. *Pharmacogenomic diversity across 12
> African populations.* (manuscript in preparation).

## Contact

- **Author:** Mamana Mbiyavanga
- **Institution:** Computational Biology Division, University of Cape Town
- **Email:** <mamana.mbiyavanga@uct.ac.za>

## License

Released under the [MIT License](LICENSE).

## Acknowledgments

- H3Africa and the AIBST consortium
- ILIFU High Performance Computing facility
