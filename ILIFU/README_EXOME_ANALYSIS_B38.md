# Exome Analysis Pipeline - GRCh38

## Overview

This directory contains the configuration and run script for running the exome analysis pipeline on GRCh38-lifted VCF files.

## Files Created

### 1. Run Script
**Location:** `/users/mamana/exome_aibst_b38/ILIFU/exome_analysis_nextflow.sh`

Executable bash script that:
- Validates requirements (Nextflow, Singularity)
- Checks config and workflow files exist
- Runs the exome analysis pipeline with appropriate settings

### 2. Configuration File
**Location:** `/users/mamana/exome_aibst_b38/ILIFU/exome_analysis_nextflow.config`

Nextflow config for GRCh38 analysis with:
- Input VCFs: `/users/mamana/exome_aibst_b38/liftover_results/final/chr%s.hg38.vcf.gz`
- GRCh38 reference genome and resources
- All GATK bundle files for hg38
- Updated PGX datasets (GRCh38 coordinates)
- Proper singularity configuration for ILIFU

### 3. Backup
**Location:** `/users/mamana/exome_aibst_b38/ILIFU/exome_analysis_nextflow.config.backup_*`

Your previous GRCh37 config has been backed up with timestamp.

## Key Configuration Changes

### Input Data
```groovy
dataset_files = [
    ['AIBST', "${data_dir}/liftover_results/final/chr%s.hg38.vcf.gz",
     "${data_dir}/data/AIBST/aibst_all_samples.csv",
     "KNK,KNL,KNM,KNP,NGH,NGI,NGY,SAV,TZA,TZB,ZWD,ZWS"],
]
```

### Reference Genome
```groovy
ref_genome = "${reference_dir}/gatk_bundle/hg38/Homo_sapiens_assembly38.fasta"
```

### Chromosomes
```groovy
chromosomes = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
```

Note: VCF files use "chr" prefix (chr1, chr2, etc.) as produced by liftover.

### Annotation Control Flags
```groovy
use_cosmic = false  // Set to true to enable COSMIC annotations (requires GRCh38 version)
```

**Note:** COSMIC annotations are disabled by default since they are not discussed in the manuscript results. The COSMIC database file available locally is GRCh37 format. To enable COSMIC annotations:
1. Set `use_cosmic = true` in the config
2. Ensure a GRCh38 COSMIC file is available at the configured path, or
3. Use the provided `liftover_cosmic.sh` script to convert the GRCh37 version to GRCh38

## Usage

### Run the Pipeline

```bash
cd /users/mamana/exome_aibst_b38
bash ILIFU/exome_analysis_nextflow.sh
```

The script will:
1. Show configuration summary
2. Verify all requirements
3. Prompt for confirmation
4. Run the pipeline with `-resume` support

### Run Directly with Nextflow

```bash
cd /users/mamana/exome_aibst_b38

nextflow run main.nf \
    -c ILIFU/exome_analysis_nextflow.config \
    -profile singularity \
    -resume
```

## Output

Results will be written to:
- Main output: `/users/mamana/exome_aibst_b38/results/`
- Reports: 
  - `execution_report.html`
  - `execution_timeline.html`
  - `pipeline_dag.html`

## Profiles Available

- **singularity** (recommended): Uses Singularity containers
- **local**: Local execution without containers
- **slurm**: For SLURM cluster execution

## Testing

For quick testing with fewer chromosomes, edit the config:

```groovy
chromosomes = [21,22]  // Test with chr21 and chr22 only
```

## Notes

1. **VCF Files Ready**: All 22 chromosome VCFs are lifted and indexed in `liftover_results/final/`
2. **Gzipped Output**: All VCFs are bgzipped (.vcf.gz) with tabix indices (.tbi)
3. **GRCh38 Resources**: Config points to all necessary GRCh38 reference files
4. **Population Groups**: Configured for AIBST populations (12 groups)

## Troubleshooting

### Check VCF Files
```bash
ls -lh /users/mamana/exome_aibst_b38/liftover_results/final/
```

### Verify Config Syntax
```bash
nextflow config ILIFU/exome_analysis_nextflow.config
```

### Check Workflow File
```bash
ls -l /users/mamana/exome_aibst_b38/main.nf
```

## Dependencies

- Nextflow >= 22.10.1
- Singularity >= 3.x
- Container: `quay.io/mypandos/pgx_tools`
