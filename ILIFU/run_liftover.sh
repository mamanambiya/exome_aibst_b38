#!/bin/bash
#
# VCF Liftover Script - GRCh37 to GRCh38 Conversion
# Uses AfriGen-D/vcf-liftover workflow from GitHub
# Author: Mamana Mbiyavanga
# Date: 2025-10-31
#
# Description:
#   Converts AIBST VCF files from GRCh37/hg19 to GRCh38/hg38 using
#   the AfriGen-D Nextflow pipeline with CrossMap.
#
# Usage:
#   bash run_liftover.sh
#
# Requirements:
#   - Nextflow >= 22.10.1
#   - Singularity (for containers)
#   - Input VCF files must be bgzipped and indexed (.tbi)
#

set -euo pipefail

# ============================================================================
# Configuration
# ============================================================================

# Input VCF directory
INPUT_VCF_DIR="/cbio/users/mamana/analysis/exome-analysis/exome_aibst/data/AIBST/VCF/CHRS"
INPUT_VCF_PATTERN="aibst_all_samples_chrm*_clean_snpeff_dbsnp.annot.vcf.gz"

# Reference genomes
GRCH37_REF="/cbio/users/mamana/resources/reference-data/reference/gatk_bundle/human_g1k_v37_decoy.fasta"
GRCH38_REF="/cbio/users/mamana/resources/reference-data/reference/gatk_bundle/hg38/Homo_sapiens_assembly38.fasta"

# Liftover chain file
CHAIN_FILE="/cbio/users/mamana/resources/reference-data/reference/liftover/GRCh37/hg19ToHg38.over.chain.gz"

# Output configuration
WORK_DIR="$(pwd)"
OUTPUT_DIR="${WORK_DIR}/liftover_results"
SAMPLESHEET="${WORK_DIR}/vcf_samplesheet.csv"

# Workflow location (local clone)
WORKFLOW_DIR="${WORK_DIR}/vcf-liftover"

# Nextflow configuration
NXF_SINGULARITY_CACHEDIR="${HOME}/.singularity"
export NXF_SINGULARITY_CACHEDIR

# Disable Nextflow metrics that require 'ps' command
NXF_DISABLE_CHECK_LATEST=true
NXF_ANSI_LOG=false
export NXF_DISABLE_CHECK_LATEST NXF_ANSI_LOG

# ============================================================================
# Functions
# ============================================================================

log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"
}

error() {
    echo "[ERROR] $*" >&2
    exit 1
}

check_requirements() {
    log "Checking requirements..."

    # Load Java 22 (required for Nextflow)
    if command -v module &> /dev/null; then
        module load java/openjdk-22.0.1
        log "✓ Loaded Java $(java -version 2>&1 | head -1)"
    fi

    # Check Nextflow
    if ! command -v nextflow &> /dev/null; then
        error "Nextflow not found. Please install Nextflow >= 22.10.1"
    fi

    # Check Singularity
    if ! command -v singularity &> /dev/null; then
        error "Singularity not found. Please load Singularity module or install it"
    fi

    # Check reference files
    [[ -f "$GRCH37_REF" ]] || error "GRCh37 reference not found: $GRCH37_REF"
    [[ -f "$GRCH38_REF" ]] || error "GRCh38 reference not found: $GRCH38_REF"
    [[ -f "$CHAIN_FILE" ]] || error "Chain file not found: $CHAIN_FILE"

    # Check input VCF directory
    [[ -d "$INPUT_VCF_DIR" ]] || error "Input VCF directory not found: $INPUT_VCF_DIR"

    # Check workflow directory
    [[ -d "$WORKFLOW_DIR" ]] || error "Workflow directory not found: $WORKFLOW_DIR (run: git clone https://github.com/AfriGen-D/vcf-liftover.git)"

    log "✓ All requirements satisfied"
}

generate_samplesheet() {
    log "Generating samplesheet: $SAMPLESHEET"

    # Create CSV header with correct column names
    echo "sample_id,vcf_path" > "$SAMPLESHEET"

    # Find all matching VCF files and add to samplesheet
    local count=0
    for vcf in ${INPUT_VCF_DIR}/${INPUT_VCF_PATTERN}; do
        if [[ -f "$vcf" ]]; then
            # Extract chromosome number for sample name
            local basename=$(basename "$vcf")
            local chr=$(echo "$basename" | sed -n 's/.*chrm\([0-9XY]*\)_.*/chr\1/p')

            # Check if .tbi index exists
            if [[ ! -f "${vcf}.tbi" ]]; then
                log "WARNING: Index file missing for $vcf - skipping"
                continue
            fi

            echo "${chr},${vcf}" >> "$SAMPLESHEET"
            count=$((count + 1))
        fi
    done

    if [[ $count -eq 0 ]]; then
        error "No VCF files found matching pattern: ${INPUT_VCF_DIR}/${INPUT_VCF_PATTERN}"
    fi

    log "✓ Generated samplesheet with $count VCF files"
    log "  Location: $SAMPLESHEET"
}

check_vcf_indices() {
    log "Checking VCF index files..."

    local missing_count=0
    while IFS=, read -r sample vcf; do
        [[ "$sample" == "sample_id" ]] && continue  # Skip header

        if [[ ! -f "${vcf}.tbi" ]]; then
            log "  Missing index: ${vcf}.tbi"
            missing_count=$((missing_count + 1))
        fi
    done < "$SAMPLESHEET"

    if [[ $missing_count -gt 0 ]]; then
        error "$missing_count VCF files are missing .tbi index files. Please index them with: tabix -p vcf <file>.vcf.gz"
    fi

    log "✓ All VCF files are properly indexed"
}

run_liftover() {
    log "Starting VCF liftover workflow..."
    log "  Workflow: $WORKFLOW_DIR"
    log "  Input: $SAMPLESHEET"
    log "  Target: GRCh38"
    log "  Output: $OUTPUT_DIR"

    nextflow run "$WORKFLOW_DIR" \
        --input "$SAMPLESHEET" \
        --target_fasta "$GRCH38_REF" \
        --chain_file "$CHAIN_FILE" \
        --source_build hg19 \
        --target_build hg38 \
        --outdir "$OUTPUT_DIR" \
        --validate_output true \
        --check_build_compatibility false \
        --singularity_cache_dir "$NXF_SINGULARITY_CACHEDIR" \
        -profile singularity \
        -resume

    local exit_code=$?

    if [[ $exit_code -eq 0 ]]; then
        log "✓ Liftover workflow completed successfully!"
        log ""
        log "Results location: $OUTPUT_DIR"
        log "  - Lifted VCFs: ${OUTPUT_DIR}/lifted_vcf/"
        log "  - Unmapped variants: ${OUTPUT_DIR}/unmapped/"
        log "  - QC reports: ${OUTPUT_DIR}/qc/"
        log "  - Execution report: ${OUTPUT_DIR}/execution_report.html"
    else
        error "Liftover workflow failed with exit code: $exit_code"
    fi
}

show_summary() {
    log ""
    log "=========================================="
    log "VCF Liftover Summary"
    log "=========================================="
    log "Input VCFs: $(grep -c '^chr' "$SAMPLESHEET" || echo 0)"
    log "Source Build: GRCh37/hg19"
    log "Target Build: GRCh38/hg38"
    log "Output Directory: $OUTPUT_DIR"
    log ""
    log "Next Steps:"
    log "1. Review execution report: ${OUTPUT_DIR}/execution_report.html"
    log "2. Check QC statistics in: ${OUTPUT_DIR}/qc/"
    log "3. Validate lifted VCFs in: ${OUTPUT_DIR}/lifted_vcf/"
    log "4. Update workflow configs to use lifted VCFs"
    log "=========================================="
}

# ============================================================================
# Main Execution
# ============================================================================

main() {
    log "=========================================="
    log "VCF Liftover - GRCh37 to GRCh38"
    log "=========================================="
    log "Working directory: $WORK_DIR"
    log ""

    check_requirements
    generate_samplesheet
    check_vcf_indices

    log ""
    read -p "Proceed with liftover? [y/N] " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        log "Aborted by user"
        exit 0
    fi

    run_liftover
    show_summary
}

# Run main function
main "$@"
