#!/bin/bash
#
# Exome Analysis Pipeline - GRCh38
# Uses lifted VCF files from GRCh37 to GRCh38 conversion
# Author: Mamana Mbiyavanga
# Date: 2025-11-06
# Updated: 2025-11-14 - Fixed dataset_qc_alt and annotate_snpeff processes
#

set -euo pipefail

# ============================================================================
# Configuration
# ============================================================================

# Directories
WORK_DIR="/cbio/users/mamana/exome_aibst_b38"
CONFIG="${WORK_DIR}/HPC/exome_analysis_nextflow.config"
WORKFLOW="${WORK_DIR}/main.nf"
OUTPUT_DIR="${WORK_DIR}/results"
LOG_DIR="${WORK_DIR}/logs"

# Nextflow configuration
NXF_SINGULARITY_CACHEDIR="/cbio/users/mamana/singularity-containers"
export NXF_SINGULARITY_CACHEDIR

# Create directories
mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

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

    # Check config file
    [[ -f "$CONFIG" ]] || error "Config file not found: $CONFIG"

    # Check workflow file
    [[ -f "$WORKFLOW" ]] || error "Workflow file not found: $WORKFLOW"

    log "✓ All requirements satisfied"
}

show_info() {
    log ""
    log "================================================================================"
    log "  Exome Analysis Pipeline - GRCh38 (OPTIMIZED)"
    log "================================================================================"
    log "Working directory: $WORK_DIR"
    log "Config: $CONFIG"
    log "Workflow: $WORKFLOW"
    log "Output: $OUTPUT_DIR"
    log "Logs: $LOG_DIR"
    log "Singularity cache: $NXF_SINGULARITY_CACHEDIR"
    log ""
    log "Recent Fixes Applied:"
    log "  ✓ dataset_qc_alt: Fixed FORMAT field errors (AD, PL, DP)"
    log "  ✓ annotate_snpeff: Optimized memory (18GB) and CPU (4 cores)"
    log "  ✓ annotate_dbsnp: Increased parallelization (maxForks: 20)"
    log "  ✓ All processes: Added automatic VCF indexing"
    log "================================================================================"
    log ""
}

run_pipeline() {
    log "Starting exome analysis pipeline..."
    log "Profiles: slurm,singularity"
    log "Resume: enabled (-resume flag)"
    log ""

    cd "$WORK_DIR"

    # Generate timestamp for logs
    TIMESTAMP=$(date '+%Y%m%d_%H%M%S')

    nextflow run "$WORKFLOW" \
        -c "$CONFIG" \
        -profile slurm,singularity \
        -resume \
        -with-report "${LOG_DIR}/execution_report_${TIMESTAMP}.html" \
        -with-timeline "${LOG_DIR}/execution_timeline_${TIMESTAMP}.html" \
        -with-dag "${LOG_DIR}/pipeline_dag_${TIMESTAMP}.html"

    local exit_code=$?

    log ""
    log "================================================================================"
    if [[ $exit_code -eq 0 ]]; then
        log "  ✓ PIPELINE COMPLETED SUCCESSFULLY"
        log "================================================================================"
        log ""
        log "Results location: $OUTPUT_DIR"
        log "  - Annotated VCFs: ${OUTPUT_DIR}/data/AIBST/VCF/"
        log "  - SnpEff statistics: ${OUTPUT_DIR}/snpeff_stats/"
        log ""
        log "Execution reports:"
        log "  - Report: ${LOG_DIR}/execution_report_${TIMESTAMP}.html"
        log "  - Timeline: ${LOG_DIR}/execution_timeline_${TIMESTAMP}.html"
        log "  - DAG: ${LOG_DIR}/pipeline_dag_${TIMESTAMP}.html"
        log ""
        log "Monitor performance:"
        log "  ./scripts/monitor_annotation_pipeline.sh"
    else
        log "  ✗ PIPELINE FAILED (Exit code: $exit_code)"
        log "================================================================================"
        log ""
        log "Check logs:"
        log "  - Nextflow log: .nextflow.log"
        log "  - Latest run: ${LOG_DIR}/execution_report_${TIMESTAMP}.html"
        error "Pipeline failed with exit code: $exit_code"
    fi
}

# ============================================================================
# Main Execution
# ============================================================================

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Exome Analysis Pipeline - GRCh38 (Optimized)

OPTIONS:
    -h, --help          Show this help message

EXAMPLES:
    # Run directly
    $0

    # Submit to SLURM cluster
    sbatch $0

EOF
    exit 0
}

main() {
    # Parse arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                usage
                ;;
            *)
                echo "Unknown option: $1"
                usage
                ;;
        esac
    done

    show_info
    check_requirements

    run_pipeline
}

# Run main function
main "$@"
