# Workflow Architecture

This document provides a detailed technical overview of the exome analysis workflow architecture, data flow, and module organization.

## Table of Contents
- [Overview](#overview)
- [Workflow DAG](#workflow-dag)
- [Module Organization](#module-organization)
- [Data Flow](#data-flow)
- [Process Details](#process-details)
- [Channel Operations](#channel-operations)
- [Error Handling](#error-handling)

## Overview

The workflow is built using **Nextflow DSL2**, which provides:
- Modular process definitions
- Explicit channel operations
- Workflow composition
- Better code reusability

### Key Design Principles

1. **Modularity:** Processes are organized into domain-specific modules
2. **Scalability:** Chunked processing for large datasets
3. **Parallelization:** Chromosome-wise and chunk-wise parallel execution
4. **Fault Tolerance:** Automatic retry with resource scaling
5. **Reproducibility:** Container-based execution environments

## Workflow DAG

### High-Level Workflow Structure

```
main workflow
├── data_proc                      # Core data processing
│   ├── check_files               # Validate input files
│   ├── get_map                   # Extract chromosome info
│   ├── generate_chunks_vcf       # Split into chunks
│   ├── split_vcf_to_chunk        # Extract chunk regions
│   ├── dataset_qc_alt            # Quality control
│   ├── annotate_snpeff           # Functional annotation
│   ├── annotate_dbsnp            # dbSNP annotation
│   ├── annotate_clinvar          # ClinVar annotation
│   ├── annotate_cosmic           # COSMIC annotation
│   ├── annotate_vcf_from_vcf     # MAF annotation
│   ├── fill_tags                 # Fill INFO tags
│   ├── concat_chunks_vcf         # Merge chunks
│   ├── singleton_dataset_chrm    # Identify singletons
│   ├── get_gene_vcf              # Extract PGx variants
│   ├── concat_vcf                # Concatenate chromosomes
│   └── merge_datasets_group      # Merge across datasets
│
├── data_proc_pop                  # Population-level processing
│   ├── split_dataset_vcf_pop     # Split by population
│   ├── get_pgx_pop               # Extract PGx per pop
│   └── concat_pop                # Concatenate per pop
│
├── group_pops_by_chrm            # Group populations by chr
│   └── merge_pop_groups          # Merge population VCFs
│
├── group_pops                     # Group populations (all chr)
│   └── cat_chrm_groups           # Concatenate chromosomes
│
├── group_2pops                    # Create pairwise groups
│   └── merge_pop_groups          # Merge 2-pop VCFs
│
├── group_pops_pgx                # Group PGx data
│   └── merge_pop_groups          # Merge PGx VCFs
│
├── group_2pops_pgx               # Pairwise PGx groups
│   └── merge_pop_groups          # Merge 2-pop PGx VCFs
│
├── fisher_group_pgx              # Fisher's exact test
│   ├── vcf_to_plink              # Convert to PLINK
│   ├── generate_pheno_fam_sample # Generate phenotype file
│   └── fisher_test_plink         # Run Fisher test
│
├── fisher_2pops_pgx              # Pairwise Fisher test
│   ├── vcf_to_plink              # Convert to PLINK
│   ├── generate_pheno_fam_sample1# Generate phenotype
│   ├── fisher_test_plink_pops    # Fisher test
│   ├── filter_fisher_test        # Filter results
│   └── combine_csv               # Combine results
│
├── fst                           # FST analysis
│   ├── fst_analysis              # Calculate FST
│   ├── get_fst_weir_estimates    # Extract estimates
│   ├── combine_fst_weir_estimates# Combine all
│   ├── get_fst_weir_estimates_cutoff # Filter by cutoff
│   ├── combine_fst_mean_weir_estimates # Mean FST
│   ├── generate_fst_matrix       # Create matrix
│   └── plot_fst_matrix           # Visualize
│
├── hdv_pgx_pop                   # HDV detection
│   ├── get_freq_3_3              # Calculate frequencies
│   ├── split_id_freq_2           # Split CHROM:POS:REF:ALT
│   ├── hdv_dataset               # Identify HDV
│   └── combine_hdvs_pops         # Combine results
│
└── csq_pgx                       # Consequence analysis
    ├── csq                       # Predict consequences
    └── plot_csq                  # Visualize results
```

## Module Organization

### File Structure

```
module/
├── vcf_qc.nf              # Quality control processes
├── annotate_vcf.nf        # Annotation processes
├── fst.nf                 # FST calculation processes
├── freq_vcf.nf            # Frequency calculation
├── hdv_vcf.nf             # HDV detection
├── novel_vcf.nf           # Novel variant discovery
├── structure.nf           # PCA and population structure
├── subset_vcf.nf          # VCF subsetting utilities
├── exome_analysis.nf      # General analysis utilities
└── csq.nf                 # Consequence prediction
```

### Module Imports

From [main.nf:4-13](main.nf#L4-L13):

```groovy
include { dataset_qc_alt; get_vcf_site; ... } from './module/vcf_qc'
include { fst_analysis; combine_fst_weir_estimates; ... } from './module/fst'
include { check_files; singleton_pop; ... } from './module/exome_analysis'
include { annotate_snpeff; annotate_dbsnp; ... } from './module/annotate_vcf.nf'
include { get_sites as get_sites1; ... } from './module/subset_vcf.nf'
include { novel_in_existing_database; ... } from './module/novel_vcf'
include { get_freq; get_freq2; ... } from './module/freq_vcf'
include { hdv; hdv_dataset; ... } from './module/hdv_vcf'
include { csq; plot_csq } from './module/csq'
include { pruning_vcf; vcf_to_plink1; ... } from './module/structure'
```

## Data Flow

### 1. Input Preparation

```
Input VCF files (per dataset, per chromosome or wildcard)
    ↓
datasets_map_data Channel
    ↓ [dataset, vcf, vcf_idx, sample_file]
get_map (extract chromosome list from VCF)
    ↓ [dataset, vcf, vcf_idx, sample_file, map_file]
generate_chunks_vcf (create genomic chunks)
    ↓ [dataset, vcf, vcf_idx, sample_file, chunks_file]
flatMap (expand to individual chunks)
    ↓ [dataset, chrm, start, end, vcf, vcf_idx, sample_file]
split_vcf_to_chunk (extract chunk regions)
    ↓ [dataset, chunk_vcf, sample_file, chrm]
```

### 2. Annotation Pipeline

```
Chunk VCFs
    ↓
dataset_qc_alt (quality filtering)
    ↓
annotate_snpeff (functional annotation)
    ↓
annotate_dbsnp (dbSNP IDs and frequencies)
    ↓
annotate_clinvar (clinical significance)
    ↓
annotate_cosmic (cancer mutations)
    ↓
annotate_vcf_from_vcf (population MAFs)
    ↓
fill_tags (calculate AN, AC, AF tags)
    ↓
groupTuple by [dataset, chrm] (group chunks by chr)
    ↓
concat_chunks_vcf (merge chunks into chromosomes)
    ↓
Annotated chromosome VCFs
```

### 3. Population Processing

```
Annotated VCFs [dataset, vcf, sample, chrm]
    ↓
flatMap with dataset_files (add population info)
    ↓ [pop, dataset, vcf, sample, chrm]
split_dataset_vcf_pop (subset by population)
    ↓ [pop, dataset, pop_vcf, pop_sample, chrm]
combine with gene_lists
    ↓ [pop, pop_vcf, pop_sample, chrm, pgx_name, pgx_bed]
get_pgx_pop (extract PGx variants)
    ↓ [pop, pop_pgx_vcf, pop_sample, chrm]
groupTuple by [pop, dataset] (group chromosomes)
    ↓
concat_pop (concatenate chromosomes)
    ↓
Per-population PGx VCFs
```

### 4. Statistical Analysis Flow

```
Population VCFs
    ↓
    ├─→ Fisher's Exact Test
    │   ├─ vcf_to_plink (convert format)
    │   ├─ generate_pheno_fam_sample (create phenotypes)
    │   └─ fisher_test_plink (compute p-values)
    │
    ├─→ FST Analysis
    │   ├─ fst_analysis (calculate Weir & Cockerham FST)
    │   ├─ get_fst_weir_estimates (extract per-variant FST)
    │   ├─ combine_fst_weir_estimates (merge results)
    │   └─ generate_fst_matrix (create population matrix)
    │
    └─→ HDV Detection
        ├─ get_freq_3_3 (calculate frequencies)
        ├─ hdv_dataset (identify highly differentiated variants)
        └─ combine_hdvs_pops (merge across populations)
```

## Process Details

### Key Process Categories

#### 1. VCF Manipulation

**split_vcf_to_chunk**
- Input: Full VCF, genomic coordinates
- Output: VCF subset for specific region
- Tool: bcftools view
- Purpose: Enable parallel processing

**concat_chunks_vcf**
- Input: Multiple VCF chunks for same chromosome
- Output: Single concatenated VCF
- Tool: bcftools concat
- Purpose: Merge parallel results

**split_dataset_vcf_pop**
- Input: Multi-population VCF
- Output: Population-specific VCF
- Tool: bcftools view -S
- Purpose: Population-level analysis

#### 2. Annotation

**annotate_snpeff**
- Input: VCF
- Output: VCF with functional annotations (ANN field)
- Tool: SnpEff
- Database: GRCh37.75
- Annotations: Gene, transcript, consequence, amino acid change

**annotate_dbsnp**
- Input: VCF
- Output: VCF with dbSNP IDs and frequencies
- Tool: bcftools annotate
- Purpose: Add rsIDs and population frequencies

**annotate_clinvar**
- Input: VCF
- Output: VCF with clinical significance
- Tool: bcftools annotate
- Purpose: Flag clinically relevant variants

**annotate_cosmic**
- Input: VCF
- Output: VCF with COSMIC IDs
- Tool: bcftools annotate
- Purpose: Flag cancer-related mutations

**annotate_vcf_from_vcf**
- Input: VCF, population MAF VCFs
- Output: VCF with population-specific frequencies
- Tool: bcftools annotate
- Annotations: KG_AF, gnomAD_AF, ExAC_AF, etc.

#### 3. Frequency Calculations

**get_freq**
- Input: VCF
- Output: Frequency file (.frq)
- Tool: vcftools --freq
- Purpose: Calculate allele frequencies

**get_freq_3_3**
- Input: VCF, population groups
- Output: Per-population frequencies
- Tool: bcftools +fill-tags, vcftools
- Purpose: Population-specific AF calculations

#### 4. Statistical Tests

**fisher_test_plink**
- Input: PLINK binary files (.bed/.bim/.fam), phenotype file
- Output: Fisher test results (.assoc.fisher)
- Tool: PLINK --fisher
- Purpose: Test allele frequency differences between groups

**fst_analysis**
- Input: VCF with two populations
- Output: Weir & Cockerham FST estimates
- Tool: vcftools --weir-fst-pop
- Purpose: Measure population differentiation

#### 5. Variant Filtering

**get_gene_vcf**
- Input: VCF, BED file
- Output: Variants overlapping BED regions
- Tool: bcftools view -R
- Purpose: Extract variants in specific genes/regions

**dataset_qc_alt**
- Input: VCF
- Output: Filtered VCF
- Tool: bcftools view
- Filters: PASS, biallelic SNPs, remove low quality

#### 6. Format Conversion

**vcf_to_plink**
- Input: VCF
- Output: PLINK binary files (.bed/.bim/.fam)
- Tool: PLINK --make-bed
- Purpose: Enable PLINK-based analyses

**vcf_to_plink1**
- Input: VCF
- Output: PLINK files for PCA
- Tool: PLINK --make-bed with filters
- Purpose: Prepare for SMARTPCA

#### 7. Population Structure

**smartpca_dataset**
- Input: PLINK eigenstrat files
- Output: PCA results (.evec, .eval)
- Tool: smartpca (EIGENSOFT)
- Purpose: Compute principal components

**plot_pca_group**
- Input: PCA results
- Output: PCA plots (PDF)
- Tool: R (ggplot2)
- Purpose: Visualize population structure

## Channel Operations

### Common Patterns

#### 1. flatMap for Data Expansion

```groovy
// Expand dataset to include all populations
dataset_dataset_pops = data.flatMap{ dataset, dataset_vcf, dataset_sample, chrm ->
    add_pop_groups_datas = []
    params.dataset_files.each{ dataset_name, dataset_vcf_, dataset_sample_, dataset_pops ->
        dataset_pops.split(',').each{ pop ->
            if(dataset == dataset_name){
                add_pop_groups_datas << [pop, dataset, file(dataset_vcf), file(dataset_sample), chrm]
            }
        }
    }
    return add_pop_groups_datas
}
```

#### 2. groupTuple for Aggregation

```groovy
// Group chunks by dataset and chromosome
concat_chunks_vcf(fill_tags.out.groupTuple(by:[0,3]))
// Input:  [dataset, chunk1_vcf, sample, chrm]
//         [dataset, chunk2_vcf, sample, chrm]
//         [dataset, chunk3_vcf, sample, chrm]
// Output: [dataset, [chunk1_vcf, chunk2_vcf, chunk3_vcf], sample, chrm]
```

#### 3. combine for Cartesian Product

```groovy
// Combine each VCF with gene lists
pop_data = split_dataset_vcf_pop.out[0]
    .combine(gene_lists)
    .map{ pop, dataset, pop_vcf, pop_sample, chrm, pgx_name, pgx_bed ->
        [ pop, file(pop_vcf), file(pop_sample), chrm, pgx_name, file(pgx_bed) ]
    }
```

#### 4. map for Transformation

```groovy
// Reorder and transform channel elements
get_freq_3_2(concat_vcf.out.map{ dataset, vcf, sample ->
    [ dataset, file(vcf), file(sample), '' ]
})
```

#### 5. filter for Selection

```groovy
// Select only AIBST dataset
hdv_all_chrm_dataset_aibst = data_proc.out.fill_tags_dataset_chrm
    .flatMap{ dataset, vcf, sample, chrm ->
        datas = []
        if(dataset == 'AIBST'){
            datas << [ dataset, file(vcf), file(sample), chrm ]
        }
        return datas
    }
```

### Channel Element Structure

Common channel tuple structures used throughout:

```groovy
// Dataset VCF
[dataset_name, vcf_file, sample_file, chromosome]

// Population VCF
[population_code, dataset_name, vcf_file, sample_file, chromosome]

// PGx VCF
[dataset_name, vcf_file, sample_file, chromosome, gene_list_name, bed_file]

// Frequency data
[dataset_name, vcf_file, sample_file, sites_file, maf_file, frq_file, chromosome]

// FST data
[population_pair, weir_fst_file, mean_fst_file]

// HDV data
[population_code, hdv_all_file, hdv_base_file, fold_threshold, ac_threshold, group_threshold]
```

## Error Handling

### Retry Strategy

Configured in process scope:

```groovy
process {
    errorStrategy = 'retry'
    maxRetries = 3
    maxErrors = 10000

    // Memory scaling on retry
    memory = { 30.GB * task.attempt }

    // CPU scaling on retry
    cpus = { 1 * task.attempt }
}
```

### Process-Specific Error Handling

```groovy
withLabel: 'smartpca' {
    // Retry on specific exit codes
    errorStrategy = {
        if (task.exitStatus in [143, 137, 255]) {
            'retry'
        } else {
            'terminate'
        }
    }
}
```

### Exit Codes

- **143**: SIGTERM (timeout or resource limit)
- **137**: SIGKILL (out of memory)
- **255**: General error

### Input Validation

```groovy
// Check files exist before processing
check_files([
    params.dbsnp_vcf,
    "${params.dbsnp_vcf}.tbi",
    params.clinvar,
    "${params.clinvar}.tbi"
])
```

## Performance Optimization

### 1. Chunked Processing

**Benefit:** Reduces memory requirements, enables parallelization

```groovy
chunk_size = 25000000  // 25MB genomic chunks
```

- Splits large chromosomes into manageable pieces
- Each chunk processed independently
- Results concatenated at chromosome level

### 2. Parallel Execution

**Chromosome-level parallelism:**
- Each chromosome processed independently
- 22 autosomes = up to 22 parallel jobs

**Chunk-level parallelism:**
- Each chunk within chromosome processed in parallel
- 100+ chunks for full genome = 100+ parallel jobs

**Population-level parallelism:**
- Each population processed independently
- 12 AIBST populations = 12 parallel jobs

### 3. Resource Allocation

**Label-based resource assignment:**

```groovy
withLabel: "small" { memory = 8.GB }      // Simple bcftools operations
withLabel: "medium" { memory = 18.GB }     // VCFtools, PLINK
withLabel: "bigmem" { memory = 28.GB }     // SnpEff annotation
withLabel: "extrabig" { memory = 200.GB }  // Full dataset merging
```

### 4. Caching with -resume

Nextflow caches successful process outputs:

```bash
# Resume from last successful step
nextflow run main.nf -resume
```

**Cache location:** `work/` directory

**When cache is valid:**
- Input files unchanged
- Process script unchanged
- Process configuration unchanged

## Workflow Composition

### Nested Workflows

Workflows can call other workflows:

```groovy
workflow {
    main:
        data_proc()  // Call data_proc workflow
        data_proc_pop(data_proc.out.fill_tags_dataset_chrm)  // Pass output to next workflow
}
```

### Workflow Outputs

Workflows can emit multiple named outputs:

```groovy
workflow data_proc {
    main:
        // Processing steps...
    emit:
        annotate_cosmic
        qc_dataset_chrm = dataset_qc_alt.out
        fill_tags_dataset_chrm = concat_chunks_vcf.out
        pgx_dataset = concat_vcf.out
        datasets_map_data
}
```

### Output Usage

```groovy
// Access workflow outputs
data_proc_pop(data_proc.out.fill_tags_dataset_chrm)
group_pops_pgx(data_proc_pop.out.pop_pgx)
```

## Best Practices

1. **Use explicit channel operations** - Avoid implicit channel creation
2. **Label processes** - Use descriptive labels for resource allocation
3. **Group related processes** - Organize into domain-specific modules
4. **Document channel structure** - Comment tuple element meanings
5. **Validate inputs early** - Use check_files at workflow start
6. **Emit meaningful outputs** - Name workflow outputs clearly
7. **Use consistent naming** - Follow naming conventions across modules
8. **Test incrementally** - Test new processes in isolation first

## Extending the Workflow

### Adding a New Process

1. Define process in appropriate module file
2. Import in main.nf
3. Add to workflow
4. Specify resource requirements
5. Test with small dataset

### Adding a New Analysis

1. Create new workflow block
2. Define inputs from existing workflows
3. Chain processes
4. Emit relevant outputs
5. Integrate into main workflow

### Example: Adding New Statistical Test

```groovy
// In module/stats.nf
process my_new_test {
    label 'medium'

    input:
    tuple val(dataset), path(vcf), path(sample)

    output:
    tuple val(dataset), path("${dataset}.result.txt")

    script:
    """
    my_statistical_tool --vcf ${vcf} --out ${dataset}.result.txt
    """
}

// In main.nf
include { my_new_test } from './module/stats.nf'

workflow new_stats_analysis {
    take: data
    main:
        my_new_test(data)
    emit:
        results = my_new_test.out
}

// Add to main workflow
workflow {
    main:
        data_proc()
        new_stats_analysis(data_proc.out.fill_tags_dataset_chrm)
}
```

## References

- [Nextflow Documentation](https://www.nextflow.io/docs/latest/)
- [Nextflow DSL2 Guide](https://www.nextflow.io/docs/latest/dsl2.html)
- [Nextflow Patterns](https://nextflow-io.github.io/patterns/)
