#!/usr/bin/env nextflow
/*
========================================================================================
    African Exome Analysis Workflow - 2025 Edition
========================================================================================
    Author: Mamana Mbiyavanga
    Affiliation: University of Cape Town
    Description: Modernized Nextflow DSL2 pipeline for comprehensive exome analysis
                 with focus on African populations and pharmacogenomics

    Version: 2.0.0
    Date: 2025-01-11

    Improvements over main.nf:
    - Modular workflow organization with clear separation of concerns
    - Simplified channel operations for better maintainability
    - Enhanced error handling and validation
    - Removed ~300 lines of commented/dead code
    - Consistent naming conventions and documentation
    - Optimized performance through better channel management

    Usage:
        nextflow run main.nf -profile singularity,slurm -c CONFIG_FILE

    For detailed documentation, see:
        - README.md: Quick start and overview
        - WORKFLOW_ARCHITECTURE.md: Technical details
        - CONFIGURATION.md: Configuration guide
========================================================================================
*/

nextflow.enable.dsl=2

/*
========================================================================================
    MODULE IMPORTS
========================================================================================
*/

// ===== Quality Control & Data Processing =====
include {
    dataset_qc_alt
    get_vcf_site
    tabix_data_proc_pop
    sites_only
} from './module/vcf_qc'

// ===== FST Analysis & Population Differentiation =====
include {
    fst_analysis
    combine_fst_weir_estimates
    combine_fst_weir_estimates as combine_fst_weir_estimates_cutoff
    combine_fst_mean_weir_estimates
    combine_fst_mean_weir_estimates as combine_fst_mean_weir_estimates_1
    combine_fst_mean_weir_estimates as combine_fst_mean_weir_estimates_2
    get_fst_weir_estimates
    get_fst_weir_estimates_cutoff
    combine_weir_fst_analysis
    plot_fst_matrix
    generate_fst_matrix
} from './module/fst'

// ===== General Utilities & File Operations =====
include {
    check_files
    singleton_pop
    split_dataset_vcf_pop
    split_dataset_vcf_pop_tabix
    concat_pop
    group_pops_data
    merge_pop_groups
    cat_chrm_groups
    combine_csv
} from './module/exome_analysis'

// ===== Annotation Pipeline =====
include {
    annotate_snpeff
    annotate_dbsnp
    annotate_clinvar
    annotate_cosmic
    combine_mafs_dataset
    fix_mafs_headers
    annotate_mafs
    filter_empty_vcf
    clean_invalid_variants
    fill_tags
    get_pop
    concat_mafs
    concat_vcf
    annotate_vcf_from_vcf
    annotate_gnomad_afr
    concat_chunks_vcf
    reheader_samples
} from './module/annotate_vcf.nf'

// ===== chrX Ploidy Correction =====
include {
    infer_sex_chrx
    fix_chrx_ploidy
} from './module/chrx_ploidy'


// ===== Pathogenicity Score Annotations =====
include {
    annotate_cadd
    annotate_alphamissense
    annotate_dbnsfp
} from './module/pathogenicity_scores'

// ===== PGx Star Allele Calling =====
include { pgx_star_alleles } from './module/pgx_calling'

// ===== HLA Typing =====
include { hla_typing } from './module/hla_typing'
// ===== VCF Subsetting & Region Extraction =====
include {
    get_sites as get_sites1
    get_sites as get_sites2
    get_gene_vcf
    get_gene_vcf1 as get_novel_vcf
    get_gene_vcf as get_pgx_pop
    get_gene_vcf_simple
    get_gene_vcf1
    vcf_to_plink
    generate_pheno_fam_sample
    generate_pheno_fam_sample1
    get_protein_coding_variants
    get_vcf_ind
    merge_groups
    generate_chunks_vcf
    split_vcf_to_chunk
    filter_alt_contigs
} from './module/subset_vcf.nf'

// ===== Novel Variant Discovery =====
include {
    novel_in_existing_database
    get_novel_sites
    novel_count
    novel_count_combine
    combine_novel_count as combine_novel_count3
    combine_novel_count as combine_novel_count5
    combine_novel_count as combine_novel_count4
    combine_novel_count as combine_novel_count1
    combine_novel_count as combine_novel_count2
    count_variants
    combine_counts
    generate_circos_dataframe
    generate_circos_plot
    singleton_dataset_chrm
} from './module/novel_vcf'

// ===== Frequency Calculations =====
include {
    get_freq
    get_freq2
    get_freq_3
    get_freq_3_1
    cat_freq
    cat_sites
    data_for_upset
    plot_upset_R
    fisher_test_plink
    fisher_test_plink_pops
    filter_fisher_test
    split_id_freq
    split_id_freq_1
    split_id_freq_2
    get_sites_simple
    get_map_vcf
    get_map
    get_freq_3_2
    get_freq_3_3
    get_freq_from_vcf_sites
    add_clinvar_to_freq
    annotate_pop_mafs
} from './module/freq_vcf'

// ===== Highly Differentiated Variants (HDV) =====
include { hdv; hdv_dataset; combine_hdvs; annotate_hdvs } from './module/hdv_vcf'
include { combine_hdvs as combine_hdvs_all } from './module/hdv_vcf'
include { combine_hdvs as combine_hdvs_base } from './module/hdv_vcf'
include { combine_hdvs as combine_hdvs_base_all } from './module/hdv_vcf'

// ===== HDV FST Analysis and Plots =====
include {
    compute_pgx_fst
    plot_fst_distribution
    generate_hdv_supp_table
} from './module/hdv_fst_plots'

// ===== Consequence Analysis =====
include {
    csq
    plot_csq
} from './module/csq'

// ===== Population Structure & PCA =====
include {
    pruning_vcf
    vcf_to_plink1
    smartpca_dataset
    plot_pca_group
    update_evec
    updated_ped
} from './module/structure'
// ===== Advanced Population Structure & Selection Scans =====
include {
    admixture
    selection_scans
} from './module/structure_advanced'

/*
========================================================================================
    HELPER FUNCTIONS
========================================================================================
*/

/**
 * Generate population to dataset mapping
 * Creates a list of [population, dataset_name] tuples for all populations
 *
 * @return List<Tuple> - List of [pop, dataset] pairs
 */
def get_dataset_pop() {
    def add_pop_groups_datas = []
    params.dataset_files.each { dataset_name, dataset_vcf_, dataset_sample_, dataset_pops ->
        // Add dataset name itself
        add_pop_groups_datas << [dataset_name, dataset_name]
        // Add all populations in this dataset
        dataset_pops.split(',').each { pop ->
            add_pop_groups_datas << [pop, dataset_name]
        }
    }
    return add_pop_groups_datas
}

/**
 * Filter data for specific datasets
 * Helper to extract only specified dataset entries
 *
 * @param data Channel - Input data channel
 * @param dataset_name String - Dataset to filter for
 * @return Filtered channel entries
 */
def filter_by_dataset(data, dataset_name) {
    return data.filter { it[0] == dataset_name }
}


/*
========================================================================================
    MAIN DATA PROCESSING WORKFLOW
========================================================================================
    Core workflow for VCF processing, annotation, and initial filtering

    Inputs: Raw VCF files from params.dataset_files
    Outputs:
        - annotated_vcfs: Fully annotated chromosome VCFs
        - qc_dataset_chrm: Quality-controlled VCFs before annotation
        - fill_tags_dataset_chrm: VCFs with filled INFO tags
        - singl_dataset_chrm: Singleton variant sites per chromosome
        - pgx_dataset_chrm: PGx variant VCFs per chromosome
        - pgx_dataset: Concatenated PGx variants (all chromosomes)
        - all_datasets_annotated: Merged datasets with annotations
        - datasets_map_data: Original dataset mapping information
========================================================================================
*/

workflow data_proc {
    main:
        // ===== Step 1: Validate Input Files =====
        log.info "="*80
        log.info "Starting Data Processing Workflow"
        log.info "="*80

        // Check critical reference files exist
        def files_to_check = [
            params.dbsnp_vcf,
            "${params.dbsnp_vcf}.tbi",
            params.clinvar,
            "${params.clinvar}.tbi"
        ]

        // Only check COSMIC if enabled
        if (params.use_cosmic) {
            files_to_check.add(params.cosmic)
            files_to_check.add("${params.cosmic}.tbi")
        }

        check_files(files_to_check)

        // ===== Step 2: Prepare Dataset Channel =====
        // Build channel with dataset VCF files and their indices
        // Filter by params.chromosomes if specified
        datasets_map = []
        def chroms_to_process = params.chromosomes as List
        log.info "Chromosomes to process: ${chroms_to_process.join(',')}"

        params.dataset_files.each { dataset_name, dataset_vcf, dataset_sample, dataset_pops ->
            def vcfs = file(dataset_vcf)

            // Handle both single files and file lists
            if (vcfs instanceof List) {
                vcfs.each { vcf_ ->
                    // Extract chromosome from filename and filter
                    def vcf_name = vcf_.getName()
                    def should_process = false

                    // Check if this file matches any requested chromosome
                    chroms_to_process.each { chr ->
                        if (vcf_name.contains("chr${chr}_") ||
                            vcf_name.contains("chr${chr}.") ||
                            vcf_name.contains("chrm${chr}_") ||
                            vcf_name.contains("chrm${chr}.") ||
                            vcf_name.contains("_${chr}_") ||
                            vcf_name.contains("_${chr}.vcf")) {
                            should_process = true
                        }
                    }

                    if (should_process) {
                        def vcf_idx = determineIndexFile(vcf_)
                        check_files([vcf_, vcf_idx, dataset_sample])
                        datasets_map << [dataset_name, file(vcf_), file(vcf_idx), file(dataset_sample)]
                        log.info "Including: ${vcf_name}"
                    } else {
                        log.info "Skipping: ${vcf_name} (not in requested chromosomes)"
                    }
                }
            } else {
                def vcf_idx = determineIndexFile(dataset_vcf)
                check_files([dataset_vcf, vcf_idx, dataset_sample])
                datasets_map << [dataset_name, file(dataset_vcf), file(vcf_idx), file(dataset_sample)]
            }
        }

        datasets_map_data = Channel.from(datasets_map)

        // ===== Optional: Rename VCF samples =====
        // If sample_rename_file is set, rename samples to match the sample mapping file
        if (params.sample_rename_file && file(params.sample_rename_file).exists()) {
            log.info "Renaming VCF samples using: ${params.sample_rename_file}"
            reheader_samples(datasets_map_data)
            datasets_ch = reheader_samples.out
        } else {
            datasets_ch = datasets_map_data
        }

        // ===== chrX Sex Inference =====
        // Infer sample sex from chrX heterozygosity rates (runs once per dataset)
        chrx_vcf_for_sex = datasets_ch
            .filter { dataset, vcf, idx, sample -> vcf.name =~ /chrX/ }
            .map { dataset, vcf, idx, sample -> [dataset, vcf, idx] }
        infer_sex_chrx(chrx_vcf_for_sex)

        // ===== Step 3: Generate Genomic Chunks =====
        // Extract chromosome information from VCF
        get_map(datasets_ch)

        // Filter alternate contigs before chunking
        filter_alt_contigs(get_map.out)

        // Generate chunk coordinates for parallel processing (main contigs only)
        generate_chunks_vcf(
            filter_alt_contigs.out.map { dataset, vcf, vcf_idx, sample_file, map_file ->
                [dataset, file(vcf), file(vcf_idx), file(sample_file), file(map_file), '', params.chunk_size]
            }
        )

        // Expand chunks into individual processing units
        chunks_datas = generate_chunks_vcf.out.flatMap { dataset, vcf, vcf_idx, sample_file, chunk_file ->
            def datas = []
            log.info "DEBUG flatMap: chunk_file=${chunk_file}, class=${chunk_file.getClass()}"
            // CephFS workaround: retry file read if empty (metadata cache lag)
            def content = ''
            for (int attempt = 0; attempt < 5; attempt++) {
                content = ['cat', chunk_file.toString()].execute().text
                if (content.trim().length() > 0) break
                log.warn "flatMap: chunk_file empty on attempt ${attempt+1}, retrying in 5s..."
                Thread.sleep(5000)
            }
            log.info "DEBUG flatMap: content length=${content.length()}, first 100 chars=${content.take(100)}"
            def lines = content.trim().split('\n')
            log.info "DEBUG flatMap: ${lines.size()} lines found"
            lines.each { line ->
                def data = line.trim().split(',')
                if (data.size() >= 3) {
                    datas << [dataset, data[0], data[1], data[2], vcf, vcf_idx, sample_file]
                }
            }
            log.info "DEBUG flatMap: returning ${datas.size()} chunks"
            return datas
        }

        // Extract chunk regions from VCF
        split_vcf_to_chunk(chunks_datas)

        // ===== Step 4: Quality Control =====
        // Filter variants: PASS only, biallelic SNPs, quality thresholds
        dataset_qc_alt(split_vcf_to_chunk.out)

        // ===== Step 5: Annotation Pipeline =====
        log.info "Starting Annotation Pipeline"

        // 5.1: Functional annotation with SnpEff
        annotate_snpeff(
            dataset_qc_alt.out.combine([[params.snpeff_human_db, params.snpeff_database]])
        )

        // 5.2: Add dbSNP IDs and frequencies
        annotate_dbsnp(
            annotate_snpeff.out.combine([[file(params.dbsnp_vcf), file("${params.dbsnp_vcf}.tbi")]])
        )

        // 5.3: Add ClinVar clinical significance annotations
        annotate_clinvar(
            annotate_dbsnp.out.combine([[file(params.clinvar), file("${params.clinvar}.tbi")]])
        )

        // 5.4: Add COSMIC cancer mutation annotations (conditional)
        if (params.use_cosmic) {
            annotate_cosmic(
                annotate_clinvar.out.combine([[file(params.cosmic), file("${params.cosmic}.tbi")]])
            )
            annotation_output = annotate_cosmic.out
        } else {
            annotation_output = annotate_clinvar.out
        }

        // ===== Step 5.5: gnomAD v3 African Allele Frequency =====
        // Add AF_afr from gnomAD v3 (global AF already in GNOMAD_ALL_MAF from dbSNP)
        def gnomad_v3_pattern = params.mafs_annotations.find { it[0] == 'gnomAD_v4' }?[1]

        if (gnomad_v3_pattern) {
            gnomad_afr_input = annotation_output.map { dataset, vcf, sample, chrm, snpeff_db ->
                def chrm_num = chrm.replaceAll('chr', '')
                def gnomad_vcf = file(gnomad_v3_pattern.replace('%s', chrm_num))
                def gnomad_idx = file(gnomad_vcf.toString() + '.tbi')
                [dataset, vcf, sample, chrm, snpeff_db, gnomad_vcf, gnomad_idx]
            }
            annotate_gnomad_afr(gnomad_afr_input)
            pop_maf_output = annotate_gnomad_afr.out
        } else {
            pop_maf_output = annotation_output
        }

        // ===== Step 5.6: CADD Pathogenicity Scores =====
        // Add CADD v1.7 PHRED-scaled deleteriousness scores
        annotate_cadd(pop_maf_output)

        // ===== Step 5.7: AlphaMissense Predictions =====
        // Add protein structure-based pathogenicity predictions
        annotate_alphamissense(annotate_cadd.out)

        // ===== Step 5.8: dbNSFP Functional Predictions =====
        // Add SIFT, PolyPhen-2, REVEL, MetaRNN, GERP++ scores
        if (params.dbnsfp_db != null && params.dbnsfp_db && file(params.dbnsfp_db).exists() && file(params.dbnsfp_db + '.tbi').exists()) {
            annotate_dbnsfp(annotate_alphamissense.out)
            pathogenicity_output = annotate_dbnsfp.out
        } else {
            pathogenicity_output = annotate_alphamissense.out
        }

        // ===== Step 6: Filter Empty VCF Chunks =====
        // Skip empty chunks (header-only VCFs with no variants)
        filter_empty_vcf(pathogenicity_output)

        // ===== Step 6.5: Clean Invalid Variants =====
        // Remove variants with missing ALT alleles (ALT=".") which are liftover artifacts
        clean_invalid_variants(filter_empty_vcf.out)

        // ===== Step 6.75: chrX Ploidy Correction =====
        // Fix male genotypes on chrX non-PAR regions (diploid -> haploid)
        // Must happen before fill_tags so AF/AC/AN are calculated correctly
        clean_invalid_variants.out.branch {
            chrx: it[3] =~ /chrX|^X$/
            other: true
        }.set { cleaned_by_chrm }

        fix_chrx_ploidy(
            cleaned_by_chrm.chrx.combine(
                infer_sex_chrx.out.map { dataset, sex_map -> [dataset, sex_map] },
                by: 0
            )
        )

        ploidy_corrected = fix_chrx_ploidy.out.mix(cleaned_by_chrm.other)

        // ===== Step 7: Fill INFO Tags =====
        // Recalculate AC, AN, AF, etc.
        fill_tags(ploidy_corrected)

        // ===== Step 8: Concatenate Chunks =====
        // Merge chunks back into chromosome-level VCFs
        concat_chunks_vcf(fill_tags.out.groupTuple(by: [0, 3]))

        // ===== Step 8: Identify Singletons =====
        // Extract singleton variants per chromosome for each dataset
        singleton_dataset_chrm(
            concat_chunks_vcf.out.map { dataset, vcf, sample, chrm ->
                [dataset, file(vcf), chrm]
            }
        )

        // ===== Step 9: Extract PGx Variants =====
        // Extract pharmacogenomic variants using gene lists
        get_gene_vcf(
            concat_chunks_vcf.out.combine(
                Channel.from(params.gene_lists).map { name, bed ->
                    [name, file(bed)]
                }
            )
        )

        // Concatenate PGx variants across chromosomes
        concat_vcf(get_gene_vcf.out.groupTuple())

        // Generate sites-only version (no genotypes)
        sites_only(concat_vcf.out)

        // Annotate per-population MAFs into VCF
        annotate_pop_mafs(concat_vcf.out)

        // Calculate frequencies for PGx variants
        get_freq_3_2(
            annotate_pop_mafs.out.map { dataset, vcf, sample ->
                [dataset, file(vcf), file(sample), '']
            }
        )
        split_id_freq_1(get_freq_3_2.out)

        // ===== Step 10: Merge Datasets =====
        // Combine all datasets for cross-dataset analysis
        merge_datasets_group(sites_only.out)

        log.info "Data Processing Workflow Completed"

    emit:
        annotation_output  // Either annotate_cosmic.out or annotate_clinvar.out depending on use_cosmic flag
        qc_dataset_chrm = dataset_qc_alt.out
        fill_tags_dataset_chrm = concat_chunks_vcf.out
        singl_dataset_chrm = singleton_dataset_chrm.out
        pgx_dataset_chrm = get_gene_vcf.out
        pgx_dataset = concat_vcf.out
        pop_mafs_vcf = annotate_pop_mafs.out
        all_datasets_annotated = merge_datasets_group.out.all_dataset_annotation
        datasets_map_data
}

/**
 * Helper function to determine index file based on VCF type
 */
def determineIndexFile(vcf_file) {
    def vcf_path = vcf_file instanceof String ? vcf_file : vcf_file.toString()
    if (file(vcf_path).getExtension() == "gz") {
        return "${vcf_path}.tbi"
    } else if (file(vcf_path).getExtension() == "bcf") {
        return "${vcf_path}.csi"
    }
    return "${vcf_path}.tbi"  // default to .tbi
}

/*
========================================================================================
    POPULATION MAF ANNOTATION WORKFLOW
========================================================================================
    Process population allele frequency annotations from multiple databases

    Prepares MAF annotation files for cross-referencing with study variants
========================================================================================
*/

workflow process_mafs_annot {
    main:
        def datas = []
        params.mafs_annotations.each { dataset, vcf ->
            def vcfs = file(vcf)

            if (vcfs instanceof List) {
                vcfs.each { vcf_ ->
                    def vcf_idx = determineIndexFile(vcf_)
                    check_files([vcf_, vcf_idx])
                    datas << [dataset, file(vcf_), file(vcf_idx)]
                }
            } else {
                def vcf_idx = determineIndexFile(vcf)
                check_files([vcf, vcf_idx])
                datas << [dataset, file(vcf), file(vcf_idx)]
            }
        }

        mafs_ = Channel.from(datas)
        get_vcf_chrom(mafs_)

    emit:
        mafs = get_vcf_chrom.out.vcf_chrm
}

/*
========================================================================================
    VCF CHROMOSOME EXTRACTION WORKFLOW
========================================================================================
    Extract chromosome-specific data from multi-chromosome VCF files
========================================================================================
*/

workflow get_vcf_chrom {
    take:
        data

    main:
        get_map_vcf(data)

        vcf_chrm = get_map_vcf.out.flatMap { dataset, vcf, vcf_idx, dataset_map ->
            def chrm_annot_data = []
            def chrms = file(dataset_map).text.split()
            chrms.each { chrm ->
                chrm_annot_data << [dataset, vcf, vcf_idx, chrm]
            }
            return chrm_annot_data
        }

    emit:
        vcf_chrm
}

/*
========================================================================================
    DATASET MERGING WORKFLOW
========================================================================================
    Merge all datasets and add comprehensive annotations
========================================================================================
*/

workflow merge_datasets_group {
    take:
        data

    main:
        // Merge all datasets into single "ALL" dataset
        merge_groups(
            data.map { dataset, vcf, sample ->
                ['ALL', file(vcf)]
            }.groupTuple()
        )

        // Calculate frequencies for merged dataset
        get_freq_from_vcf_sites(
            merge_groups.out.map { dataset, vcf ->
                [dataset, file(vcf), '', '']
            }
        )

        // Add ClinVar, PharmGKB, and GWAS annotations
        check_files([params.clinvar1])
        add_clinvar_to_freq(
            get_freq_from_vcf_sites.out.map { dataset, vcf_file, sample_file, chrm, freq_file ->
                [dataset, freq_file, file(params.clinvar1), file(params.pharmgkb_file), file(params.gwas_file)]
            }
        )

    emit:
        all_dataset_annotation = add_clinvar_to_freq.out
}

/*
========================================================================================
    POPULATION-LEVEL DATA PROCESSING WORKFLOW
========================================================================================
    Split datasets by population and extract population-specific variants

    Takes chromosome-level annotated VCFs and:
    1. Splits by population
    2. Extracts PGx variants per population
    3. Concatenates chromosomes for each population
========================================================================================
*/

workflow data_proc_pop {
    take:
        data

    main:
        log.info "Starting Population-Level Processing"

        // ===== Step 1: Prepare Gene Lists =====
        gene_lists = Channel.from(params.gene_lists)
            .map { name, bed -> [name, file(bed)] }

        // ===== Step 2: Split Datasets by Population =====
        // Expand dataset data to include all populations
        dataset_dataset_pops = data.flatMap { dataset, dataset_vcf, dataset_sample, chrm ->
            def add_pop_groups_datas = []
            params.dataset_files.each { dataset_name, dataset_vcf_, dataset_sample_, dataset_pops ->
                dataset_pops.split(',').each { pop ->
                    if (dataset == dataset_name) {
                        add_pop_groups_datas << [pop, dataset, file(dataset_vcf), file(dataset_sample), chrm]
                    }
                }
            }
            return add_pop_groups_datas
        }

        // Split VCF by population
        split_dataset_vcf_pop(dataset_dataset_pops)

        // ===== Step 3: Extract PGx Variants per Population =====
        pop_data = split_dataset_vcf_pop.out[0]
            .combine(gene_lists)
            .map { pop, dataset, pop_vcf, pop_sample, chrm, pgx_name, pgx_bed ->
                [pop, file(pop_vcf), file(pop_sample), chrm, pgx_name, file(pgx_bed)]
            }

        get_pgx_pop(pop_data)

        // Add dataset information back
        pop_pgx_data_chrm = get_pgx_pop.out
            .combine(get_dataset_pop(), by: 0)
            .map { pop, pop_vcf, pop_sample, chrm, dataset ->
                [pop, dataset, pop_vcf, pop_sample, chrm]
            }

        // ===== Step 4: Concatenate Chromosomes per Population =====
        concat_pop(pop_pgx_data_chrm.groupTuple())

        pop_pgx_data = concat_pop.out.map { pop, dataset, pop_vcf, pop_sample ->
            [pop, file(pop_vcf), file(pop_sample), dataset]
        }

        log.info "Population-Level Processing Completed"

    emit:
        pop_exome_chrm = split_dataset_vcf_pop.out
        pop_exome = concat_pop.out
        pop_pgx = pop_pgx_data
        pop_pgx_chrm = pop_pgx_data_chrm
}

/*
========================================================================================
    POPULATION STRUCTURE ANALYSIS - PCA
========================================================================================
    Principal Component Analysis for population structure
========================================================================================
*/

workflow pca {
    take:
        data

    main:
        // LD pruning
        pruning_vcf(
            data.map { group_3, chrm_3, group_vcf_3, group_sample_3 ->
                [group_3, file(group_vcf_3), file(group_sample_3)]
            }
        )

        // Group pruned chromosomes
        groups_data = pruning_vcf.out.groupTuple(by: [0])
        cat_chrm_groups(groups_data)

        // Convert to PLINK format
        vcf_to_plink1(cat_chrm_groups.out)

        // Update PED files with population information
        updated_ped(vcf_to_plink1.out)

        // Run SMARTPCA
        smartpca_dataset(updated_ped.out)

        // Update evec file with population labels
        update_evec(smartpca_dataset.out)

        // Generate PCA plots
        plot_pca_group(update_evec.out)

    emit:
        data
        plink_data = vcf_to_plink1.out
}

/*
========================================================================================
    PCA FOR PGX VARIANTS
========================================================================================
    PCA analysis specifically for pharmacogenomic variants
========================================================================================
*/

workflow pca_pgx {
    take:
        data

    main:
        pruning_vcf(data)
        vcf_to_plink1(pruning_vcf.out)
        updated_ped(vcf_to_plink1.out)
        smartpca_dataset(updated_ped.out)
        update_evec(smartpca_dataset.out)
        plot_pca_group(update_evec.out)

    emit:
        data
        plink_data = vcf_to_plink1.out
}

/*
========================================================================================
    To be continued in next section...
========================================================================================
*/

/*
========================================================================================
    HDV (HIGHLY DIFFERENTIATED VARIANTS) ANALYSIS WORKFLOWS
========================================================================================
*/

/**
 * HDV Analysis for All Chromosomes (Dataset-level)
 * Identifies variants with significantly different frequencies across populations
 */
workflow hdv_all_chrm_dataset {
    take:
        data

    main:
        // Calculate allele frequencies
        get_freq(data)
        split_id_freq(get_freq.out)

        // Concatenate frequencies across chromosomes
        cat_freq(
            split_id_freq.out.map { dataset, vcf, sample, frq1, frq2, frq3, chrm ->
                [dataset, file(frq2), chrm]
            }.groupTuple()
        )

    emit:
        data
}

/**
 * HDV Analysis for PGx Variants (Dataset-level)
 * Identifies highly differentiated pharmacogenomic variants
 */
workflow hdv_pgx_dataset {
    take:
        data

    main:
        // Calculate frequencies with p-value file
        get_freq_3_1(
            data.map { dataset, vcf, sample, chrm, group, pvalue_file ->
                [dataset, file(vcf), file(sample), chrm, file(pvalue_file)]
            }
        )
        split_id_freq(get_freq_3_1.out)

        // HDV detection parameters
        def fold = 2           // Fold change threshold
        def ac = 3             // Allele count threshold
        def group_treshold = 1 // Minimum number of populations
        def test_pops = 'AGVP_AF,KG_AF,KG_AFR_AF,gnomAD_AF,gnomAD_AFR_AF,TOPMED_AF'

        // Identify HDV variants
        hdv_dataset(
            split_id_freq.out.map { dataset, freq_file, pvalue_file, chrm ->
                [dataset, file(freq_file), file(pvalue_file), chrm, fold, test_pops, ac, group_treshold, params.fst_cutoff]
            }
        )

        // Combine HDV results across chromosomes
        combine_hdvs(
            hdv_dataset.out.groupTuple(by: [0]).map { dataset, hdvs ->
                [dataset, hdvs, "dataset_${fold}fold_ac${ac}_g${group_treshold}"]
            }
        )

    emit:
        get_freq_3_1
        hdvs = hdv_dataset.out
}

/**
 * Combine HDV Results Across Populations
 * Aggregates highly differentiated variants from multiple populations
 */
workflow combine_hdvs_pops {
    take:
        data

    main:
        // Combine ALL HDV results
        combine_hdvs(
            data.map { pop, hdv_all, hdv_base, fold, ac, group_treshold ->
                if (file(hdv_all).countLines() > 1 && ac == 3 && fold == 2) {
                    return [pop, hdv_all, fold, ac, group_treshold]
                }
            }
            .groupTuple(by: [0, 2, 3, 4])
            .map { pop, hdvs_all, fold, ac, group_treshold ->
                [pop, hdvs_all, fold, ac, group_treshold, "pop_${fold}fold_ac${ac}_g${group_treshold}_all"]
            }
        )

        // Aggregate across all populations
        combine_hdvs_all(
            combine_hdvs.out
                .map { pop, hdv, fold, ac, group_treshold ->
                    ['ALL', pop, hdv, fold, ac, group_treshold]
                }
                .groupTuple(by: [0, 3, 4, 5])
                .map { dataset, pops, hdvs, fold, ac, group_treshold ->
                    [dataset, hdvs, fold, ac, group_treshold, "All_hdvs_pops_${fold}fold_ac${ac}_g${group_treshold}_all"]
                }
        )

        // Combine BASE HDV results
        combine_hdvs_base(
            data.map { pop, hdv_all, hdv_base, fold, ac, group_treshold ->
                if (file(hdv_base).countLines() > 1 && ac == 3 && fold == 2) {
                    return [pop, hdv_base, fold, ac, group_treshold]
                }
            }
            .groupTuple(by: [0, 2, 3, 4])
            .map { pop, hdvs_base, fold, ac, group_treshold ->
                [pop, hdvs_base, fold, ac, group_treshold, "pop_${fold}fold_ac${ac}_g${group_treshold}_base"]
            }
        )

        // Aggregate base HDVs across populations
        combine_hdvs_base_all(
            combine_hdvs_base.out
                .map { pop, hdv, fold, ac, group_treshold ->
                    ['ALL', pop, hdv, fold, ac, group_treshold]
                }
                .groupTuple(by: [0, 3, 4, 5])
                .map { dataset, pops, hdvs, fold, ac, group_treshold ->
                    [dataset, hdvs, fold, ac, group_treshold, "All_hdvs_pops_${fold}fold_ac${ac}_g${group_treshold}_base"]
                }
        )

        // Concatenate all HDV results
        hdvs_pops = combine_hdvs.out
            .concat(combine_hdvs_all.out, combine_hdvs_base.out, combine_hdvs_base_all.out)

    emit:
        hdvs_pops
}

/**
 * Annotate HDV Results with Functional Information
 */
workflow annotate_hdvs_all {
    take:
        data

    main:
        annotate_hdvs(
            data.map { dataset, dataset_hdv, fold, ac, group_treshold, group, group_annot ->
                [dataset, file(dataset_hdv), file(group_annot)]
            }
        )

    emit:
        annotate_hdvs.out
}

/**
 * HDV Analysis for PGx Variants (Population-level)
 * Most comprehensive HDV analysis combining FST and Fisher tests
 */
workflow hdv_pgx_pop {
    take:
        data
        pop_mafs_vcf  // annotated PGx VCF with per-population MAFs [dataset, vcf, sample]

    main:
        log.info "Starting Population-Level HDV Analysis for PGx Variants"

        // Combine per-pop data with the annotated PGx VCF (for cross-pop MAFs)
        freq_input = data.combine(
            pop_mafs_vcf.map { dataset, vcf, sample ->
                [file(vcf), file("${vcf}.tbi")]
            }
        )

        // Calculate frequencies
        get_freq_3_3(freq_input)
        split_id_freq_2(get_freq_3_3.out)

        // HDV parameters
        def fold_ch = Channel.of(2)
        def ac_ch = Channel.of(3)
        def group_treshold_ch = Channel.of(1, 2, 3, 4, 5)  // Multiple thresholds
        def fold_ac_group_treshold = fold_ch.combine(ac_ch).combine(group_treshold_ch)

        def base_pops = 'KNK,KNL,NGY,SAV,TZA,TZB,ZWD,ZWS'
        def test_pops = params.FST_TEST

        // Combine data with parameter combinations
        hdv_pop_data = split_id_freq_2.out
            .combine(fold_ac_group_treshold)
            .map { dataset, freq_file, chrm, pvalue_file, fst_file, fold_val, ac_val, group_treshold_val ->
                [dataset, file(freq_file), file(pvalue_file), file(fst_file), chrm, fold_val, test_pops, base_pops, ac_val, group_treshold_val, params.fst_cutoff]
            }

        // Identify HDV variants
        hdv_dataset(hdv_pop_data)

        // Combine results across populations
        combine_hdvs_pops(hdv_dataset.out)

        log.info "Population-Level HDV Analysis Completed"

    emit:
        hdvs_pops = combine_hdvs_pops.out.hdvs_pops
}

/*
========================================================================================
    VARIANT COUNTING & NOVEL DISCOVERY WORKFLOWS
========================================================================================
*/

/**
 * Count Variants for PGx Regions (Population-level)
 */
workflow count_pgx_pop {
    take:
        data

    main:
        count_data = data.map { dataset, chrm, pop, pop_vcf, pop_sample, dataset_vcf, dataset_singl_sites, dataset_singls ->
            [pop, file(pop_vcf), dataset, file(dataset_singl_sites), chrm]
        }
        count_dataset(count_data)

    emit:
        counts = count_dataset
}

/**
 * Count All Variants (Dataset-level)
 */
workflow count_all_dataset {
    take:
        data

    main:
        count_data = data.map { dataset, chrm, dataset_vcf, dataset_sample, dataset_vcf_, dataset_singl_sites, dataset_singls ->
            [dataset, file(dataset_vcf), dataset, file(dataset_singl_sites), chrm]
        }
        count_dataset(count_data)
        novel_dataset(count_data)

    emit:
        counts = count_dataset
}

/**
 * Count PGx Variants (Dataset-level)
 */
workflow count_pgx_dataset {
    take:
        data

    main:
        count_data = data.map { dataset, chrm, dataset_vcf, dataset_sample, dataset_vcf_, dataset_singl_sites, dataset_singls ->
            [dataset, file(dataset_vcf), dataset, file(dataset_singl_sites), chrm]
        }
        count_dataset(count_data)
        novel_dataset(count_data)

        // Concatenate novel variants
        novel_dataset_vcfs = novel_dataset.out.novel_dataset_vcf.groupTuple(by: 0)
        concat_vcf(novel_dataset_vcfs)

    emit:
        count_dataset
        novel_vcf = novel_dataset.out.novel_dataset_vcf
}

/**
 * Core Variant Counting Workflow
 * Counts variants per dataset/population/individual
 */
workflow count_dataset {
    take:
        data

    main:
        // Filter for AIBST dataset or individual samples
        chrm_dataset_data = data.flatMap { dataset, vcf, sample_id, singletons, chrm ->
            def datas = []
            // For datasets or individuals in AIBST
            if (dataset == 'AIBST' || sample_id == 'AIBST') {
                datas << [dataset, file(vcf), sample_id, file(singletons), chrm]
            }
            return datas
        }

        count_variants(chrm_dataset_data)

        // Aggregate counts across chromosomes
        dataset_counts = count_variants.out
            .flatMap { dataset, csv, sample_id, chrm ->
                [[dataset, dataset, file(csv), chrm]]
            }
            .groupTuple()
            .flatMap { dataset, sample_ids, count_csvs, chrms ->
                [[dataset, count_csvs, true, 'all']]
            }

        combine_counts(dataset_counts)

    emit:
        data
}

/**
 * Novel Variant Discovery Workflow
 * Identifies variants not present in reference databases
 */
workflow novel_dataset {
    take:
        data

    main:
        count_data = data.map { dataset, vcf, sample_id, singletons, chrm ->
            [dataset, file(vcf), sample_id, chrm]
        }

        // Identify variants not in other databases
        novel_in_existing_database_wf(count_data)
        combine_novel_count1_wf(novel_in_existing_database_wf.out.novel_in_dbs)

        // Extract novel variant sites
        novel_data = novel_in_existing_database_wf.out.novel_in_dbs
            .map { dataset, vcf, not_in_dbs, novel_sites, chrm ->
                if (!(file(novel_sites).isEmpty())) {
                    return [dataset, file(vcf), '', chrm, 'novel', file(novel_sites)]
                }
            }

        // Generate novel VCFs
        get_novel_vcf(novel_data)

        // Count novel variants
        count_data_novel = get_novel_vcf.out
            .map { dataset, dataset_novel_vcf, label, chrm ->
                [dataset, file(dataset_novel_vcf), dataset, label, chrm]
            }
            .combine(data, by: [0, 4])
            .map { dataset, chrm, dataset_novel_vcf, dataset_id1, label, dataset_vcf, dataset_id2, dataset_singletons ->
                [dataset, file(dataset_novel_vcf), dataset_id2, file(dataset_singletons), chrm]
            }

        count_variants(count_data_novel)

        // Aggregate novel counts
        dataset_novel_counts = count_variants.out
            .flatMap { dataset2, csv2, sample_id2, chrm2 ->
                def datas = []
                if (dataset2 == 'AIBST' || sample_id2 == 'AIBST') {
                    datas << [dataset2, sample_id2, file(csv2), chrm2]
                }
                return datas
            }
            .groupTuple(by: [0, 1])
            .flatMap { dataset3, sample_id, count_csvs3, chrms3 ->
                [[dataset3, count_csvs3, true, sample_id]]
            }

        combine_counts(dataset_novel_counts)

        // Combine totals
        combine_counts_data = combine_counts.out
            .groupTuple(by: [3])
            .map { datasets, count_totals, count_summaries, group ->
                [group, count_totals, true, group]
            }
        combine_combine_counts(combine_counts_data)

    emit:
        data
        novel_in_dbs = novel_in_existing_database_wf.out.novel_in_dbs
        novel_pop_vcf = count_variants.out
        novel_dataset_vcf = get_novel_vcf.out
}

/**
 * Helper: Identify Variants Not in Reference Databases
 */
workflow novel_in_existing_database_wf {
    take:
        data

    main:
        all_chrm_dataset_data = data.flatMap { dataset1, vcf1, sample_id1, chrm1 ->
            def datas = []
            if (dataset1 == 'AIBST' || sample_id1 == 'AIBST') {
                datas << [dataset1, file(vcf1), chrm1]
            }
            return datas
        }

        novel_in_existing_database(all_chrm_dataset_data)

    emit:
        novel_in_dbs = novel_in_existing_database.out
}

/**
 * Helper: Combine Novel Counts
 */
workflow combine_novel_count1_wf {
    take:
        data

    main:
        combine_data = data
            .flatMap { dataset1, vcf1, count1, csv1, chrm ->
                [[dataset1, file(count1)]]
            }
            .groupTuple()
            .flatMap { dataset2, counts2 ->
                [[dataset2, counts2, 'false', 'all']]
            }

        combine_novel_count1(combine_data)

        // Combine totals
        combine_counts_data = combine_novel_count1.out
            .groupTuple(by: [3])
            .map { datasets, count_totals, count_summaries, group ->
                [group, count_totals, true, group]
            }
        combine_combine_counts(combine_counts_data)

    emit:
        novel_counts = combine_novel_count1.out
}

/**
 * Helper: Aggregate Combined Counts
 */
workflow combine_combine_counts {
    take:
        data

    main:
        combine_counts(data)

    emit:
        data
}

/**
 * Novel Variants in All Chromosomes (Population-level)
 */
workflow novel_all_chrm_pop {
    take:
        data

    main:
        novel_in_existing_database(
            data.map { pop, dataset, vcf, sample ->
                [pop, file(vcf), dataset]
            }
        )

        get_freq(
            novel_in_existing_database.out.map { pop, vcf, count, csv, chrm ->
                [pop, file(vcf), file(count), chrm]
            }
        )

        get_novel_sites(
            get_freq.out.map { dataset, vcf, sample, sites, mafs, frq, chrm ->
                [dataset, file(vcf), file(mafs), '']
            }
        )

        get_gene_vcf(
            get_novel_sites.out.map { dataset, vcf, sites, chrm ->
                [dataset, file(vcf), '', chrm, '', file(sites)]
            }
        )

        novel_count(
            get_gene_vcf.out.map { dataset, vcf, sample, chrm ->
                [dataset, file(vcf), chrm]
            }
        )

        novel_count_combine(
            novel_count.out
                .map { dataset2, csv2, chrm2 ->
                    ['AIBST', dataset2, file(csv2), chrm2]
                }
                .groupTuple()
                .map { dataset, pops, count_csvs, chrms ->
                    [dataset, pops, count_csvs, chrms, 'all']
                }
        )

    emit:
        data
}

/**
 * Novel PGx Variants (Dataset-level)
 */
workflow novel_pgx_dataset {
    take:
        data

    main:
        // Identify novel sites
        novel_in_existing_database(
            data.map { dataset, vcf, csv, chrm ->
                [dataset, file(vcf), chrm]
            }
        )

        combine_novel_count5(
            novel_in_existing_database.out
                .map { dataset, vcf, count, csv, chrm ->
                    [dataset, file(count)]
                }
                .groupTuple()
                .map { dataset, counts ->
                    [dataset, counts, 'false', 'pgx']
                }
        )

        // Calculate frequencies
        get_freq(
            novel_in_existing_database.out.map { pop, vcf, count, csv, chrm ->
                [pop, file(vcf), file(count), chrm]
            }
        )

        get_novel_sites(
            get_freq.out.map { dataset, vcf, sample, sites, mafs, frq, chrm ->
                [dataset, file(vcf), file(mafs), '']
            }
        )

        // Extract novel variants
        get_gene_vcf(
            get_novel_sites.out.map { dataset, vcf, sites, chrm ->
                [dataset, file(vcf), '', chrm, '', file(sites)]
            }
        )

        // Count novel variants
        novel_count(
            get_gene_vcf.out.map { dataset, vcf, sample, chrm ->
                [dataset, file(vcf), chrm]
            }
        )

        novel_count_combine(
            novel_count.out
                .map { dataset2, csv2, chrm2 ->
                    ['AIBST', dataset2, file(csv2), chrm2]
                }
                .groupTuple()
                .map { dataset, pops, count_csvs, chrms ->
                    [dataset, pops, count_csvs, chrms, 'pgx']
                }
        )

    emit:
        data
}

/*
========================================================================================
    CONSEQUENCE ANALYSIS WORKFLOW
========================================================================================
*/

/**
 * Variant Consequence Analysis for PGx Variants
 * Predicts and summarizes functional consequences
 */
workflow csq_pgx {
    take:
        data

    main:
        csq(data)
        plot_csq(csq.out)

    emit:
        data
}

/*
========================================================================================
    POPULATION GROUPING WORKFLOWS
========================================================================================
    These workflows group populations for comparative analyses
========================================================================================
*/

/**
 * Group Populations for PGx Analysis
 * Groups populations based on dataset_groups configuration
 */
workflow group_pops_pgx {
    take:
        data

    main:
        groups_data_pgx = data.flatMap { pop, pop_vcf, pop_sample, chrm ->
            def datas = []
            params.dataset_groups.each { group, group_pops ->
                def pops = group_pops.split(',')
                if (pop in pops) {
                    datas << [group, pop, file(pop_vcf), file(pop_sample)]
                }
            }
            return datas
        }
        .groupTuple()
        .map { group, pops, pop_vcfs, pop_samples ->
            if (pops.size() > 1) {
                [group, '', pop_vcfs, pop_samples]
            }
        }

        merge_pop_groups(groups_data_pgx)

    emit:
        groups_pgx = merge_pop_groups.out
}

/**
 * Group Populations (All Chromosomes)
 * Concatenates chromosome-level VCFs for each population
 */
workflow group_pops {
    take:
        data

    main:
        pop_data = data.groupTuple(by: [0, 1])
            .map { pop, dataset, pop_vcfs, pop_samples, chrms ->
                [pop, pop_vcfs, pop_samples]
            }

        cat_chrm_groups(pop_data)

    emit:
        pop_data = cat_chrm_groups.out
}

/**
 * Group Populations by Chromosome
 * Groups populations while maintaining chromosome separation
 */
workflow group_pops_by_chrm {
    take:
        data

    main:
        groups_data_chrm = data.flatMap { pop, dataset, pop_vcf, pop_sample, chrm ->
            def datas = []
            params.dataset_groups.each { group, group_pops ->
                def pops = group_pops.split(',')
                if (pop in pops) {
                    datas << [chrm, group, pop, file(pop_vcf), file(pop_sample)]
                }
            }
            return datas
        }
        .groupTuple(by: [0, 1])
        .flatMap { chrm, group, pops, pop_vcfs, pop_samples ->
            def datas1 = []
            if (pops.size() > 1) {
                datas1 << [group, chrm, pop_vcfs, pop_samples]
            }
            return datas1
        }

        merge_pop_groups(groups_data_chrm)

    emit:
        groups_chrm = merge_pop_groups.out
}

/**
 * Create Pairwise Population Combinations
 * Generates all unique population pairs for comparative analysis
 */
workflow group_2pops {
    take:
        data

    main:
        def test_pops = Channel.from(params.FST_TEST.split(','))
        def tmp1 = []

        // Generate unique population pairs
        my_pops = test_pops.combine(test_pops)
            .flatMap { pops_2 ->
                def my_pops_ = []
                if (pops_2[0] != pops_2[1]) {
                    if ([pops_2[1], pops_2[0]] !in tmp1) {
                        tmp1 << [pops_2[0], pops_2[1]]
                        my_pops_ << ["${pops_2[0]}_${pops_2[1]}"]
                    }
                }
                return my_pops_
            }

        // Create combined VCFs for each pair
        data_2pops = data.combine(data)
        pop2_per_group_cha = data_2pops
            .map { pop1, pop1_vcf, pop1_sample, pop2, pop2_vcf, pop2_sample ->
                ["${pop1}_${pop2}", 'ALL', [file(pop1_vcf), file(pop2_vcf)], [file(pop1_sample), file(pop2_sample)]]
            }
            .combine(my_pops, by: 0)

        merge_pop_groups(pop2_per_group_cha)

    emit:
        pops2_data = pop2_per_group_cha
        pops_vcf = merge_pop_groups.out
}

/**
 * Create Pairwise Population Combinations for PGx
 * Similar to group_2pops but maintains dataset information
 */
workflow group_2pops_pgx {
    take:
        data

    main:
        def test_pops = Channel.from(params.FST_TEST.split(','))
        def tmp1 = []

        my_pops = test_pops.combine(test_pops)
            .flatMap { pops_2 ->
                def my_pops_ = []
                if (pops_2[0] != pops_2[1]) {
                    if ([pops_2[1], pops_2[0]] !in tmp1) {
                        tmp1 << [pops_2[0], pops_2[1]]
                        my_pops_ << ["${pops_2[0]}_${pops_2[1]}"]
                    }
                }
                return my_pops_
            }

        data_2pops = data.combine(data)
        pop2_per_group_cha = data_2pops
            .map { pop1, pop1_vcf, pop1_sample, dataset1, pop2, pop2_vcf, pop2_sample, dataset2 ->
                ["${pop1}_${pop2}", 'ALL', [file(pop1_vcf), file(pop2_vcf)], [file(pop1_sample), file(pop2_sample)]]
            }
            .combine(my_pops, by: 0)

        merge_pop_groups(pop2_per_group_cha)

    emit:
        pops2_data = pop2_per_group_cha
        pops_vcf = merge_pop_groups.out
}

/*
========================================================================================
    STATISTICAL ANALYSIS WORKFLOWS
========================================================================================
*/

/**
 * Fisher's Exact Test for Population Groups (PGx)
 * Tests for significant allele frequency differences between populations
 */
workflow fisher_group_pgx {
    take:
        data

    main:
        // Filter for specified test groups
        fisher_group_pgx_cha = data.flatMap { group_5, chrm_5, group_vcf_5, group_sample_5 ->
            def datas = []
            def fisher_test_group = params.FISHER_TEST.split(',')
            fisher_test_group.each { group_f ->
                if (group_f == group_5) {
                    datas << [group_5, file(group_vcf_5), file(group_sample_5), chrm_5]
                }
            }
            return datas
        }

        // Convert to PLINK format
        vcf_to_plink(fisher_group_pgx_cha)

        // Generate phenotype files
        generate_pheno_fam_sample(vcf_to_plink.out)

        // Run Fisher's exact test
        fisher_test_plink(
            generate_pheno_fam_sample.out.map { dataset, ped, bim, fam, pheno, chrm ->
                [dataset, file(ped), file(bim), file(fam), file(pheno), chrm, 'pgx']
            }
        )

    emit:
        pvalues = fisher_test_plink.out
}

/**
 * Fisher's Exact Test for Population Groups (All Variants)
 */
workflow fisher_group {
    take:
        data

    main:
        fisher_group_cha = data.flatMap { group_4, chrm_4, group_vcf_4, group_sample_4 ->
            def datas = []
            def fisher_test_group = params.FISHER_TEST.split(',')
            fisher_test_group.each { group_f ->
                if (group_f == group_4) {
                    datas << [group_4, file(group_vcf_4), file(group_sample_4), chrm_4]
                }
            }
            return datas
        }

        vcf_to_plink(fisher_group_cha)
        generate_pheno_fam_sample(vcf_to_plink.out)

        fisher_test_plink(
            generate_pheno_fam_sample.out.map { dataset, ped, bim, fam, pheno, chrm ->
                [dataset, file(ped), file(bim), file(fam), file(pheno), chrm, '']
            }
        )

        combine_csv(fisher_test_plink.out.groupTuple(by: [0]))

    emit:
        pvalues = fisher_test_plink.out
}

/**
 * Fisher's Exact Test for Pairwise Populations (PGx)
 */
workflow fisher_2pops_pgx {
    take:
        data

    main:
        vcf_to_plink(
            data.map { pops, dataset, pops_vcf, pops_samples ->
                [pops, pops_vcf, pops_samples, dataset]
            }
        )

        generate_pheno_fam_sample1(vcf_to_plink.out)

        fisher_test_plink_pops(
            generate_pheno_fam_sample1.out.map { dataset, ped, bim, fam, pheno, chrm ->
                [dataset, file(ped), file(bim), file(fam), file(pheno), chrm, 'pgx']
            }
        )

        filter_fisher_test(fisher_test_plink_pops.out)

        combine_csv(
            filter_fisher_test.out
                .groupTuple(by: [2])
                .map { pops, fisher_files, dataset ->
                    [dataset, fisher_files, pops]
                }
        )

    emit:
        data
        pvalues = combine_csv.out
}

/**
 * Fisher's Exact Test for Pairwise Populations (All Variants)
 */
workflow fisher_2pops {
    take:
        data

    main:
        vcf_to_plink(
            data.map { pops, dataset, pops_vcf, pops_samples ->
                [pops, pops_vcf, pops_samples, dataset]
            }
        )

        generate_pheno_fam_sample1(vcf_to_plink.out)

        fisher_test_plink_pops(
            generate_pheno_fam_sample1.out.map { dataset, ped, bim, fam, pheno, chrm ->
                [dataset, file(ped), file(bim), file(fam), file(pheno), chrm, 'pgx']
            }
        )

        filter_fisher_test(fisher_test_plink_pops.out)

        combine_csv(
            filter_fisher_test.out
                .groupTuple(by: [2])
                .map { pops, fisher_files, dataset ->
                    [dataset, fisher_files, pops]
                }
        )

    emit:
        data
        pvalues = combine_csv.out
}

/**
 * FST Analysis Workflow
 * Calculates population differentiation using Weir & Cockerham's FST
 */
workflow fst {
    take:
        data

    main:
        log.info "Starting FST Analysis"

        // Calculate FST for all population pairs
        fst_analysis(data)

        // Extract per-variant FST estimates
        get_fst_weir_estimates(fst_analysis.out)

        // Combine all FST results
        combine_fst_weir_estimates(
            get_fst_weir_estimates.out
                .groupTuple(by: [0])
                .map { dataset, pops, weir_fsts, mean_fsts ->
                    [dataset, weir_fsts]
                }
        )

        // Extract FST values above cutoff
        get_fst_weir_estimates_cutoff(
            fst_analysis.out.map { dataset, pops, fst_weir_estimate, fst_log ->
                [dataset, pops, fst_weir_estimate, fst_log, params.fst_cutoff]
            }
        )

        combine_fst_weir_estimates_cutoff(
            get_fst_weir_estimates_cutoff.out
                .groupTuple(by: [0])
                .map { dataset, pops, weir_fsts, mean_fsts ->
                    [dataset, weir_fsts]
                }
        )

        // Calculate mean FST per population
        combine_fst_mean_weir_estimates_1(
            get_fst_weir_estimates.out
                .map { dataset, pops, weir_fst, mean_fst ->
                    [pops.split('__')[0], weir_fst]
                }
                .mix(
                    get_fst_weir_estimates.out.map { dataset, pops, weir_fst, mean_fst ->
                        [pops.split('__')[1], weir_fst]
                    }
                )
                .groupTuple(by: [0])
        )

        // Calculate overall mean FST
        combine_fst_mean_weir_estimates(
            get_fst_weir_estimates.out
                .groupTuple(by: [0])
                .map { dataset, pops, weir_fst, mean_fst ->
                    [dataset, mean_fst]
                }
        )

        // Generate FST matrix for visualization
        generate_fst_matrix(combine_fst_mean_weir_estimates.out)

        log.info "FST Analysis Completed"

    emit:
        data
        fst_2pops_all = combine_fst_weir_estimates.out
        fst_2pops = combine_fst_mean_weir_estimates_1.out
}

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
    Orchestrates all analysis workflows in the correct order
========================================================================================
*/

workflow {
    main:
        log.info ""
        log.info "="*80
        log.info "  African Exome Analysis Workflow - 2025 Edition"
        log.info "="*80
        log.info "  Project: ${workflow.projectDir}"
        log.info "  Command: ${workflow.commandLine}"
        log.info "  Profile: ${workflow.profile}"
        log.info "  Work Dir: ${workflow.workDir}"
        log.info "  Output: ${params.outDir}"
        log.info "="*80
        log.info ""

        // ===== PHASE 1: Core Data Processing =====
        log.info "PHASE 1: Core Data Processing"
        data_proc()

        data_proc_pop(data_proc.out.fill_tags_dataset_chrm)

        // ===== PHASE 2: Population Grouping =====
        log.info "PHASE 2: Population Grouping"

        // Group by chromosome for QC and general analysis
        group_pops_by_chrm(data_proc_pop.out.pop_exome_chrm)

        // Group all chromosomes per population
        group_pops(data_proc_pop.out.pop_exome_chrm)

        // Create pairwise combinations for FST
        group_2pops(group_pops.out.pop_data)

        // Group PGx variants
        group_pops_pgx(data_proc_pop.out.pop_pgx)
        group_2pops_pgx(data_proc_pop.out.pop_pgx)

        // ===== PHASE 3: Statistical Analyses =====
        log.info "PHASE 3: Statistical Analyses"

        // Fisher's Exact Test
        fisher_group_pgx(group_pops_pgx.out.groups_pgx)
        fisher_2pops_pgx(group_2pops_pgx.out.pops_vcf)

        // FST Analysis
        fst(group_2pops.out.pops2_data)

        // ===== PHASE 4: HDV Detection =====
        log.info "PHASE 4: Highly Differentiated Variants Detection"

        // Combine Fisher and FST results for comprehensive HDV detection
        hdv_pgx_pop_cha = data_proc_pop.out.pop_pgx
            .map { pop, pop_vcf, pop_sample, dataset ->
                [pop, file(pop_vcf), dataset, '']
            }
            .combine(fisher_2pops_pgx.out.pvalues)
            .combine(fst.out.fst_2pops, by: [0])
            .flatMap { pop, pop_vcf_chrm, dataset, chrm, group, pvalue_file, fst_file ->
                def datas = []
                def pops = params.FST_TEST.split(',')
                if (pop in pops) {
                    datas << [pop, file(pop_vcf_chrm), group, chrm, file(pvalue_file), file(fst_file)]
                }
                return datas
            }

        hdv_pgx_pop(hdv_pgx_pop_cha, data_proc.out.pop_mafs_vcf)

        // ===== PHASE 5: Consequence Analysis =====
        log.info "PHASE 5: Variant Consequence Analysis"

        csq_data = data_proc.out.pgx_dataset
            .combine([file(params.so_term_file)])
            .flatMap { dataset3, vcf3, sample3, so_term ->
                def datas3 = []
                if (dataset3 == 'AIBST') {
                    datas3 << [dataset3, file(vcf3), file(sample3), file(so_term)]
                }
                return datas3
            }
        csq_pgx(csq_data)

        // ===== PHASE 6: Optional Analyses =====
        log.info "PHASE 6: Optional Analyses"

        // PCA Analysis
        pca(group_pops_by_chrm.out.groups_chrm)
        // ADMIXTURE Analysis (conditional)
        if (params.run_admixture) {
            log.info "  Running ADMIXTURE (K=" + params.admixture_k_min + "-" + params.admixture_k_max + ")"
            // Run admixture only on AIBST (all 127 samples), not AIBST_NO_ZWD
            admixture_input = pca.out.plink_data
                .filter { dataset, ped, map, sample -> dataset == 'AIBST' }
            admixture(admixture_input)
        }
        // Selection Scans (conditional)
        if (params.run_selection_scans) {
            log.info "  Running Selection Scans (iHS, XP-EHH, Tajima's D)"
            // Prepare input: [pop, chr, vcf]
            selection_input = data_proc_pop.out.pop_exome_chrm
                .map { pop, chr, vcf ->
                    tuple(pop, chr, file(vcf))
                }
            selection_scans(selection_input)
        }

        // Variant Counting
        count_all_dataset(data_proc.out.fill_tags_dataset_chrm.combine(data_proc.out.singl_dataset_chrm, by:[0,3]))
        count_pgx_dataset(data_proc.out.pgx_dataset_chrm.combine(data_proc.out.singl_dataset_chrm, by:[0,3]))

        // ===== PHASE 7: Pharmacogenomics from BAMs =====
        // Star allele calling using existing joint VCF + BAMs for depth
        if (params.run_pgx_calling || params.run_hla_typing) {
            log.info "PHASE 7: Pharmacogenomics (Star Alleles + HLA Typing)"

            // Read BAMs from sarek work directory (unique per sample)
            bam_ch = Channel
                .fromFilePairs("${params.sarek_bam_dir}/**/*.recal.{bam,bam.bai}", size: 2, flat: true)
                .map { prefix, bam, bai ->
                    def sample_id = bam.name.replaceAll(/\.recal\.bam$/, "")
                    [sample_id, bam, bai]
                }
                .unique { it[0] }  // One BAM per sample

            if (params.run_pgx_calling) {
                log.info "  Running PyPGx star allele calling (joint VCF + BAM depth)"
                def joint_vcf = file(params.joint_vcf)
                def joint_vcf_tbi = file("${params.joint_vcf}.tbi")
                def fasta = file(params.fasta)
                def fasta_fai = file("${params.fasta}.fai")
                pgx_star_alleles(joint_vcf, joint_vcf_tbi, bam_ch, fasta, fasta_fai)
            }

            if (params.run_hla_typing) {
                log.info "  Running HLA typing (OptiType + arcasHLA)"
                hla_typing(bam_ch)
            }
        }

        log.info ""
        log.info "="*80
        log.info "  Workflow Completed Successfully!"
        log.info "  Results: ${params.outDir}"
        log.info "="*80
        log.info ""
}

/*
========================================================================================
    WORKFLOW COMPLETION HANDLER
========================================================================================
*/

workflow.onComplete {
    log.info ""
    log.info "="*80
    log.info "  Pipeline execution summary"
    log.info "="*80
    log.info "  Completed at: ${workflow.complete}"
    log.info "  Duration    : ${workflow.duration}"
    log.info "  Success     : ${workflow.success}"
    log.info "  Exit status : ${workflow.exitStatus}"
    log.info "  Error report: ${workflow.errorReport ?: 'None'}"
    log.info "="*80
    log.info ""
}

workflow.onError {
    log.error ""
    log.error "="*80
    log.error "  Pipeline execution failed!"
    log.error "="*80
    log.error "  Error: ${workflow.errorMessage}"
    log.error "  Error report: ${workflow.errorReport}"
    log.error "="*80
    log.error ""
}

/*
========================================================================================
    END OF WORKFLOW
========================================================================================
*/
