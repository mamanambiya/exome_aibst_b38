#!/usr/bin/env nextflow
nextflow.enable.dsl=2

"""
Author: Mamana M.
Affiliation: University of Cape Town
Aim: A simple nextflow workflow to download and process 1000 Genomes project data.
Date: Mon June 8 14:03:11 CET 2015
Run: 

Latest modification:
  - TODO
"""

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

//// Consequences

process csq {
    tag "csq_${dataset}"
    label "medium"
    // publishDir "${params.work_dir}/data/PGX_ONLY/CSQ/${dataset}", mode:'symlink'
    // publishDir "${params.work_dir}/REPORTS/${dataset}/PGX_ONLY/CSQ", mode:'copy', pattern: "*_summaryVariantsPerGene.tsv"
    publishDir "${params.outDir}/csq/pgx_only/${dataset}", mode:'copy', overwrite: true

    input:
        tuple val(dataset), file(dataset_vcf), file(dataset_sample), file(so_term_file)

    output:
        tuple val(dataset), file(dataset_sample), file("${dataset}.singletons.vcf.gz"), file("${dataset}_SO.terms.MAF.summary"), file("${dataset}.singletons.uniqueID"), file("${dataset}.nonsing-0.01MAF.uniqueID"), file("${dataset}.0.01-0.05MAF.uniqueID"), file("${dataset}.greater.0.05.uniqueID")

    script:
        """
        # Create different vcfs for different frequency classes
        # Singletons
        vcftools --gzvcf ${dataset_vcf} --max-mac 1 --recode-INFO-all --recode --stdout | \
        bgzip -c > ${dataset}.singletons.vcf.gz
        # Non-singletons to MAF 0.01
        vcftools --gzvcf ${dataset_vcf} --recode --mac 2 --max-maf 0.01 --recode-INFO-all --stdout | \
        bgzip -c > ${dataset}.nonsing-0.01MAF.vcf.gz
        # MAF 0.01-0.05
        vcftools --gzvcf ${dataset_vcf} --recode --maf 0.01 --max-maf 0.05 --recode-INFO-all --stdout | \
        bgzip -c > ${dataset}.0.01-0.05MAF.vcf.gz
        # MAF > 0.05
        vcftools --gzvcf ${dataset_vcf} --recode --maf 0.05 --recode-INFO-all --stdout | \
        bgzip -c > ${dataset}.greater.0.05MAF.vcf.gz
        # Count number of appearances of each of the Sequence Ontology (SO) terms
        while read SO; do
            zcat ${dataset_vcf} | grep \$SO | wc -l
        done < ${so_term_file} > ${dataset}_allCount
        while read SO; do
            zcat ${dataset}.singletons.vcf.gz | grep \$SO | wc -l
        done < ${so_term_file} > ${dataset}_singletonsCount
        while read SO; do
            zcat ${dataset}.nonsing-0.01MAF.vcf.gz | grep \$SO | wc -l
        done < ${so_term_file} > ${dataset}_nonsing-0.01MAFCount
        while read SO; do
            zcat ${dataset}.0.01-0.05MAF.vcf.gz | grep \$SO | wc -l
        done < ${so_term_file} > ${dataset}_0.01-0.05MAFCount
        while read SO; do
            zcat ${dataset}.greater.0.05MAF.vcf.gz | grep \$SO | wc -l
        done < ${so_term_file} > ${dataset}_greater0.05MAFCount
        echo "terms 	countAll	countSingletons	countnonSing-0.01MAF	count0.01-0.05MAF	countGreater0.05" > ${dataset}_SO.terms.MAF.summary
        paste ${so_term_file} ${dataset}_allCount ${dataset}_singletonsCount  ${dataset}_nonsing-0.01MAFCount ${dataset}_0.01-0.05MAFCount ${dataset}_greater0.05MAFCount >> ${dataset}_SO.terms.MAF.summary
        # Create UniqueID for each of the vcf frequency classes
        zcat ${dataset}.singletons.recode.vcf.gz | grep -v ^# | awk '{print \$1":"\$2":"\$5}' >  ${dataset}.singletons.uniqueID
        zcat ${dataset}.nonsing-0.01MAF.recode.vcf.gz | grep -v ^# | awk '{print \$1":"\$2":"\$5}' >  ${dataset}.nonsing-0.01MAF.uniqueID
        zcat ${dataset}.0.01-0.05MAF.recode.vcf.gz | grep -v ^# | awk '{print \$1":"\$2":"\$5}' >  ${dataset}.0.01-0.05MAF.uniqueID
        zcat ${dataset}.greater.0.05MAF.recode.vcf.gz | grep -v ^# | awk '{print \$1":"\$2":"\$5}' >  ${dataset}.greater.0.05.uniqueID
        """
}

//// Step 16.2: Plot consequences

process plot_csq {
    tag "plot_csq_${dataset}"
    label "medium"
    // publishDir "${params.work_dir}/data/PGX_ONLY/CSQ/${dataset}", mode:'symlink'
    // publishDir "${params.outDir}/REPORTS/${dataset}/PGX_ONLY/CSQ", mode:'copy', overwrite: true
    publishDir "${params.outDir}/csq/pgx_only/${dataset}", mode:'copy', overwrite: true

    input:
        tuple val(dataset), file(dataset_sample), file(dataset_singletons_vcf_gz), file(dataset_SO_terms_MAF_summary), file(dataset_singletons_uniqueID), file(dataset_nonsing_0_01MAF_uniqueID), file(dataset_0_01_0_05MAF_uniqueID), file(dataset_greater_0_05_uniqueID)
    output:
        tuple val(dataset), file(dataset_pgxFunctionalClassesCounts_tiff), file(dataset_pgxFunctionalClassesByFrequency_tiff)
    script:
        dataset_SO_terms_MAF_summary = dataset_SO_terms_MAF_summary
        dataset_pgxFunctionalClassesCounts_tiff = "${dataset}_pgxFunctionalClassesCounts.pdf"
        dataset_pgxFunctionalClassesByFrequency_tiff = "${dataset}_pgxFunctionalClassesByFrequency.pdf"
        template "step16_2_consequence.R"
    }