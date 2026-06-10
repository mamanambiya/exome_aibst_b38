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

process pruning_vcf1 {
    tag "pruning_${dataset}"
    label "bigmem"

    input:
        tuple val(dataset), file(dataset_vcf), file(dataset_sample)
    output:
        tuple val(dataset), file("${base}.ped"), file("${base}.map"), file(dataset_sample)
    script:
        base = dataset_vcf.getSimpleName()
        """
        tabix -f ${dataset_vcf}
        plink2 \
                --vcf ${dataset_vcf} \
                --indep-pairwise 50 5 0.5 \
                --allow-no-sex \
                --make-bed \
                --snps-only --biallelic-only strict \
                --out ${base}
        plink2 \
            --vcf ${dataset_vcf} \
            --extract ${base}.prune.in \
            --make-bed --snps-only --biallelic-only strict \
            --out ${base}
        """
}

process pruning_vcf {
    tag "pruning_${dataset}"
    label "bigmem"

    input:
        tuple val(dataset), file(dataset_vcf), file(dataset_sample)
    output:
        tuple val(dataset), file(vcf_pruned), file(dataset_sample)
    script:
        base = dataset_vcf.getSimpleName()
        vcf_pruned = "${base}_pruned.bcf"
        """
        tabix -f ${dataset_vcf}
        bcftools norm -c x -f ${params.ref_genome} ${dataset_vcf} |\
        bcftools view -m2 -M2 -v snps |\
        bcftools annotate --set-id '%CHROM\\_%POS' | \
        bcftools +prune -l 0.5 -w 10kb -n 10 -e'F_MISSING>=0.02' -Ob -o ${vcf_pruned}
        """
}

process vcf_to_plink1 {
    tag "vcf_to_plink_${dataset}"
    label "bigmem"

    input:
        tuple val(dataset), file(dataset_vcf), file(dataset_sample)
    output:
        tuple val(dataset), file("${base}.ped"), file("${base}.map"), file(dataset_sample)
    script:
        base = dataset_vcf.getSimpleName()
        """
        plink --vcf ${dataset_vcf} --vcf-half-call missing --recode --out ${base}
        """
}

// reduce_sample id size to less than 20 and set pheno (col 6) of PED to 1 or it get ignored by smartpca
process updated_ped {
    tag "reduce_sample_id_size_${dataset}"
    label "bigmem"

    input:
        tuple val(dataset), file(dataset_ped), file(dataset_map), file(dataset_sample)
    output:
        tuple val(dataset), file(new_dataset_ped), file(dataset_map), file(dataset_sample)
    script:
        new_dataset_ped = "${dataset_ped.baseName}_new.ped"
        in_ped = dataset_ped
        out_ped = new_dataset_ped
        template "reduce_sample_id_size.py"       
}

process smartpca_dataset {
    tag "smartpca_dataset_${dataset}"
    label "smartpca"

    input:
        tuple val(dataset), file(dataset_ped), file(dataset_map), file(dataset_sample)
    output:
        tuple val(dataset), file(dataset_evec), file(dataset_eval), file(dataset_grmjunk), file(dataset_sample)
    script:
        base = dataset_ped.getSimpleName()
        dataset_evec = "${base}.evec"
        dataset_eval = "${base}.eval"
        dataset_grmjunk = "${base}.evec_grmjunk"
        """
        ## Create parameter file for smartpca
        echo -e \
        "genotypename:   ${dataset_ped}
        snpname:         ${dataset_map}
        indivname:       ${dataset_ped}
        evecoutname:     ${dataset_evec}
        evaloutname:     ${dataset_eval}
        altnormstyle:    NO
        numoutevec:      5
        numoutlieriter:  5
        familynames:     NO
        outliermode:     2
        numthreads:      ${task.cpus}
        grmoutname:      ${dataset_grmjunk}" > ${dataset}.EIGENSTRAT.par
        ## Run smartpca
        smartpca \
                -p ${dataset}.EIGENSTRAT.par \
                > ${dataset}.EIGENSTRAT.log
        """
}

process update_evec {
    tag "update_evec_${group_name}"
    label "bigmem"
    input:
        tuple val(group_name), file(group_evec), file(group_eval), file(group_grmjunk), file(group_sample)
    output:
        tuple val(group_name), file(group_evec_update), file(group_eval), file(group_grmjunk), file(group_sample)
    script:
        group_evec_update = "${file(group_evec).baseName}_update.evec"
        evec_file = group_evec
        evec_out = group_evec_update
        annot_file = group_sample
        template "update_evec.py"
}

"""
Step: Plot PCA for group
"""
process plot_pca_group {
    tag "plot_pca_group_${group_name}"
    label "rplot"
    publishDir "${params.outDir}/pca", overwrite: true, mode:'copy', pattern: '*tiff'
    
    input:
        tuple val(group_name), file(group_evec), file(group_eval), file(group_grmjunk), file(group_sample)
    output:
        tuple val(group_name), file(group_evec), file(group_sample), file("${output_pdf}*tiff")
    script:
        output_pdf = "${group_name}"
        input_evec = group_evec
        template "plot_pca.R"
}