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

process dataset_qc_alt {
    tag "target_qc_${dataset_name}_${chrm}"
    label "medium"

    input:
        tuple val(dataset_name), file(dataset_vcf), file(dataset_sample), val(chrm)
    
    output:
        tuple val(dataset_name), file(vcf_out), file(dataset_sample), val(chrm)
    
    script:
        base = file(dataset_vcf.baseName).baseName
        vcf_out = "${base}_clean.vcf.gz"
        """
        # Note: Source VCFs are pre-cleaned (_fixed.hg38.vcf.gz) with invalid 1KG_MAF tag removed
        # Remove problematic FORMAT fields that have wrong number of values for multiallelic sites
        # Split multiallelic sites with bcftools norm
        # Filter out invalid genotype indices (GT~"2") that cause SnpEff to fail with header corruption
        #
        # After bcftools norm splits multiallelic sites, some genotypes reference allele indices
        # that no longer exist (e.g., GT=0/2 when only 1 ALT allele exists after split).
        # These MUST be removed BEFORE SnpEff annotation to prevent corrupted output.
        bcftools annotate -x FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/PL,FORMAT/DP ${dataset_vcf} -Ou |\
        bcftools norm \
            --check-ref x \
            -m-any \
            --fasta-ref ${params.ref_genome} \
            -Ou |\
        bcftools view -e 'ALT="." || GT~"2"' -Oz -o ${vcf_out}
        """
}

process get_vcf_site {
    tag "get_vcf_site_${dataset}_${chrm}"
    label "medium"

    input:
        tuple val(dataset), file(vcf), file(sample), val(chrm)
    
    output:
        tuple val(dataset), file(vcf), file(sample), file(vcf_sites), val(chrm)
    
    script:
        base = file(vcf.baseName).baseName
        vcf_sites = "${base}.sites"
        """
        echo -e '#CHRM:POS' > ${vcf_sites}
        bcftools query \
        -f "%CHROM:%POS\\n" \
        ${vcf} >> ${vcf_sites}
        """
}

process sites_only {
    tag "sites_only_${dataset}"
    label "bigmem"

    input:
        tuple val(dataset), file(vcf), file(sample)

    output:
        tuple val(dataset), file(sites_vcf), file(sample)

    script:
        sites_vcf = "${vcf.getSimpleName()}_sites.bcf"
        """
        tabix ${vcf}
        bcftools view ${vcf} --drop-genotypes --threads ${task.cpus} -Ob -o ${sites_vcf}
        tabix ${sites_vcf}
        """
}

process tabix_data_proc_pop {
    tag "target_qc_${pop}_${group}"
    label "medium"

    input:
        tuple val(pop), file(pop_vcf), file(pop_sample), val(group)
    
    output:
        tuple val(pop), file(pop_vcf), file(pop_sample), val(group)
    
    script:
        """
        tabix ${pop_vcf}
        """
}