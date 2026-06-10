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

// def add_pop_groups(it1, dataset_files){ 
//     add_pop_groups_datas = []
//     dataset = it1[0]
//     dataset_vcf = it1[1]
//     dataset_sample = it1[2]
//     chrm = it1[3]
//     dataset_files.each{ dataset_name, dataset_vcf_, dataset_sample_, dataset_pops ->
//         dataset_pops.split(',').each{ pop ->
//             if(dataset == dataset_name){
//                 add_pop_groups_datas = [pop, dataset, file(dataset_vcf), file(dataset_sample), chrm]
//             }
//         }
//     }
//     return add_pop_groups_datas
// }

process concat_mafs {
    tag "concat_mafs_${dataset}"
    label "medium"
    
    input:
        tuple val(dataset), val(vcfs), val(samples), val(sites), val(pops), val(mafs)
    
    output:
        tuple val(dataset), file(vcf), file(sample), file(site), val(pops), val(mafs)
    
    script:
        vcf = "${dataset}.vcf.gz"
        sample = "${dataset}.sample"
        site = "${dataset}.sites"
        """
        head -n 1 ${sites[0]} > ${site}
        tail -q -n +2 ${sites.join(' ')} >> ${site}
        ln -s ${vcfs[0]} ${vcf}
        ln -s ${samples[0]} ${sample}
        """
}

process combine_mafs_dataset {
    tag "mafs_${dataset_name}_${chrm}_${mafs_dataset}"
    label "extrabig"
    
    input:
        tuple val(dataset_name), file(dataset_vcf), file(dataset_sample), file(sites), val(chrm), val(mafs_dataset), val(mafs_files)
    
    output:
        tuple val(dataset_name), file(dataset_vcf), file(dataset_sample), val(chrm), file(outAnnot), file(outHdr)
    
    script:
        inTSV = mafs_files
        outAnnot = "${file(dataset_vcf.baseName).baseName}_mafs.tsv"
        outHdr = "${file(dataset_vcf.baseName).baseName}_mafs.hdr"
        annot = "True"
        template "freq_annotation.py"
}

process fix_mafs_headers {
    tag "fix_mafs_${dataset_name}_${chrm}"
    label "extrabig"
    
    input:
        tuple val(dataset_name), file(dataset_vcf), file(dataset_sample), val(chrm), file(annot), file(hdr)
    
    output:
        tuple val(dataset_name), file(dataset_vcf), file(dataset_sample), val(chrm), file(outAnnot), file(outHdr)
    
    script:
        outAnnot = "${annot.baseName}_fixed.tsv"
        outHdr = "${hdr.baseName}_fixed.hdr"
        template "fix_mafs_headers.py"
}

process annotate_mafs {
    tag "mafs_${dataset_name}_${chrm}"
    label "medium"
    // maxForks 10  // removed to allow full parallelism
    
    input:
        tuple val(dataset_name), file(dataset_vcf), file(dataset_sample), val(chrm), file(maf_annot), file(header_annot)
    
    output:
        tuple val(dataset_name), file(outVCF), file(dataset_sample), val(chrm), val(snpeff_database)
    
    script:
        base = "${dataset_name}_${maf_annot.baseName}_${chrm}"
        outVCF = "${base}.vcf.gz"
        """
        columns=\$(awk 'NR==1{print \$0}' ${maf_annot} | sed -e 's/\\s\\+/,/g')
        tail -n+2 ${maf_annot} | sort -k1,1n -k2,2n -V > ${maf_annot}.annot.tsv
        bgzip ${maf_annot}.annot.tsv
        tabix -s1 -b2 -e3 ${maf_annot}.annot.tsv.gz
        bcftools sort ${dataset_vcf} -T . -Oz -o ${dataset_vcf}.sorted.vcf.gz
        tabix ${dataset_vcf}.sorted.vcf.gz
        bcftools annotate -a ${maf_annot}.annot.tsv.gz -h ${header_annot} -c \${columns} -Oz -o ${outVCF} ${dataset_vcf}.sorted.vcf.gz
        rm ${dataset_vcf}.sorted.vcf.gz
        tabix ${outVCF}
    """
}

process annotate_vcf_from_vcf {
    tag "mafs_${dataset_name}_${annot_name}_${chrm}"
    label "medium"
    
    input:
        tuple val(chrm), val(dataset_name), file(dataset_vcf), file(dataset_sample), val(annot_name), file(annot_vcf), file(annot_vcf_idx)
    
    output:
        tuple val(dataset_name), file(outVCF), file(dataset_sample), val(chrm), val(snpeff_database)
    
    script:
        base = "${dataset_name}_${annot_name}_${chrm}"
        outVCF = "${base}.vcf.gz"
        """
        tabix ${dataset_vcf}
        bcftools annotate -a ${annot_vcf} -c INFO -Oz -o ${outVCF} ${dataset_vcf}
        """
}


/*
    Annotate VCF with gnomAD v4.1 allele frequencies.
    Adds global, exome, genome, and population-specific AFs.
    Uses -R regions file (padded WES BED) for speed.
*/
process annotate_gnomad_afr {
    tag "gnomad_afr_${dataset_name}_${chrm}"
    label "medium"
    cpus 4

    input:
        tuple val(dataset_name), file(dataset_vcf), file(dataset_sample), val(chrm), val(snpeff_database), file(gnomad_vcf), file(gnomad_vcf_idx)

    output:
        tuple val(dataset_name), file(outVCF), file(dataset_sample), val(chrm), val(snpeff_database)

    script:
        base = "${dataset_name}_gnomad_afr_${chrm}"
        outVCF = "${base}.vcf.gz"
        regions_arg = params.wes_regions_bed ? "-R ${params.wes_regions_bed}" : ""
        """
        tabix -f ${dataset_vcf}
        bcftools annotate \
            -a ${gnomad_vcf} \
            -c INFO/AF_joint,INFO/AF_exomes,INFO/AF_genomes,INFO/AF_joint_afr,INFO/AF_joint_amr,INFO/AF_joint_eas,INFO/AF_joint_nfe,INFO/AF_joint_sas \
            ${regions_arg} \
            --threads ${task.cpus} \
            -Oz -o ${outVCF} \
            ${dataset_vcf}
        tabix ${outVCF}
        """
}


/*
    Annotate vcf file with snpEff
*/
process annotate_snpeff{
    tag "snpeff_${dataset_name}_${chrm}"
    label "medium"
    cpus 4
    publishDir "${params.outDir}/snpeff_stats/${dataset_name}", pattern: "*.{html,csv}", mode: 'copy'

    input:
        tuple val(dataset_name), file(vcf_file), file(dataset_sample), val(chrm), val(snpeff_human_db), val(snpeff_database)

    output:
        tuple val(dataset_name), file(vcf_out), file(dataset_sample), val(chrm), val(snpeff_database)

    script:
        base = "${file(vcf_file.baseName).baseName}"
        vcf_out = "${base}_snpeff.vcf.gz"
        mem_gb = task.memory ? task.memory.toGiga() : 18
        """
        set -eo pipefail
        tabix ${vcf_file}
        snpEff -Xmx${mem_gb}g \
            ${snpeff_human_db} -lof \
            -stats ${base}_snpeff.html \
            -csvStats ${base}_snpeff.csv \
            -dataDir ${snpeff_database} \
            ${vcf_file} -v | \
        bgzip -@ ${task.cpus} -c  > ${vcf_out}
        tabix ${vcf_out}
        """
}

process annotate_dbsnp{
    tag "dbsnp_${dataset_name}_${chrm}"
    label "medium"
    cpus 4
    publishDir "${params.data_dir}/data/AIBST/VCF/", mode: 'copy'

    input:
    tuple val(dataset_name), file(vcf_file), file(dataset_sample), val(chrm), val(snpeff_database), file(dbsnp_vcf), file(dbsnp_vcf_tbi)

    output:
    tuple val(dataset_name), file(vcf_out), file(dataset_sample), val(chrm), val(snpeff_database)

    script:
    vcf_out = "${file(vcf_file.baseName).baseName}_dbsnp.vcf.gz"
    """
    set -eo pipefail
    tabix ${vcf_file}
    SnpSift \
        annotate \
        ${dbsnp_vcf} \
        -dataDir ${snpeff_database} \
        ${vcf_file} | \
    bgzip -@ ${task.cpus} -c > ${vcf_out}
    tabix ${vcf_out}
    """
}

/*
Step 7: Annotate with snpEff using clinvar
*/
process annotate_clinvar {
    tag "clinvar_${dataset_name}_${chrm}"
    label "medium"
    cpus 4
    publishDir "${params.data_dir}/data/AIBST/VCF/", pattern: "*_clinvar.vcf.gz*", mode: 'copy'

    input:
    tuple val(dataset_name), file(vcf_file), file(dataset_sample), val(chrm), val(snpeff_database), file(clinvar), file("${clinvar}.tbi")

    output:
    tuple val(dataset_name), file(vcf_out), file(dataset_sample), val(chrm), val(snpeff_database)

    script:
    vcf_out = "${file(vcf_file.baseName).baseName}_clinvar.vcf.gz"
    """
    set -eo pipefail
    tabix ${vcf_file}
    SnpSift \
        annotate \
        ${clinvar} \
        -dataDir ${snpeff_database} \
        ${vcf_file} | \
    bgzip -@ ${task.cpus} -c > ${vcf_out}
    tabix ${vcf_out}
    """
}



/*
Annotate with snpEff using cosmic
*/
process annotate_cosmic {
    tag "cosmic_${dataset_name}_${chrm}"
    label "small"
    publishDir "${params.data_dir}/data/AIBST/VCF/", pattern: "*_cosmic.vcf.gz*", mode: 'copy'

    input:
    tuple val(dataset_name), file(vcf_file), file(dataset_sample), val(chrm), val(snpeff_database), file(cosmic), file("${cosmic}.tbi")

    output:
    tuple val(dataset_name), file(vcf_out), file(dataset_sample), val(chrm)

    script:
    vcf_out = "${file(vcf_file.baseName).baseName}_cosmic.vcf.gz"
    """
    set -eo pipefail
    SnpSift \
        annotate \
        ${cosmic} \
        -dataDir ${snpeff_database} \
        ${vcf_file} |
    bgzip -c > ${vcf_out}
    tabix ${vcf_out}
    """
}

process filter_empty_vcf {
    tag "filter_empty_${dataset_name}_${chrm}"
    label "small"
    errorStrategy 'ignore'

    input:
        tuple val(dataset_name), file(vcf_file), file(dataset_sample), val(chrm), val(snpeff_database)

    output:
        tuple val(dataset_name), file(vcf_file), file(dataset_sample), val(chrm), val(snpeff_database), optional: true

    script:
        """
        # Check if VCF has any variant lines (non-header lines)
        VARIANT_COUNT=\$(zcat ${vcf_file} | grep -v '^#' | wc -l)

        if [ "\$VARIANT_COUNT" -eq 0 ]; then
            echo "VCF is empty (only headers), skipping..."
            exit 1
        fi
        """
}

process clean_invalid_variants {
    tag "clean_invalid_${dataset_name}_${chrm}"
    label "small"

    input:
        tuple val(dataset_name), file(vcf_file), file(dataset_sample), val(chrm), val(snpeff_database)

    output:
        tuple val(dataset_name), file(vcf_out), file(dataset_sample), val(chrm), val(snpeff_database)

    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_cleaned.vcf.gz"
        """
        # Remove malformed variants from liftover that cause bcftools +fill-tags to fail:
        # - ALT="." (monomorphic sites)
        # - Genotypes referencing allele index >= 2 (malformed - more alleles in GT than in ALT field)
        bcftools view -e 'ALT="." || GT~"2"' ${vcf_file} -Oz -o ${vcf_out}
        """
}

process fill_tags {
    tag "fill_tags_${dataset_name}_${chrm}"
    label "small"
    cpus 4

    input:
        tuple val(dataset_name), file(vcf_file), file(dataset_sample), val(chrm), val(snpeff_database)

    output:
        tuple val(dataset_name), file(vcf_out), file(dataset_sample), val(chrm)
        // tuple val(dataset_name), file(vcf_out), val(chrm)

    script:
        vcf_out = "${file(vcf_file.baseName).baseName}_maf.vcf.gz"
        """
        bcftools +fill-tags ${vcf_file} --threads ${task.cpus} -Oz -o ${vcf_out}
        """
}

process get_pop {
    tag "get_pop_${dataset}"
    label "small"

    input:
        tuple val(dataset), file(vcf_out), file(dataset_sample), val(chrm)
    output:
        tuple val(dataset), file(vcf_out), file(dataset_sample_pop), val(chrm)
    script:
        dataset_sample_pop = "${dataset_sample.getSimpleName()}.sample.pop"
        """
        tail -q -n+2 ${dataset_sample} | sed -e 's/ /\\t/g' | awk '{print \$1}' > ${dataset_sample_pop}
        """
}

// Concatenate chromosome VCFs by population
process concat_vcf {
    tag "concat_dataset_${dataset}"
    label "medium"
    
    input:
        tuple val(dataset), val(vcfs), val(samples), val(chrms)
    
    output:
        tuple val(dataset), file(dataset_vcf), val(dataset_sample)
    
    script:
        dataset_sample = file(samples[0])
        dataset_vcf = "${dataset}.vcf.gz"
        """
        bcftools concat ${vcfs.join(' ')} |\
        bcftools sort -T . -Oz -o ${dataset_vcf}
        tabix ${dataset_vcf}
        """
}

process concat_chunks_vcf {
    tag "concat_dataset_${dataset}_${chrm}"
    label "medium"

    input:
        tuple val(dataset), val(vcfs), val(samples), val(chrm)

    output:
        tuple val(dataset), file(dataset_vcf), file(dataset_sample), val(chrm)

    script:
        dataset_sample = "${dataset}_samples.csv"
        dataset_vcf = "${dataset}.vcf.gz"
        """
        # Create symlink to sample file in work directory
        ln -s ${samples[0]} ${dataset_sample}

        bcftools concat ${vcfs.join(' ')} |\
        bcftools sort -T . -Oz -o ${dataset_vcf}
        tabix ${dataset_vcf}
        """
}
process reheader_samples {
    tag "reheader_${dataset}"
    label "small"
    publishDir "${params.outDir}/reheadered", mode: "symlink"

    input:
        tuple val(dataset), path(vcf), path(vcf_idx), path(sample_file)

    output:
        tuple val(dataset), path("${dataset}_renamed.vcf.gz"), path("${dataset}_renamed.vcf.gz.tbi"), path(sample_file)

    script:
        """
        bcftools reheader -s ${params.sample_rename_file} ${vcf} -o ${dataset}_renamed.vcf.gz
        tabix -p vcf ${dataset}_renamed.vcf.gz
        """
}
