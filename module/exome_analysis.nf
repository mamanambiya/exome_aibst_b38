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

scratch_dir = params.scratch_dir

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

// check if files exist [name, file1, file2, ...]
def check_files(file_list) {
    file_list.each { myfile ->
        if (!file(myfile).exists() && !file(myfile).isFile()) exit 1, "|-- ERROR: File ${myfile} not found. Please check your config file."
    }
}


// Step: Split dataset per population
process split_dataset_vcf_pop {
    tag "split_dataset_${dataset}_${pop}_${chrm}"
    label "medmem"

    input:
        tuple val(pop), val(dataset), file(dataset_chrm_vcf), file(dataset_sample), val(chrm)
    
    output:
        tuple val(pop), val(dataset), file(pop_chrm_vcf), file(pop_sample_update), val(chrm)
    
    script:
        pop_chrm_vcf = "${pop}_${dataset}_${chrm}.vcf.gz"
        pop_sample = "${pop}.sample"
        pop_sample_update = "${pop}_update.sample"
        """
        awk ' \$2==\"${pop}\" { \$2=\$2"\\t${dataset}"; print \$0}' ${dataset_sample} > ${pop_sample_update}
        grep ${pop} ${dataset_sample} | awk '{print \$1}' > ${pop_sample}
        ## Keep only samples for population and Recalculate AC, AN, AF
        tabix ${dataset_chrm_vcf}
        bcftools view --samples-file ${pop_sample} --min-ac 1:minor ${dataset_chrm_vcf} | \
        bcftools +fill-tags | \
        bgzip -c > ${pop_chrm_vcf}
        tabix ${pop_chrm_vcf}
        """
}

// Step: Split dataset per population
process split_dataset_vcf_pop_tabix {
    tag "split_dataset_${dataset}_${pop}_${chrm}"
    label "medmem"

    input:
        tuple val(pop), val(dataset), file(dataset_chrm_vcf), file(dataset_sample), val(chrm)
    
    output:
        tuple val(pop), val(dataset), file(dataset_chrm_vcf), file(dataset_sample), val(chrm)
    
    script:
        """
        tabix ${dataset_chrm_vcf}
        """
}


// Step: Get singleton for each population
process singleton_pop {
    tag "singleton_pop_${pop}_${dataset}"
    label "medmem"

    input:
        tuple val(pop), val(dataset), file(pop_vcf), file(pop_sample)
    
    output:
        tuple val(pop), val(dataset), file(pop_vcf), file(pop_sample), file("${base}.singletons.bed")
    
    script:
        base = file(pop_vcf.baseName).baseName
        """
        vcftools --gzvcf ${pop_vcf} --singletons --out ${base}_tmp
        awk '{print \$1"\\t"\$2-1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5}' ${base}_tmp.singletons > ${base}.singletons.bed
        """
}

// Concatenate chromosome VCFs by population
process concat_pop {
    tag "concat_dataset_${pop}_${dataset}"
    label "medmem"
    
    input:
        tuple val(pop), val(datasets), val(vcfs), val(samples), val(chrms)
    
    output:
        tuple val(pop), val(dataset), file(pop_vcf), val(pop_sample)
        // tuple pop, file(pop_vcf)
    
    script:
        dataset = datasets[0]
        pop_sample = file(samples[0])
        pop_vcf = "${pop}_${dataset}.vcf.gz"
        """
        bcftools concat ${vcfs.join(' ')} -Oz -o ${pop}.tmp1.vcf.gz
        ## Recalculate AC, AN, AF
        bcftools +fill-tags ${pop}.tmp1.vcf.gz -Oz -o ${pop}.tmp2.vcf.gz
        bcftools sort ${pop}.tmp2.vcf.gz -T . -Oz -o ${pop_vcf}
        tabix ${pop_vcf}
        rm ${pop}.tmp*.vcf.gz
        """
}

def group_pops_data(it, groups) {
    group_pops_datas = []
    pop = it[0]
    pop_vcf = it[2]
    pop_sample = it[3]
    groups.each{ group, pops ->
        pops = pops.split(',')
        if( pop in pops){
            group_pops_datas << [ group, pop_vcf, pop_sample]
        }
    }
    return group_pops_datas
}

/*
    Step: Merge populations into group
*/
process merge_pop_groups {
    tag "merge_pop_groups_${group}_${chrm}"
    label "bigmem"
    
    input:
        tuple val(group), val(chrm), val(pop_vcfs), val(pop_samples)
    
    output:
        tuple val(group), val(chrm), file(vcf_out), file(sample_out)
    
    script:
        vcf_out = "${group}_${chrm}.vcf.gz"
        sample_out = "${group}_${chrm}.sample"
        """
        # Create local symlinks and index them to avoid remote directory access issues
        local_vcfs=""
        i=0
        for vcf in ${pop_vcfs.join(' ')}; do
            local_name="input_\${i}.vcf.gz"
            ln -sf "\${vcf}" "\${local_name}"
            tabix -f -p vcf "\${local_name}"
            local_vcfs="\${local_vcfs} \${local_name}"
            i=\$((i+1))
        done

        bcftools merge \
            \${local_vcfs} \
            -Oz -o ${group}.tmp1.vcf.gz
        ## Recalculate AC, AN, AF
        bcftools +fill-tags ${group}.tmp1.vcf.gz -Oz -o ${group}.tmp2.vcf.gz
        bcftools query -l ${group}.tmp2.vcf.gz | sort > sample.txt
        bcftools sort ${group}.tmp2.vcf.gz -T . |\
        bcftools view -S sample.txt -Oz -o ${vcf_out}
        cat ${pop_samples.join(' ')} > ${sample_out}
        rm ${group}.tmp*.vcf.gz
        """
}

process cat_chrm_groups {
    tag "cat_chrm_groups_${group}"
    label "extrabig"
    
    input:
        tuple val(group), val(group_vcfs), val(group_samples)
    
    output:
        tuple val(group), file(vcf_out), file(sample_out)
    
    script:
        vcf_out = "${group}.vcf.gz"
        sample_out = "${group}.sample"
        """
        bcftools concat ${group_vcfs.join(' ')} |\
        bcftools +fill-tags |\
        bcftools sort -T . -Oz -o ${vcf_out}
        tabix -f ${vcf_out}
        cat ${group_samples[0]} > ${sample_out}
        """
}

process combine_csv {
    tag "combine_csv_${dataset}"
    publishDir "${params.outDir}/REPORTS/fisher/PGX_ONLY/${dataset}/", mode:'copy'
    label "medmem"

    input:
        tuple val(dataset), val(fisher_csvs), val(chrms)

    output:
        tuple val(dataset), file(combined_fisher_csv)

    script:
        combined_fisher_csv = "${dataset}.all.assoc.fisher.csv"
        """
        head -n1 ${fisher_csvs[0]} > ${combined_fisher_csv}
        tail -q -n+2 ${fisher_csvs.join(' ')} >> ${combined_fisher_csv}
        """
}

process combine_csv_simple {
    tag "cat_${dataset}"
    label "bigmem"
    // publishDir "${params.outDir}/freq/${dataset}/", mode: 'copy'

    input:
        tuple val(dataset), val(csv_files)

    output:
        tuple val(dataset), file(csv_out)

    script:
        csv_out = "${dataset}_all.csv"
        """
        echo -e 'POP1\\tPOP1\\tWeighted FST' > ${csv_out}
        cat ${csv_files.join(' ')} >> ${site}
        """
}
