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

process get_sites{
  tag "sites_${dataset}_${chrm}"
  label "bigmem"

  input:
    tuple val(dataset), file(vcf_file), val(chrm)

  output:
    tuple val(dataset), file(vcf_file), val(chrm)
    tuple val(dataset), file(vcf_file), file(sites), val(chrm)

  script:
    sites = "${file(vcf_file.baseName).baseName}.sites.tsv"
    """
    tabix ${vcf_file}
    zcat ${vcf_file} | \
    SnpSift extractFields - -e "." ID CHROM POS REF ALT AC MAF AGVP_AF SAHGP_AF KG_AF KG_AFR_AF gnomAD_AF ExAC_AF ExAC_AFR_AF gnomAD_AF gnomAD_AFR_AF TOPMED CLNDISDB CLNDN CLNSIG "ANN[0].EFFECT" "ANN[0].GENE" >> ${sites}
    """
}

process get_gene_vcf {
    tag "get_gene_vcf_${dataset}_${chrm}_${gene_list_name}"
    label "bigmem"

    input:
        tuple val(dataset), file(dataset_vcf), file(dataset_sample), val(chrm), val(gene_list_name), file(gene_list_bed)

    output:
        tuple val(dataset), file(gene_list_vcf), file(dataset_sample), val(chrm)

    script:
        base = "${file(dataset_vcf.baseName).baseName}_${gene_list_name}"
        gene_list_vcf = "${base}.vcf.gz"
        """
        awk '{print \$1"\\t"\$2-10000"\\t"\$3+10000"\\t"\$4}' ${gene_list_bed} > ${gene_list_name}_tab.bed
        tabix ${dataset_vcf}
        bcftools view --regions-file ${gene_list_name}_tab.bed ${dataset_vcf} | \
        bcftools sort -T . -Oz -o ${gene_list_vcf}
        bcftools annotate --set-id "%CHROM\\:%POS" ${gene_list_vcf} | \
        bcftools query -f '%ID\\n' >> ${base}.sites
        rm -f ${dataset}_tmp.vcf.gz
        tabix ${gene_list_vcf}
        """
}

process get_vcf_ind {
  tag "vcf_ind_${sample_id}_${dataset}"
  label "small"

  input:
        tuple val(sample_id), file(dataset_vcf), val(dataset), val(chrm)
  output:
        tuple val(dataset), file(sample_vcf), val(sample_id), val(chrm)
  script:
      base = "${dataset}_${sample_id}_${chrm}"
      sample_vcf = "${base}.vcf.gz"
      """
      ## Generate vcf for single sample and only heterozygote variants
      bcftools view -s ${sample_id} ${dataset_vcf} | bcftools view -i '(GT="het" & GT!~"\\.") || (GT="AA" & GT!~"\\.")' | bgzip -c > ${sample_vcf}
      """
    //    ## Generate VCF of protein conding variants
    //   bcftools view -s ${sample} ${vcf_file_coding} | bcftools view -i '(GT="het" & GT!~"\\.") || (GT="AA" & GT!~"\\.")' | bgzip -c > ${base}_${sample}_coding.vcf.gz
    //   awk '/${sample}/ {print \$1"\t"\$2}' ${dataset_singletons} > ${sample}.singletons
}

process get_protein_coding_variants {
    tag "get_protein_coding_${dataset}_${chrm}"
    label "bigmem"

    input:
        tuple val(dataset), file(dataset_vcf), file(dataset_sample), val(chrm)

    output:
        tuple val(dataset), file(dataset_exome_vcf), file(dataset_sample), val(chrm)

    script:
        base = "${file(dataset_vcf.baseName).baseName}"
        dataset_exome_vcf = "${base}_exome.vcf.gz"
        """
        ## Generate VCF of protein conding variants
        zcat ${dataset_vcf} | \
        SnpSift filter " \
            (ANN[*].EFFECT has 'synonymous_variant') || \
            (ANN[*].EFFECT has 'missense_variant') || \
            (ANN[*].EFFECT has 'stop_gained') || \
            (ANN[*].EFFECT has 'frameshift_variant') || \
            (ANN[*].EFFECT has 'stop_lost') || \
            (ANN[*].EFFECT has 'inframe_insertion') || \
            (ANN[*].EFFECT has 'inframe_deletion') || \
            (ANN[*].EFFECT has 'coding_sequence_variant') " | \
        bgzip -c > ${dataset_exome_vcf}
        tabix ${dataset_exome_vcf}
        """
}

process get_gene_vcf1 {
    tag "get_gene_${dataset}_${gene_list_name}_${chrm}"
    label "bigmem"

    input:
        tuple val(dataset), file(dataset_vcf), file(dataset_sample), val(chrm), val(gene_list_name), file(gene_list_file)

    output:
        tuple val(dataset), file(gene_list_vcf), val(gene_list_name), val(chrm)

    script:
        base = "${file(dataset_vcf.baseName).baseName}_${gene_list_name}"
        gene_list_vcf = "${base}.vcf.gz"
        """
        tabix ${dataset_vcf}
        bcftools view --regions-file ${gene_list_file} ${dataset_vcf} | \
        bcftools sort -T . -Oz -o ${gene_list_vcf}
        """
}

process generate_chunks_vcf {
    tag "generate_chunks_${dataset}"
    label "bigmem"

    input:
        tuple val(dataset), file(vcf), file(vcf_idx), file(sample_file), file(mapFile), val(chrms), val(chunk_size)
    output:
        tuple val(dataset), file(vcf), file(vcf_idx), file(sample_file), file(chunkFile)
    script:
        chromosomes = chrms
        chunk = ''
        chunkFile = "chunks.txt"
        template "generate_chunks.py"
}

process split_vcf_to_chunk {
    tag "split_${dataset}_${chrm}:${start}-${end}"
    label "bigmem"

    input:
        tuple val(dataset), val(chrm), val(start), val(end), file(dataset_vcf), file(dataset_vcf_idx), file(dataset_sample)
    output:
        tuple val(dataset), file(vcf_chunk_out), file(dataset_sample), val(chrm)
    script:
        vcf_chunk_out = "${dataset_vcf.getSimpleName()}_${chrm}_${start}-${end}.bcf"
        """
        bcftools view --regions ${chrm}:${start}-${end} ${dataset_vcf} --threads ${task.cpus} -Ob -o ${vcf_chunk_out}
        tabix ${vcf_chunk_out}
        """
}


process get_gene_vcf_simple {
    tag "get_gene_${dataset_name}_${chrm}"
    label "bigmem"

    input:
        tuple val(dataset_name), file(dataset_vcf), file(gene_list_name), val(chrm)

    output:
        tuple val(dataset_name), file(gene_list_vcf), val(chrm)

    script:
        base = "${file(dataset_vcf.baseName).baseName}"
        gene_list_vcf = "${base}.vcf.gz"
        """
        tabix ${dataset_vcf}
        bcftools view --regions-file ${gene_list_name} ${dataset_vcf} | \
        bcftools sort -T . -Oz -o ${gene_list_vcf}
        """
}

process vcf_to_plink {
    tag "vcf_to_plink_${dataset}_${chrm}"
    label "bigmem"

    input:
        tuple val(dataset), file(dataset_vcf), file(dataset_sample), val(chrm)
    output:
        tuple val(dataset), file(plink_bed), file(plink_bim), file(plink_fam), file(dataset_sample), val(chrm)
    script:
        plink_bed = "${dataset}_${chrm}.bed"
        plink_bim = "${dataset}_${chrm}.bim"
        plink_fam = "${dataset}_${chrm}.fam"
        """
        bcftools annotate \
            -x INFO,^FORMAT \
            --set-id +'%CHROM\\_%POS\\_%REF\\_%ALT' ${dataset_vcf} -Oz -o ${dataset_vcf}_tmp.vcf.gz
        plink2 \
            --vcf ${dataset_vcf}_tmp.vcf.gz \
            --make-bed \
            --max-alleles 2 --vcf-half-call missing \
            --out ${dataset}_${chrm}
        """
}

process generate_pheno_fam_sample {
    tag "pheno_${dataset}_${chrm}"
    label "small"

    input:
        tuple val(dataset), file(plink_bed), file(plink_bim), file(plink_fam), file(dataset_sample), val(chrm)
    output:
        tuple val(dataset), file(plink_bed), file(plink_bim), file(plink_fam), file(pheno), val(chrm)
    script:
        inSample = dataset_sample
        inFam = plink_fam
        pheno = "${dataset}_${chrm}_pheno"
        template "generate_pheno_fam_sample.py"
}

process generate_pheno_fam_sample1 {
    tag "pheno_${dataset}_${chrm}"
    label "small"

    input:
        tuple val(dataset), file(plink_bed), file(plink_bim), file(plink_fam), file(dataset_sample), val(chrm)
    output:
        tuple val(dataset), file(plink_bed), file(plink_bim), file(plink_fam), file(pheno), val(chrm)
    script:
        inSample = dataset_sample
        inFam = plink_fam
        pheno = "${dataset}_${chrm}_pheno"
        template "generate_pheno_fam_sample.1.py"
}


process merge_groups {
    tag "merge_groups_${dataset}"
    label "largemem"

    input:
        tuple val(dataset), val(vcfs)
    output:
        tuple val(dataset), file(vcf_out)
    script:
        vcf_out = "${dataset}.bcf"
        if(vcfs.size() > 1){
            """
            bcftools merge ${vcfs.join(' ')} |\
            bcftools sort -T . -Ob -o ${vcf_out}
            """
        }
        else{
            """
            bcftools sort ${vcfs.join(' ')} -T . -Ob -o ${vcf_out}
            """
        }

}

process filter_alt_contigs {
    tag "filter_alt_${dataset}"
    label "medium"

    input:
        tuple val(dataset), file(dataset_vcf), file(vcf_idx), file(dataset_sample), file(dataset_map)

    output:
        tuple val(dataset), file(dataset_vcf), file(vcf_idx), file(dataset_sample), file(filtered_map)

    script:
        base = file(dataset_vcf.baseName).baseName
        filtered_map = "${base}_filtered.map"
        """
        # Filter out alternate contigs and non-standard chromosomes
        # Keep only chr1-22 and chrX (matching gnomAD coverage)
        grep -v '_' ${dataset_map} | grep -P '^chr(\\d+|X)\\t' > ${filtered_map} || touch ${filtered_map}
        """
}