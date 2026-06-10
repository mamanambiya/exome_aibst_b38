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

process fst_analysis {
   tag "fst_${pops}_${dataset}"
   label "small"
  //  publishDir "${params.work_dir}/REPORTS/${dataset}/PGX_ONLY/FST", mode:'copy', pattern:'*.tiff'
   input:
       tuple val(pops), val(dataset), val(pops_vcf_list), val(pops_sample_list)
   output:
       tuple val(dataset), val("${pop1}__${pop2}"), file("${fst_basename}.weir.fst"), file(fst_log)
   script:
      pop1 = pops.split('_')[0]
      pop2 = pops.split('_')[1]
      pop1_vcf = file(pops_vcf_list[0])
      pop2_vcf = file(pops_vcf_list[1])
      fst_basename = "${dataset}.${pop1}__${pop2}"
      fst_log = "${dataset}.${pop1}__${pop2}.fst.log"
      """
      tabix -f ${pop1_vcf}
      tabix -f ${pop2_vcf}
      bcftools norm -m -any ${pop1_vcf} |\
      bcftools view --min-ac 3 -Ob -o ${pop1}.temp.bcf
      tabix -f ${pop1}.temp.bcf
      bcftools norm -m -any ${pop2_vcf} |\
      bcftools view --min-ac 2 -Ob -o ${pop2}.temp.bcf
      tabix -f ${pop2}.temp.bcf
      bcftools merge \
          ${pop1}.temp.bcf ${pop2}.temp.bcf |\
      bcftools +fill-tags -Oz -o ${fst_basename}.fst.temp.vcf.gz
      ## Compute Fst for pair of populations, here the log file will be used to extract the weighted mean fst
      awk '{print \$1}' ${pops_sample_list[0]} > ${pops_sample_list[0]}.updated
      awk '{print \$1}' ${pops_sample_list[1]} > ${pops_sample_list[1]}.updated
      vcftools \
          --gzvcf ${fst_basename}.fst.temp.vcf.gz \
          --weir-fst-pop ${pops_sample_list[0]}.updated \
          --weir-fst-pop ${pops_sample_list[1]}.updated \
          --out ${fst_basename} \
          2> ${fst_log}
      rm -f *.temp.*cf
      """
}


// SnpSift filter "ANN[*].EFFECT has 'synonymous_variant'" |\


// process combine_fst_weir_estimates {
//    tag "comb_fst_${dataset}_${pops}"
//    label "small"
// //    publishDir "${params.work_dir}/data/PGX_ONLY/FST/${dataset}", mode:'symlink'
//    input:
//        tuple val(dataset), val(pops), file(fst_weir), file(fst_log)
//    output:
//        tuple val(dataset), val(pops), file(fst_out)
//    script:
//        pops_ = pops.split('_')
//        fst_out = "${dataset}_${pops}.weighted.fst.estimates.txt"
//        """
//         # Get mean Fst estimates for each of the comparisons
//         weir_fst=`grep "Weir and Cockerham weighted Fst estimate" ${fst_log} | awk -F":" '{print \$NF}'`
//         echo \"${pops_[0]} ${pops_[1]} \$weir_fst\" >> ${fst_out}
//        """
// }

// process get_fst_weir_estimates {
//    tag "weir_fst_${dataset}_${pops}"
//    label "small"
// //    publishDir "${params.work_dir}/data/PGX_ONLY/FST/${dataset}", mode:'symlink'
//    input:
//        tuple val(dataset), val(pops), file(fst_weir), file(fst_log)
//    output:
//        tuple val(dataset), val(pops), file(fst_weir_out), file(mean_fst_out)
//    script:
//        pops_ = pops.split('__')
//        mean_fst_out = "${dataset}_${pops}.mean.weighted.fst.estimates.txt"
//        fst_weir_out = "${dataset}_${pops}.weighted.fst.estimates.txt"
//        """
//         # Get mean Fst estimates for each of the comparisons
//         weir_fst=`grep "Weir and Cockerham weighted Fst estimate" ${fst_log} | awk -F":" '{print \$NF}'`
//         echo \"${pops_[0]} ${pops_[1]} \$weir_fst\" >> ${mean_fst_out}

//         # Get mean Fst estimates for each of the for each SNP
//         echo -e 'POPS\\tCHROM\\tPOS\\tFST' > ${fst_weir_out}
//         awk 'NR>1{print "${pops_[0]}_${pops_[1]}\\t"\$0}' ${fst_weir} >> ${fst_weir_out}
//        """
// }

process get_fst_weir_estimates {
    tag "weir_fst_${dataset}_${pops}"
    label "bigmem"
    publishDir "${params.outDir}/fst/pgx_only/${dataset}", mode:'copy', overwrite: true
    input:
       tuple val(dataset), val(pops), file(fst_weir), file(fst_log)
    output:
       tuple val(dataset), val(pops), file(fst_weir_out), file(mean_fst_out)
    script:
       mean_fst_out = "${dataset}_${pops}.mean.weighted.fst.estimates.txt"
       fst_weir_out = "${dataset}_${pops}.weighted.fst.estimates.txt"
       template "fst_process.py"
}

process get_fst_weir_estimates_cutoff {
    tag "weir_fst_${dataset}_${pops}"
    label "bigmem"
    publishDir "${params.outDir}/fst/pgx_only/${dataset}", mode:'copy', overwrite: true
    input:
       tuple val(dataset), val(pops), file(fst_weir), file(fst_log), val(fst_cutoff)
    output:
       tuple val(dataset), val(pops), file(fst_weir_out), file(mean_fst_out)
    script:
       mean_fst_out = "${dataset}_${pops}.mean.weighted.fst.estimates.txt"
       fst_weir_out = "${dataset}_${pops}.weighted.fst.estimates_${fst_cutoff.replace('.','_')}.txt"
       template "fst_process_cutoff.py"
}

process combine_fst_weir_estimates {
    tag "cat_fst_${dataset}"
    label "bigmem"
    publishDir "${params.outDir}/fst/pgx_only/${dataset}", mode:'copy', overwrite: true

    input:
        tuple val(dataset), val(csv_files)

    output:
        tuple val(dataset), file(csv_out)

    script:
        csv_out = "${dataset}_2pops_weir_fst_all.csv"
        """
        head -n1 ${csv_files[0]} > ${csv_out}
        tail -q -n+2 ${csv_files.join(' ')} >> ${csv_out}
        """
}

process combine_fst_mean_weir_estimates {
    tag "cat_${dataset}"
    label "bigmem"
    publishDir "${params.outDir}/fst/pgx_only/${dataset}", mode:'copy', overwrite: true
    // publishDir "${params.outDir}/freq/${dataset}/", mode: 'copy'

    input:
        tuple val(dataset), val(csv_files)

    output:
        tuple val(dataset), file(csv_out)

    script:
        csv_out = "${dataset}_2pops_mean_weir_fst_all.csv"
        """
        echo -e 'POP1\\tPOP1\\tWeighted_FST' > ${csv_out}
        cat ${csv_files.join(' ')} >> ${csv_out}
        """
}

process combine_weir_fst_analysis {
   tag "comb_weir_fst_${dataset}"
   label "bigmem"
//    publishDir "${params.outDir}/REPORTS/fst/PGX_ONLY/${dataset}/", mode:'copy'
    publishDir "${params.outDir}/fst/pgx_only/${dataset}", mode:'copy', overwrite: true

   input:
       tuple val(dataset), val(weir_fst_files), val(fst_cutoff)
   output:
       tuple val(dataset), file("${weir_fst_out}.csv")
   script:
       weir_fst_out = "${dataset}.weighted_weir-fst_estimates_${fst_cutoff.replace(".", '')}"
       template "combine_weir_fst.py"
}

/// Step 13.4: Generate FST matrix
process generate_fst_matrix{
   tag "generate_fst_matrix_${dataset}"
   label "bigmem"
   publishDir "${params.outDir}/fst/pgx_only/${dataset}", mode:'copy', overwrite: true
   input:
       tuple val(dataset), file(weighted_fst_estimates)
   output:
       tuple val(dataset), file(fst_matrix_output)
   script:
        fts_2by2_input = weighted_fst_estimates
        fst_matrix_output = "${dataset}_weighted_fst_estimates.matrix.tsv"
        template "convert_fst_result_to_matrix.py"
}

process plot_fst_matrix{
   tag "plot_fst_matrix_${dataset}"
   label "bigmem"
   publishDir "${params.outDir}/fst/pgx_only/${dataset}", mode:'copy', overwrite: true
   input:
       tuple val(dataset), file(fst_matrix)
   output:
       tuple val(dataset), file(fst_matrix_plot)
   script:
       fst_matrix_plot = "${dataset}_weighted_fst_estimates_png.pdf"
       template "plot_fst_matrix.py"
}


// process plot_fst_matrix{
//    tag "plot_fst_matrix_${dataset}"
//    label "rplot"
//    publishDir "${params.outDir}/fst/pgx_only/${dataset}", mode:'link', overwrite: true
//    input:
//        tuple val(dataset), file(weighted_fst_estimates)
//    output:
//        tuple val(dataset), file(weighted_fst_estimates_png)
//    script:
//        weighted_fst_estimates_png = "${dataset}_weighted_fst_estimates_png.tiff"
//        template "step13_4_plot_fst_matrix.R"
// }

// process combine_weir_fst_analysis {
//    tag "comb_weir_fst_${dataset}"
//    label "extrabig"
//    publishDir "${params.outDir}/REPORTS/fst/PGX_ONLY/${dataset}/", mode:'copy'

//    input:
//        tuple val(dataset), val(weir_fst_files)
//    output:
//        tuple val(dataset), file("${weir_fst_out}.csv"), file("${weir_fst_out}_HighDiff.csv"), file("${weir_fst_out}_HighDiff_pos.csv")//, file("${fst_out}_HighDiff.vcf.gz"), file("${fst_out}_HighDiff.ann.csv") 
//    script:
//        weir_fst_out = "${dataset}.weighted_weir-fst_estimates"
//        template "combine_weir_fst.py"
// }


// # Select columns where at least one of the Fst comparisons is highly differentiated (i.e. >0.5)

//  vcftools \
//            --gzvcf ${dataset_vcf} \
//            --positions ${"${fst_out}_HighDiff_pos.csv"} \
//            --recode-INFO-all --recode --stdout \
//            | gzip -c > ${"${fst_out}_HighDiff.vcf.gz"}
//        zcat ${fst_out}_HighDiff.vcf.gz \
//            | SnpSift extractFields - -e "." CHROM POS ID REF ALT AC AF MAF "ANN[0].GENE" AGVP_AF KG_AFR_AF ExAC_AFR_AF gnomAD_AFR_AF "ANN[0].EFFECT" CLNDN CDS GWASCAT_TRAIT GWASCAT_P_VALUE GWASCAT_PUBMED_ID\
//            | awk '{print \$1":"\$2"\\t"\$0}' > ${fst_out}_HighDiff.ann.csv