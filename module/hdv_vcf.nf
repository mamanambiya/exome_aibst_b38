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

process hdv {
    tag "hdv_${dataset}"
    label "extrabig"
    maxForks 50

    input:
      tuple val(dataset), file(frq)

    output:
      tuple val(dataset), file(outFile) , file("${outFile}.sites.tsv"), file("${outFile}.all.tsv"), file("${outFile}.genes.tsv")

    script:
      outFile = "${frq.baseName}_hdv.tsv"
      annot = frq
      template "highly_diff_variants.py"
}

process hdv_dataset {
    tag "hdv_${dataset}_${chrm}_${fold}_ac${ac}_${group_treshold}"
    publishDir "${params.outDir}/hdv/pgx_only/${dataset}", mode:'copy', overwrite: true
    label "bigmem"
    maxForks 100

    input:
      tuple val(dataset), file(frq), file(pvalues), file(fst_file), val(chrm), val(fold), val(test_pops), val(base_pops), val(ac), val(group_treshold), val(fst_cutoff)

    output:
      tuple val(dataset), file("${outFile}_all.tsv"), file("${outFile}_base.tsv"), val(fold), val(ac), val(group_treshold)
      // tuple val(dataset), file("${outFile}_top.tsv") , file("${outFile}_all_sites.tsv"), file("${outFile}_all.tsv"), file("${outFile}_all_genes.tsv"), file("${outFile}_CLNDN.tsv")

    script:
      outFile = "${dataset}_${chrm}_hdv-dataset-${fold}_ac${ac}_${group_treshold}p"
      annot = frq
      template "highly_diff_variants_dataset.py"
}

process combine_hdvs {
    tag "combine_hdv_${dataset}_${type}"
    publishDir "${params.outDir}/hdv/pgx_only/${dataset}", mode:'copy', overwrite: true
    label "bigmem"

    input:
        tuple val(dataset), val(hdvs), val(fold), val(ac), val(group_treshold), val(type)

    output:
        tuple val(dataset), file(hdvs_out), val(fold), val(ac), val(group_treshold)

    script:
        hdvs_out = "${dataset}_hdv-${type}_all.csv"
        """
        # HDV All
        head -n1 ${hdvs[0]} > ${hdvs_out}
        tail -q -n+2 ${hdvs.join(' ')} >> ${hdvs_out}
        """
}

process annotate_hdvs {
    tag "annotate_hdv_${dataset}"
    publishDir "${params.outDir}/hdv/pgx_only/${dataset}", mode:'copy', overwrite: true
    label "bigmem"

    input:
        tuple val(dataset), file(hdv), file(annot)

    output:
        tuple val(dataset), file(hdv_out)

    script:
        hdv_out = "${hdv.getSimpleName()}_annotated.csv"
        template "annot_hdv.py"
}