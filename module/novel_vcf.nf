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

process novel_in_existing_database{
    tag "no_dbs_${dataset}_${chrm}"
    label "bigmem"
    publishDir "${params.outDir}/novel/${dataset}", mode:'copy', pattern:'*.csv'

    input:
      tuple val(dataset), file(vcf_file), val(chrm)

    output:
      tuple val(dataset), file(vcf_tmp), file(counts), file(novel_sites), val(chrm)

    script:
      base = "${dataset}_${chrm}"
      vcf_tmp = "${base}_ids.vcf.gz"
      novel_sites = "${base}_novel_sites.csv"
      counts = "${base}_not_in_other_dbs.csv"
      agvp_sites = params.agvp_b38_sites ?: ''
      """
      #!/bin/bash
      set -euo pipefail

      tabix -f ${vcf_file}

      TOTAL=\$(bcftools view -H ${vcf_file} | wc -l)

      echo -e "Group;${dataset}" > ${counts}
      echo -e "Chromosome;${chrm}" >> ${counts}
      echo -e "Total polymorphic variants;\${TOTAL}" >> ${counts}

      # 1. Not in dbSNP b156 (RS field missing = not in dbSNP, 1KG, TOPMED, ExAC, ESP)
      NOT_DBSNP=\$(bcftools query -f '%INFO/RS\n' ${vcf_file} | awk '\$1=="."' | wc -l)
      echo -e "Not in dbSNP b156;\${NOT_DBSNP}" >> ${counts}

      # 2. Check gnomAD v4.1 (AF_gnomAD_joint or AF field from gnomAD annotation)
      # After gnomAD annotation step, variants have gnomAD AF in INFO
      NOT_GNOMAD=\$(bcftools query -f '%INFO/AF_joint_afr\n' ${vcf_file} 2>/dev/null | awk '\$1=="." || \$1=="0"' | wc -l || echo "0")
      echo -e "Not in gnomAD v4.1;\${NOT_GNOMAD}" >> ${counts}

      # 3. Check AGVP (position lookup against lifted-over b38 sites)
      if [ -n "${agvp_sites}" ] && [ -f "${agvp_sites}" ]; then
        # Create position list from VCF
        bcftools query -f '%CHROM:%POS\n' ${vcf_file} > vcf_positions.txt
        # Count positions NOT in AGVP
        NOT_AGVP=\$(comm -23 <(sort vcf_positions.txt) <(sort ${agvp_sites}) | wc -l)
        echo -e "Not in AGVP;\${NOT_AGVP}" >> ${counts}
      else
        NOT_AGVP=\${TOTAL}
        echo -e "Not in AGVP;\${NOT_AGVP} (AGVP sites file not provided)" >> ${counts}
      fi

      # 4. ClinVar check
      NOT_CLINVAR=\$(bcftools query -f '%INFO/CLNSIG\n' ${vcf_file} 2>/dev/null | awk '\$1=="."' | wc -l || echo "0")
      echo -e "Not in ClinVar;\${NOT_CLINVAR}" >> ${counts}

      # Determine truly novel: not in ANY database
      # Extract all fields and check each variant
      bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/RS\t%INFO/CLNSIG\n' ${vcf_file} 2>/dev/null | \
      awk -F'\t' -v agvp_file="${agvp_sites}" '
      BEGIN {
        # Load AGVP sites if available
        if (agvp_file != "" && agvp_file != "null") {
          while ((getline line < agvp_file) > 0) {
            agvp[line] = 1
          }
          close(agvp_file)
          has_agvp = 1
        } else {
          has_agvp = 0
        }
        truly_novel = 0
        total = 0
      }
      {
        total++
        in_dbsnp = (\$5 != ".")
        in_clinvar = (\$6 != ".")

        pos_key = \$1":"\$2
        in_agvp = (has_agvp && (pos_key in agvp))

        # Novel = not in dbSNP AND not in AGVP
        # (gnomAD check done separately post-annotation; dbSNP covers 1KG/TOPMED/ExAC/gnomAD-via-FREQ)
        if (!in_dbsnp && !in_clinvar && !in_agvp) {
          truly_novel++
          print \$1"\t"\$2"\t"\$3"\t"\$4
        }
      }
      END {
        print "Truly novel (absent from dbSNP+ClinVar+AGVP);" truly_novel > "/dev/stderr"
      }
      ' > ${novel_sites} 2>> ${counts}

      # Create ID-annotated VCF for downstream
      bcftools annotate --set-id +"%CHROM\\_%POS\\_%REF\\_%FIRST_ALT" ${vcf_file} | \
        bcftools +fill-tags | \
        bgzip -c > ${vcf_tmp}
      tabix -f ${vcf_tmp}
      """
}


process combine_novel_count {
  errorStrategy "ignore"
  tag "combine_novel_count_${group}_${dataset}"
  label "bigmem"
  publishDir "${params.outDir}/novel/${dataset}/combined", mode:'copy', pattern:'*.csv'

  input:
    tuple val(dataset), val(counts), val(header), val(group)
  
  output:
    tuple val(dataset), file("${count}.total.csv"), file("${count}.summary.csv"), val(group)
  
  script:
    counts = counts.join(' ')
    count = "${dataset}_${group}_not_in_other_db"
    template "combine_novel_count.py"
}

process combine_counts {
  tag "combine_counts_${group}_${dataset}"
  label "bigmem"
  publishDir "${params.outDir}/counts/${dataset}/combined", mode:'copy', pattern:'*.csv'

  input:
    tuple val(dataset), val(counts), val(header), val(group)
  
  output:
    tuple val(dataset), file("${count}.total.csv"), path("${count}.summary.csv"), val(group)
  
  script:
    counts = counts.join(' ')
    count = "${dataset}.counts.${group}"
    template "combine_novel_count.py"
}
 

/*
  Hard novel variants. Checking sites in other databases not in dbSNP such as gnomad, exac, topmed
*/
process get_novel_sites {
  tag "get_novel_sites_${dataset}_${chrm}"
  label "bigmem"

  input:
    tuple val(dataset), file(vcf_file), file(vcf_sites_bed), val(chrm)
  
  output:
    tuple val(dataset), file(vcf_file), file(novel_sites_file), val(chrm)
  
  script:
    base = file(vcf_file.baseName).baseName
    novel_sites_file = "${base}_novel.csv"
    template "get_novel_from_tsv.py"
}

// Step: Get singleton for each population
process singleton_dataset_chrm {
    tag "singleton_pop_${dataset}_${chrm}"
    label "bigmem"

    input:
        tuple val(dataset), file(dataset_vcf), val(chrm)
    
    output:
        tuple val(dataset), file(dataset_vcf), file(sites), val(chrm), file(bed)
    
    script:
        base = file(dataset_vcf.baseName).baseName
        bed = "${base}.singletons.csv"
        sites = "${base}.singletons.sites.csv"
        """
        vcftools --gzvcf ${dataset_vcf} --singletons --out ${base}_tmp
        awk '{print \$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5}' ${base}_tmp.singletons > ${bed}
        awk '{print \$1"\\t"\$2}' ${base}_tmp.singletons | tail -q -n+2 > ${sites}
        """
}

process count_variants{
    tag "count_variants_${dataset}_${chrm}"
    label "bigmem"
    publishDir "${params.outDir}/counts/", mode:'copy', pattern:'*.csv'

    input:
      tuple val(dataset), file(vcf_file), val(sample_id), file(singletons), val(chrm)

    output:
      tuple val(dataset), file(count_variants), val(sample_id), val(chrm)

    script:
      base = "${dataset}_${sample_id}_${chrm}"
      count_variants = "${base}_counts.csv"
      """
      tabix ${vcf_file}
      echo -e "Group ; ${dataset}_${chrm}" > ${count_variants}
      ## Get total number positions of variants
      echo -e "Total variants ; \$(zcat ${vcf_file} | grep -v "^#" | wc -l)" >> ${count_variants}
      ## Total SNV variants
      echo -e "SNV positions ; \$(bcftools view -v snps ${vcf_file} | grep -v "^#" | wc -l )" >> ${count_variants}
      ## Total INDELS variants
      echo -e "INDELS positions ; \$(bcftools view -v indels ${vcf_file} | grep -v "^#" | wc -l )" >> ${count_variants}
      ## Total biallelic variants
      echo -e "Biallelic positions ; \$(bcftools view -m2 -M2 ${vcf_file} | grep -v "^#" | wc -l )" >> ${count_variants}
      ## Total biallelic SNV variants
      echo -e "Biallelic SNV positions ; \$(bcftools view -m2 -M2 -v snps ${vcf_file} | grep -v "^#" | wc -l )" >> ${count_variants}
      ## Total biallelic INDELS variants
      echo -e "Biallelic INDELS positions ; \$(bcftools view -m2 -M2 -v indels ${vcf_file} | grep -v "^#" | wc -l )" >> ${count_variants}
      ## Total multiallelic variants
      echo -e "Multiallelic positions ; \$(bcftools view -m3 ${vcf_file} | grep -v "^#" | wc -l )" >> ${count_variants}
      ## Total multiallelic SNV variants
      echo -e "Multiallelic SNV positions ; \$(bcftools view -m3 -v snps ${vcf_file} | grep -v "^#" | wc -l )" >> ${count_variants}
      ## Total multiallelic INDELS variants
      echo -e "Multiallelic INDELS positions ; \$(bcftools view -m3 -v indels ${vcf_file} | grep -v "^#" | wc -l )" >> ${count_variants}
      ## Total Singleton positions
      echo -e "Singletons ; \$(bcftools view ${vcf_file} --regions-file ${singletons} | grep -v "^#" | wc -l )" >> ${count_variants}
      ## Rare variants
      echo -e "Rare variants ; \$(bcftools view -i 'MAF<=0.01' ${vcf_file} | grep -v "^#" | wc -l)" >> ${count_variants}
      ## Nonsynonymous/missense
      echo -e "Nonsynonymous SNV ; \$(zcat ${vcf_file} | SnpSift filter "(ANN[*].EFFECT has 'missense_variant')" | grep -v "^#" | wc -l)" >> ${count_variants}
      ## Synonymous
      echo -e "Synonymous SNV ; \$(zcat ${vcf_file} | SnpSift filter "(ANN[*].EFFECT has 'synonymous_variant')" | grep -v "^#" | wc -l)" >> ${count_variants}
      ## Stop gained
      echo -e "Stop gained SNV ; \$(zcat ${vcf_file} | SnpSift filter "(ANN[*].EFFECT has 'stop_gained')" | grep -v "^#" | wc -l)" >> ${count_variants}
      ## Stop lost
      echo -e "Stop lost SNV ; \$(zcat ${vcf_file} | SnpSift filter "(ANN[*].EFFECT has 'stop_lost')" | grep -v "^#" | wc -l)" >> ${count_variants}
      ## Splicing
      echo -e "Splicing SNV ; \$(zcat ${vcf_file} | SnpSift filter "(ANN[*].EFFECT has 'splice_site_region')" | grep -v "^#" | wc -l)" >> ${count_variants}
      ## LOF
      echo -e "LOF SNV ; \$(zcat ${vcf_file} | SnpSift filter "(LOF[*].PERC >= 0.5)" | grep -v "^#" | wc -l)" >> ${count_variants}
      """
}

process novel_count{
    tag "novel_count_${dataset}_${chrm}"
    label "bigmem"
    publishDir "${params.outDir}/novel/", mode:'copy', pattern:'*.csv'

    input:
      tuple val(dataset), file(vcf_file), val(chrm)

    output:
      tuple val(dataset), file(count_novel), val(chrm)

    script:
    base = file(vcf_file.baseName).baseName
      count_novel = "${base}_count.csv"
      """
      tabix ${vcf_file}
      echo -e "Group ; ${dataset}_${chrm}" > ${count_variants}
      ## Get total number positions of variants
      echo -e "Novel variants ; \$(zcat ${vcf_file} | grep -v "^#" | wc -l)" >> ${count_novel}
      ## Novel rare variants
      echo -e "Novel Rare variants ; \$(bcftools view -i 'MAF<=0.01' ${vcf_file} | grep -v "^#" | wc -l)" >> ${count_novel}
      ## Nonsynonymous/missense
      echo -e "Novel Nonsynonymous SNV ; \$(zcat ${vcf_file} | SnpSift filter "(ANN[*].EFFECT has 'missense_variant')" | grep -v "^#" | wc -l)" >> ${count_novel}
      ## Synonymous
      echo -e "Novel Synonymous SNV ; \$(zcat ${vcf_file} | SnpSift filter "(ANN[*].EFFECT has 'synonymous_variant')" | grep -v "^#" | wc -l)" >> ${count_novel}
      ## Stop gained
      echo -e "Novel Stop gained SNV ; \$(zcat ${vcf_file} | SnpSift filter "(ANN[*].EFFECT has 'stop_gained')" | grep -v "^#" | wc -l)" >> ${count_novel}
      ## Stop lost
      echo -e "Novel Stop lost SNV ; \$(zcat ${vcf_file} | SnpSift filter "(ANN[*].EFFECT has 'stop_lost')" | grep -v "^#" | wc -l)" >> ${count_novel}
      ## Splicing
      echo -e "Novel Splicing SNV ; \$(zcat ${vcf_file} | SnpSift filter "(ANN[*].EFFECT has 'splice_site_region')" | grep -v "^#" | wc -l)" >> ${count_novel}
      ## LOF
      echo -e "Novel LOF SNV ; \$(zcat ${vcf_file} | SnpSift filter "(LOF[*].PERC >= 0.5)" | grep -v "^#" | wc -l)" >> ${count_novel}
      """
}


process novel_count_combine{
    tag "novel_count_combine_${dataset}_${suffix}"
    label "small"
    publishDir "${params.outDir}/novel/${dataset}/combined/", mode:'copy', pattern:'*.csv'

    input:
      tuple val(dataset), val(pops), val(count_csvs), val(chrms), val(suffix)

    output:
      tuple val(dataset), file(csv_out)

    script:
      csv_counts = count_csvs.join(',')
      csv_labels = pops.join(',')
      csv_chrms = chrms.join(',')
      csv_out = "${dataset}_${suffix}_combine_counts.csv"
      template "combine_counts.py"
}

process generate_circos_dataframe{
    tag "generate_circos_dataframe_${dataset}_${prefix}"
    label "bigmem"
    publishDir "${params.outDir}/REPORTS/hdv/PGX_ONLY/${dataset}/", mode:'copy'

    input:
      tuple val(dataset), val(sites), val(prefix)

    output:
      tuple val(dataset), file(csv_out)

    script:
      sites = sites.join(';')
      csv_out = "${dataset}_${prefix}_circos.csv"
      template "make_circos_data.py"
}

process generate_circos_plot{
    tag "generate_circos_plot_${dataset}"
    label "rplot"
    publishDir "${params.outDir}/REPORTS/hdv/PGX_ONLY/${dataset}/", mode:'copy'

    input:
      tuple val(dataset), val(circos_csv)

    output:
      tuple val(dataset), file(circos_plot)

    script:
      circos_plot = "${dataset}.circos.novel.pdf"
      template "novel_circos.R"
}


