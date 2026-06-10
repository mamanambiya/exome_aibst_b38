#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
  Novel variant detection - v2 (2026-03-19)
  Uses dbSNP b156 RS and FREQ fields + gnomAD v4 AF to determine novelty.
  
  A variant is "truly novel" if:
    1. Not in dbSNP b156 (no RS field)
    2. Not in 1000 Genomes (no 1000Genomes entry in FREQ)
    3. Not in gnomAD (no GnomAD entry in FREQ AND no AF_gnomAD_AFR from separate annotation)
    4. Not in TOPMED (no TOPMED entry in FREQ)
    5. Not in ClinVar (no CLNSIG field)
*/

process novel_in_existing_database_v2 {
    tag "novel_v2_${dataset}_${chrm}"
    label "bigmem"
    publishDir "${params.outDir}/novel/${dataset}", mode:'copy', pattern:'*.csv'

    input:
      tuple val(dataset), file(vcf_file), val(chrm)

    output:
      tuple val(dataset), file(novel_sites), file(counts), val(chrm)

    script:
      base = "${dataset}_${chrm}"
      novel_sites = "${base}_novel_sites.csv"
      counts = "${base}_not_in_other_dbs.csv"
      """
      #!/bin/bash
      set -euo pipefail

      tabix -f ${vcf_file}
      
      TOTAL=\$(bcftools view -H ${vcf_file} | wc -l)
      
      echo -e "Group;${dataset}" > ${counts}
      echo -e "Chromosome;${chrm}" >> ${counts}
      echo -e "Total polymorphic variants;\${TOTAL}" >> ${counts}
      
      # 1. Not in dbSNP b156 (RS field is missing)
      NOT_DBSNP=\$(bcftools query -f '%INFO/RS\n' ${vcf_file} | awk '\$1=="."' | wc -l)
      echo -e "Not in dbSNP b156;\${NOT_DBSNP}" >> ${counts}
      
      # 2-4. Parse FREQ field for 1000G, gnomAD, TOPMED presence
      # Also check RS for dbSNP and CLNSIG for ClinVar
      bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/RS\t%INFO/FREQ\t%INFO/CLNSIG\t%INFO/AF_gnomAD_AFR\n' ${vcf_file} 2>/dev/null | \
      awk -F'\t' '
      BEGIN {
        no_dbsnp=0; no_1kg=0; no_gnomad=0; no_topmed=0; no_clinvar=0; no_gnomad_afr=0; truly_novel=0; total=0
      }
      {
        total++
        in_dbsnp = (\$5 != ".")
        in_1kg = (\$6 ~ /1000Genomes/)
        in_gnomad = (\$6 ~ /GnomAD/)
        in_topmed = (\$6 ~ /TOPMED/)
        in_clinvar = (\$7 != ".")
        in_gnomad_afr = (\$8 != "." && \$8 != "0")
        
        if (!in_dbsnp) no_dbsnp++
        if (!in_1kg) no_1kg++
        if (!in_gnomad) no_gnomad++
        if (!in_topmed) no_topmed++
        if (!in_clinvar) no_clinvar++
        if (!in_gnomad_afr) no_gnomad_afr++
        
        # Truly novel: not in any database
        if (!in_dbsnp && !in_1kg && !in_gnomad && !in_topmed && !in_clinvar && !in_gnomad_afr) {
          truly_novel++
          print \$1"\t"\$2"\t"\$3"\t"\$4
        }
      }
      END {
        print "Not in 1000 Genomes;" no_1kg > "/dev/stderr"
        print "Not in gnomAD (dbSNP FREQ);" no_gnomad > "/dev/stderr"
        print "Not in TOPMED;" no_topmed > "/dev/stderr"
        print "Not in ClinVar;" no_clinvar > "/dev/stderr"
        print "Not in gnomAD v4 AFR;" no_gnomad_afr > "/dev/stderr"
        print "Truly novel (absent from all);" truly_novel > "/dev/stderr"
      }
      ' > ${novel_sites} 2>> ${counts}
      """
}
