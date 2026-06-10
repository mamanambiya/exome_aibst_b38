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

process get_freq_1 {
    tag "frq_${dataset}_${chrm}"
    publishDir "${params.outDir}/freq/${dataset}/", mode:'copy', pattern:'*.frq'
    label "bigmem"

    input:
        tuple val(dataset), file(vcf_file), file(sample_file), val(chrm)

    output:
        tuple val(dataset), file(vcf_file), file(sample_file), file(vcf_sites), file(frq), file(frq_1), val(chrm)

    script:
        base = file(vcf_file.baseName).baseName
        frq = "${base}.frq"
        frq_1 = "${base}.mafs"
        vcf_sites = "${base}.sites"
        """
        tabix ${vcf_file}
        echo -e '#CHRM:POS\\tCHROM\\tPOS\\tGENE\\tEFFECT\\t${dataset}_AC_\\t${dataset}_MAF_\\tgnomAD_AF\\tKG_AF\\tAF_TGP\\tAF_EXAC\\tAF_ESP\\tTOPMED_AF\\tKNL_AF\\tKNK_AF\\tKNM_AF\\tKNP_AF\\tNGI_AF\\tNGH_AF\\tNGY_AF\\tTZA_AF\\tTZB_AF\\tSAV_AF\\tZWD_AF\\tZWS_AF\\tACB_AF\\tASW_AF\\tCEU_AF\\tESN_AF\\tGBR_AF\\tGWD_AF\\tLWK_AF\\tMSL_AF\\tYRI_AF\\tZulu_AF\\tAmhara_AF\\tBaganda_AF\\tGumuz_AF\\tOromo_AF\\tSomali_AF\\tWolayta_AF\\tCLNDISDB\\tCLNDN\\tCLNSIG' > ${frq}
        bcftools view \
            -m2 -M2 -v snps \
            ${vcf_file} | \
        bcftools +fill-tags | \
        bcftools annotate \
            --set-id +"%CHROM:%POS" | \
        SnpSift extractFields - -e "0.0" \
            ID CHROM POS ANN[0].GENE ANN[0].EFFECT AC MAF AGVP_AF KG_AF KG_AFR_AF gnomAD_AF gnomAD_AFR_AF TOPMED KNL_AF KNK_AF KNM_AF KNP_AF NGI_AF NGH_AF NGY_AF TZA_AF TZB_AF SAV_AF ZWD_AF ZWS_AF ACB_AF ASW_AF CEU_AF ESN_AF GBR_AF GWD_AF LWK_AF MSL_AF YRI_AF Zulu_AF Amhara_AF Baganda_AF Gumuz_AF Oromo_AF Somali_AF Wolaya_AF CLNDISDB CLNDN CLNSIG | \
            tail -n+2 >> ${frq}
        echo -e '#CHRM:POS\\t${dataset}_AF' > ${frq_1}
        awk -F'\\t' '{print \$2":"\$3"\\t"\$7}' ${frq} | tail -n+2 >> ${frq_1}
        echo -e '#CHRM:POS' > ${vcf_sites}
        bcftools query \
            -f "%CHROM:%POS\\n" \
            ${vcf_file} >> ${vcf_sites}
            
        """
}

process split_id_freq {
    tag "split_id_freq_${dataset}_${chrm}"
    publishDir "${params.outDir}/freq/${dataset}/", mode:'copy', pattern:'*.frq.tsv'
    label "bigmem"

    input:
        tuple val(dataset), file(vcf_file), file(sample_file), file(vcf_sites), file(frq), file(frq_1), val(chrm), file(pvalues_file)

    output:
        tuple val(dataset), file(out), file(pvalues_file), val(chrm)

    script:
        out = "${frq}.csv"
        template "split_rsid.py"
}

process split_id_freq_2 {
    tag "split_id_freq_${dataset}_${chrm}"
    publishDir "${params.outDir}/freq/${dataset}/", mode:'copy', pattern:'*.frq.tsv'
    label "bigmem"

    input:
        tuple val(dataset), file(frq_file), val(chrm), file(pvalue_file), file(fst_file)

    output:
        tuple val(dataset), file(out), val(chrm), file(pvalue_file), file(fst_file)

    script:
        out = "${frq_file}.csv"
        template "split_rsid.py"
}

process split_id_freq_1 {
    tag "split_id_freq_${dataset}_${chrm}"
    publishDir "${params.outDir}/freq/${dataset}/", mode:'copy', pattern:'*.csv'
    label "bigmem"

    input:
        tuple val(dataset), file(vcf_file), file(sample_file), val(chrm), file(frq_file)

    output:
        tuple val(dataset), file(vcf_file), file(sample_file), val(chrm), file(out)

    script:
        out = "${frq_file}.csv"
        template "split_rsid.py"
}

process get_freq {
    tag "frq_${dataset}_${chrm}"
    publishDir "${params.outDir}/freq/${dataset}/", mode:'copy', pattern:'*.frq'
    label "bigmem"

    input:
        tuple val(dataset), file(vcf_file), file(sample_file), val(chrm)

    output:
        tuple val(dataset), file(vcf_file), file(sample_file), file(vcf_sites), file(frq), file(frq_1), val(chrm)

    script:
        base = file(vcf_file.baseName).baseName
        frq = "${base}.frq"
        frq_1 = "${base}.mafs"
        vcf_sites = "${base}.sites"
        """
        # v2: added gnomAD_AFR_AF column
        tabix -f ${vcf_file}
        echo -e '#CHRM:POS\\tCHROM\\tPOS\\tREF\\tALT\\tGENE\\tEFFECT\\t${dataset}_AC_\\t${dataset}_MAF_\\tAGVP_AF\\tKG_AF\\tKG_AFR_AF\\tgnomAD_AF\\tgnomAD_AFR_AF\\tTOPMED_AF\\tKNL_AF\\tKNK_AF\\tKNM_AF\\tKNP_AF\\tNGI_AF\\tNGH_AF\\tNGY_AF\\tTZA_AF\\tTZB_AF\\tSAV_AF\\tZWD_AF\\tZWS_AF\\tACB_AF\\tASW_AF\\tCEU_AF\\tESN_AF\\tGBR_AF\\tGWD_AF\\tLWK_AF\\tMSL_AF\\tYRI_AF\\tZulu_AF\\tAmhara_AF\\tBaganda_AF\\tGumuz_AF\\tOromo_AF\\tSomali_AF\\tWolayta_AF\\tCLNDISDB\\tCLNDN\\tCLNSIG' > ${frq}
        bcftools view \
            -m2 -M2 -v snps \
            ${vcf_file} | \
        bcftools +fill-tags | \
        bcftools annotate \
            --set-id +"%CHROM:%POS" | \
        SnpSift extractFields - -e "0.0" \
            ID CHROM POS REF ALT ANN[0].GENE ANN[0].EFFECT AC MAF AGVP_AF KG_AF KG_AFR_AF gnomAD_AF gnomAD_AFR_AF TOPMED KNL_AF KNK_AF KNM_AF KNP_AF NGI_AF NGH_AF NGY_AF TZA_AF TZB_AF SAV_AF ZWD_AF ZWS_AF ACB_AF ASW_AF CEU_AF ESN_AF GBR_AF GWD_AF LWK_AF MSL_AF YRI_AF Zulu_AF Amhara_AF Baganda_AF Gumuz_AF Oromo_AF Somali_AF Wolayta_AF CLNDISDB CLNDN CLNSIG | \
            tail -n+2 >> ${frq}
        echo -e '#CHRM:POS\\t${dataset}_AF' > ${frq_1}
        awk -F'\\t' '{print \$2":"\$3"\\t"\$7}' ${frq} | tail -n+2 >> ${frq_1}
        echo -e '#CHRM:POS' > ${vcf_sites}
        bcftools query \
            -f "%CHROM:%POS\\n" \
            ${vcf_file} >> ${vcf_sites}  
        """
}

process get_freq2 {
    tag "frq_${dataset}_${chrm}"
    label "bigmem"

    input:
        tuple val(dataset), file(vcf_file), file(sample_file), val(chrm)

    output:
        tuple val(dataset), file(vcf_file), file(sample_file), file(vcf_sites), file(frq), file(frq_1), val(chrm)

    script:
        base = file(vcf_file.baseName).baseName
        frq = "${dataset}_${chrm}.frq"
        frq_1 = "${dataset}_${chrm}.mafs"
        vcf_sites = "${dataset}_${chrm}.sites"
        """
        tabix ${vcf_file}
        echo -e '#CHRM:POS\\tCHROM\\tPOS\\tGENE\\tEFFECT\\tAC\\tMAF\\tAGVP_AF\\tSAHGP_AF\\tKG_AF\\tKG_AFR_AF\\tExAC_AF\\tExAC_AFR_AF\\tgnomAD_AF\\tgnomAD_AFR_AF\\tTOPMED\\tCLNDISDB\\tCLNDN\\tCLNSIG' > ${frq}
        bcftools view \
            -m2 -M2 -v snps \
            ${vcf_file} | \
        bcftools +fill-tags | \
        bcftools annotate \
            --set-id "%CHROM:%POS" | \
        SnpSift extractFields - -e "0.0" \
            ID CHROM POS ANN[0].GENE ANN[0].EFFECT AC MAF AGVP_AF SAHGP_AF KG_AF KG_AFR_AF ExAC_AF ExAC_AFR_AF gnomAD_AF gnomAD_AFR_AF TOPMED CLNDISDB CLNDN CLNSIG | \
            tail -n+2 >> ${frq}
        echo -e '#CHRM:POS\\t${dataset}_AF' > ${frq_1}
        awk -F'\\t' '{print \$2":"\$3"\\t"\$7}' ${frq} | tail -n+2 >> ${frq_1}
        echo -e '#CHRM:POS' > ${vcf_sites}
        bcftools query \
            -f "%CHROM:%POS\\n" \
            ${vcf_file} >> ${vcf_sites}
            
        """
}

process get_freq_3 {
    tag "frq_${dataset}_${chrm}"
    publishDir "${params.outDir}/freq/${dataset}/", mode:'copy', pattern:'*.frq'
    label "bigmem"

    input:
        tuple val(dataset), file(vcf_file), file(sample_file), val(chrm), file(pvalue_file)

    output:
        tuple val(dataset), file(vcf_file), file(sample_file), file(vcf_sites), file(frq), file(frq_1), val(chrm), file(pvalue_file)

    script:
        base = file(vcf_file.baseName).baseName
        frq = "${base}.frq"
        frq_1 = "${base}.mafs"
        vcf_sites = "${base}.sites"
        """
        tabix -f ${vcf_file}
        echo -e '#CHRM:POS\\tCHROM\\tPOS\\tREF\\tALT\\tGENE\\tEFFECT\\t${dataset}_AC_\\t${dataset}_MAF_\\t${dataset}_AF_\\tAGVP_AF\\tKG_AF\\tKG_AFR_AF\\tgnomAD_AF\\tgnomAD_AFR_AF\\tTOPMED_AF\\tKNL_AF\\tKNK_AF\\tKNM_AF\\tKNP_AF\\tNGI_AF\\tNGH_AF\\tNGY_AF\\tTZA_AF\\tTZB_AF\\tSAV_AF\\tZWD_AF\\tZWS_AF\\tACB_AF\\tASW_AF\\tCEU_AF\\tESN_AF\\tGBR_AF\\tGWD_AF\\tLWK_AF\\tMSL_AF\\tYRI_AF\\tCLNDISDB\\tCLNDN\\tCLNSIG\\tCLNDISDB\\tCLNDN\\tCLNSIG' > ${frq}
        bcftools view \
            -m2 -M2 -v snps \
            ${vcf_file} | \
        bcftools +fill-tags | \
        bcftools annotate \
            --set-id +"%CHROM:%POS" | \
        SnpSift extractFields - -e "0.0" \
            ID CHROM POS REF ALT ANN[0].GENE ANN[0].EFFECT AC MAF AF AGVD_AGVD_AF 1KG_1KG_AF 1KG_afr_AFRAF gnomad_b37_gnomad_b37_AF gnomad_b37_afr_AF topmed_topmed_TOPMED AIBST_KNL_AF AIBST_KNK_AF AIBST_KNM_AF AIBST_KNP_AF AIBST_NGI_AF AIBST_NGH_AF AIBST_NGY_AF AIBST_TZA_AF AIBST_TZB_AF AIBST_SAV_AF AIBST_ZWD_AF AIBST_ZWS_AF 1KG_ACB_AF 1KG_ASW_AF 1KG_CEU_AF 1KG_ESN_AF 1KG_GBR_AF 1KG_GWD_AF 1KG_LWK_AF 1KG_MSL_AF 1KG_YRI_AF AGVD_Amhara_AF AGVD_Baganda_AF AGVD_Gumuz_AF AGVD_Oromo_AF AGVD_Somali_AF AGVD_Wolayta_AF AGVD_Zulu_AF CLNDISDB CLNDN CLNSIG | \
            tail -n+2 >> ${frq}
        echo -e '#CHRM:POS\\t${dataset}_AF' > ${frq_1}
        awk -F'\\t' '{print \$2":"\$3"\\t"\$7}' ${frq} | tail -n+2 >> ${frq_1}
        echo -e '#CHRM:POS' > ${vcf_sites}
        bcftools query \
            -f "%CHROM:%POS\\n" \
            ${vcf_file} >> ${vcf_sites}  
        """
}

process get_freq_3_3 {
    tag "frq_${dataset}_${chrm}"
    publishDir "${params.outDir}/freq/${dataset}/", mode:'copy', pattern:'*.frq'
    label "bigmem"

    input:
        tuple val(dataset), file(vcf_file), val(group), val(chrm), file(pvalue_file), file(fst_file), file(pop_mafs_vcf), file(pop_mafs_tbi)

    output:
        tuple val(dataset), file(frq_file), val(chrm), file(pvalue_file), file(fst_file)

    script:
        base = file(vcf_file.baseName).baseName
        frq_file = "${base}.frq"
        """
        set -o pipefail
        tabix -f ${vcf_file}

        # Inject cross-population MAFs from the annotated PGx VCF
        bcftools annotate \
            -a ${pop_mafs_vcf} \
            -c INFO/AIBST_MAF,INFO/KNL_MAF,INFO/KNK_MAF,INFO/KNM_MAF,INFO/KNP_MAF,INFO/NGI_MAF,INFO/NGH_MAF,INFO/NGY_MAF,INFO/TZA_MAF,INFO/TZB_MAF,INFO/SAV_MAF,INFO/ZWD_MAF,INFO/ZWS_MAF \
            -Oz -o annotated.vcf.gz \
            ${vcf_file}
        tabix -f annotated.vcf.gz

        echo -e '#CHRM:POS\\tCHROM\\tPOS\\tREF\\tALT\\tGENE\\tEFFECT\\t${dataset}_AC_\\t${dataset}_MAF_\\t${dataset}_AF_\\tAIBST_AF\\tAF_TGP\\tAF_joint_afr\\tAF_EXAC\\tAF_ESP\\tKNL_AF\\tKNK_AF\\tKNM_AF\\tKNP_AF\\tNGI_AF\\tNGH_AF\\tNGY_AF\\tTZA_AF\\tTZB_AF\\tSAV_AF\\tZWD_AF\\tZWS_AF\\tCLNDISDB\\tCLNDN\\tCLNSIG' > ${frq_file}
        bcftools view \
            -m2 -M2 -v snps \
            annotated.vcf.gz | \
        bcftools +fill-tags | \
        bcftools annotate \
            --set-id +"%CHROM:%POS" | \
        SnpSift extractFields - -e "0.0" \
            ID CHROM POS REF ALT ANN[0].GENE ANN[0].EFFECT AC MAF AF AIBST_MAF AF_TGP AF_joint_afr AF_EXAC AF_ESP KNL_MAF KNK_MAF KNM_MAF KNP_MAF NGI_MAF NGH_MAF NGY_MAF TZA_MAF TZB_MAF SAV_MAF ZWD_MAF ZWS_MAF CLNDISDB CLNDN CLNSIG | \
            tail -n+2 >> ${frq_file}
        """
}

process get_freq_3_1 {
    tag "frq_${dataset}_${chrm}"
    publishDir "${params.outDir}/freq/${dataset}/", mode:'copy', pattern:'*.frq'
    label "bigmem"

    input:
        tuple val(dataset), file(vcf_file), file(sample_file), val(chrm), file(pvalue_file)

    output:
        tuple val(dataset), file(vcf_file), file(sample_file), file(vcf_sites), file(frq), file(frq_1), val(chrm), file(pvalue_file)

    script:
        base = file(vcf_file.baseName).baseName
        frq = "${base}.frq"
        frq_1 = "${base}.mafs"
        vcf_sites = "${base}.sites"
            //  bcftools view \
            // -m2 -M2 -v snps | \
        """
        tabix -f ${vcf_file}
        echo -e '#CHRM:POS\\tCHROM\\tPOS\\tREF\\tALT\\tGENE\\tEFFECT\\t${dataset}_AC_\\t${dataset}_MAF_\\t${dataset}_AF_\\tAIBST_AF\\tAGVP_AF\\tKG_AF\\tKG_AFR_AF\\tgnomAD_AF\\tgnomAD_AFR_AF\\tTOPMED_AF\\tKNL_AF\\tKNK_AF\\tKNM_AF\\tKNP_AF\\tNGI_AF\\tNGH_AF\\tNGY_AF\\tTZA_AF\\tTZB_AF\\tSAV_AF\\tZWD_AF\\tZWS_AF\\tACB_AF\\tASW_AF\\tCEU_AF\\tESN_AF\\tGBR_AF\\tGWD_AF\\tLWK_AF\\tMSL_AF\\tYRI_AF\\tAmhara_AF\\tBaganda_AF\\tGumuz_AF\\tOromo_AF\\tSomali_AF\\tWolayta_AF\\tZulu_AF\\tCLNDISDB\\tCLNDN\\tCLNSIG' > ${frq}
        bcftools norm -m- \
            ${vcf_file} | \
        bcftools +fill-tags | \
        bcftools annotate \
            --set-id +"%CHROM:%POS" | \
        SnpSift extractFields - -e "0.0" \
            ID CHROM POS REF ALT ANN[0].GENE ANN[0].EFFECT AC MAF AF AIBST_AIBST_AF AGVD_AGVD_AF 1KG_1KG_AF 1KG_afr_AFRAF gnomad_b37_gnomad_b37_AF gnomad_b37_afr_AF topmed_topmed_TOPMED AIBST_KNL_AF AIBST_KNK_AF AIBST_KNM_AF AIBST_KNP_AF AIBST_NGI_AF AIBST_NGH_AF AIBST_NGY_AF AIBST_TZA_AF AIBST_TZB_AF AIBST_SAV_AF AIBST_ZWD_AF AIBST_ZWS_AF 1KG_ACB_AF 1KG_ASW_AF 1KG_CEU_AF 1KG_ESN_AF 1KG_GBR_AF 1KG_GWD_AF 1KG_LWK_AF 1KG_MSL_AF 1KG_YRI_AF AGVD_Amhara_AF AGVD_Baganda_AF AGVD_Gumuz_AF AGVD_Oromo_AF AGVD_Somali_AF AGVD_Wolayta_AF AGVD_Zulu_AF CLNDISDB CLNDN CLNSIG | \
            tail -n+2 >> ${frq}
        echo -e '#CHRM:POS\\t${dataset}_AF' > ${frq_1}
        awk -F'\\t' '{print \$2":"\$3"\\t"\$7}' ${frq} | tail -n+2 >> ${frq_1}
        echo -e '#CHRM:POS' > ${vcf_sites}
        bcftools query \
            -f "%CHROM:%POS\\n" \
            ${vcf_file} >> ${vcf_sites}  
        """
}

process get_freq_from_vcf_sites {
    tag "frq_${dataset}_${chrm}"
    publishDir "${params.outDir}/freq/${dataset}/", mode:'copy', pattern:'*.frq'
    label "bigmem"

    input:
        tuple val(dataset), file(vcf_file), file(sample_file), val(chrm)

    output:
        tuple val(dataset), file(vcf_file), file(sample_file), val(chrm), file(frq)

    script:
        base = file(vcf_file.baseName).baseName
        frq = "${base}.frq"
        """
        set -o pipefail
        tabix -f ${vcf_file}
        echo -e "#CHRM:POS\\tCHROM\\tPOS\\tREF\\tALT\\tGENE\\tEFFECT\\tAIBST_AF\\tAIBST_AC\\tAF_joint_afr\\tAF_TGP\\tAF_EXAC\\tAF_ESP\\tCLNDISDB\\tCLNDN\\tCLNSIG\\tAF_gnomAD\\tAF_gnomAD_exomes\\tAF_gnomAD_genomes\\tAF_gnomAD_AFR\\tAF_gnomAD_AMR\\tAF_gnomAD_EAS\\tAF_gnomAD_NFE\\tAF_gnomAD_SAS" > ${frq}
        bcftools norm -m- \
            ${vcf_file} | \
        bcftools annotate \
            --set-id +"%CHROM:%POS" | \
        SnpSift extractFields - -e "0.0" \
            ID CHROM POS REF ALT ANN[0].GENE ANN[0].EFFECT AF AC AF_joint_afr AF_TGP AF_EXAC AF_ESP CLNDISDB CLNDN CLNSIG AF_joint AF_exomes AF_genomes AF_joint_afr AF_joint_amr AF_joint_eas AF_joint_nfe AF_joint_sas | \
            tail -n+2 >> ${frq}
        """
}

process get_freq_3_2 {
    tag "frq_${dataset}_${chrm}"
    publishDir "${params.outDir}/freq/${dataset}/", mode:'copy', pattern:'*.frq'
    label "bigmem"

    input:
        tuple val(dataset), file(vcf_file), file(sample_file), val(chrm)

    output:
        tuple val(dataset), file(vcf_file), file(sample_file), val(chrm), file(frq)

    script:
        base = file(vcf_file.baseName).baseName
        frq = "${base}.frq"
        """
        set -o pipefail
        tabix -f ${vcf_file}
        echo -e "#CHRM:POS\\tCHROM\\tPOS\\tGENE\\tEFFECT\\t${dataset}_AC_\\t${dataset}_MAF_\\tAIBST_MAF\\tAF_TGP\\tAF_joint_afr\\tAF_EXAC\\tAF_ESP\\tKNL_MAF\\tKNK_MAF\\tKNM_MAF\\tKNP_MAF\\tNGI_MAF\\tNGH_MAF\\tNGY_MAF\\tTZA_MAF\\tTZB_MAF\\tSAV_MAF\\tZWD_MAF\\tZWS_MAF\\tCLNDISDB\\tCLNDN\\tCLNSIG\\tREF\\tALT\\tAF_gnomAD\\tAF_gnomAD_exomes\\tAF_gnomAD_genomes\\tAF_gnomAD_AFR\\tAF_gnomAD_AMR\\tAF_gnomAD_EAS\\tAF_gnomAD_NFE\\tAF_gnomAD_SAS" > ${frq}
        bcftools norm -m- \
            ${vcf_file} | \
        bcftools +fill-tags | \
        bcftools annotate \
            --set-id +"%CHROM:%POS" | \
        SnpSift extractFields - -e "0.0" \
            ID CHROM POS ANN[0].GENE ANN[0].EFFECT AC MAF AIBST_MAF AF_TGP AF_joint_afr AF_EXAC AF_ESP KNL_MAF KNK_MAF KNM_MAF KNP_MAF NGI_MAF NGH_MAF NGY_MAF TZA_MAF TZB_MAF SAV_MAF ZWD_MAF ZWS_MAF CLNDISDB CLNDN CLNSIG REF ALT AF_joint AF_exomes AF_genomes AF_joint_afr AF_joint_amr AF_joint_eas AF_joint_nfe AF_joint_sas | \
            tail -n+2 >> ${frq}
        """
}


process get_sites_simple {
    tag "sites_${dataset}_${chrm}"
    label "bigmem"

    input:
        tuple val(dataset), file(vcf_file), file(novel_sites), val(chrm)

    output:
        tuple val(dataset), file(vcf_file), file(vcf_sites), val(chrm)

    script:
        base = file(vcf_file.baseName).baseName
        vcf_sites = "${dataset}_${chrm}.sites"
        """
        tabix ${vcf_file}
        echo -e '#CHRM:POS' > ${vcf_sites}
        awk '{print \$1":"\$2}' ${novel_sites} >> ${vcf_sites}
        """
}

process data_for_upset{
    tag "freq_for_upset_${dataset}"
    label "bigmem"

    input:
      tuple val(dataset), val(site_files)

    output:
      tuple val(dataset), file(csv_out)

    script:
        sites = site_files.join(',')
        csv_out = "${dataset}.upset.csv"
        template "data_for_upset.py"
}

process plot_upset_R{
    tag "plot_upset_R_${dataset}"
    label "bigmem"
    publishDir "${params.outDir}/pgx", mode: 'copy'

    input:
        tuple val(dataset), file(csv_file)

    output:
        tuple val(dataset), file(upset_plot)

    script:
        upset_plot = "${dataset}.upset.pdf"
        template "novel_upset.R"
}

process get_freq1 {
    tag "frq_${dataset}_${chrm}"
    label "bigmem"

    input:
        tuple val(dataset), file(vcf_file), file(sample_file), val(chrm)

    output:
        tuple val(dataset), file(frq), val(chrm)
        tuple val(dataset), file(vcf_file), file(sample_file), file(vcf_sites), file(frq_1), val(chrm)

    script:
        base = file(vcf_file.baseName).baseName
        frq = "${base}.frq"
        frq_1 = "${base}.mafs"
        vcf_sites = "${base}.sites"
        """
        tabix ${vcf_file}
        echo -e 'ID\\tCHROM\\tPOS\\tGENE\\tEFFECT\\tCLNDISDB\\tCLNDN\\tCLNSIG\\tAC\\tMAF\\tAGVP_AF\\tSAHGP_AF\\tKG_AF\\tKG_AFR_AF\\tExAC_AF\\tExAC_AFR_AF\\tgnomAD_AF\\tgnomAD_AFR_AF\\tTOPMEDKNK_AF\\tKNL_AF\\tKNM_AF\\tKNP_AF\\tNGH_AF\\tNGI_AF\\tNGY_AF\\tSAV_AF\\tTZA_AF\\tTZB_AF\\tZWD_AF\\tZWS_AF' > ${frq}
        bcftools view \
            -m2 -M2 -v snps \
            ${vcf_file} | \
        bcftools +fill-tags | \
        bcftools annotate \
            --set-id +"%CHROM\\_%POS\\_%REF\\_%FIRST_ALT" | \
        SnpSift extractFields - -e "0.0" \
            ID CHROM POS ANN[0].GENE ANN[0].EFFECT CLNDISDB CLNDN CLNSIG AC MAF AGVP_AF SAHGP_AF KG_AF KG_AFR_AF ExAC_AF ExAC_AFR_AF gnomAD_AF gnomAD_AFR_AF TOPMED KNK_AF KNL_AF KNM_AF KNP_AF NGH_AF NGI_AF NGY_AF SAV_AF TZA_AF TZB_AF ZWD_AF ZWS_AF | \
            tail -n+2 >> ${frq}
        echo -e '#CHRM:POS\\t${dataset}_AF' > ${frq_1}
        awk -F'\\t' '{print \$2":"\$3"\\t"\$10}' ${frq} | tail -n+2 >> ${frq_1}
        echo -e '#CHRM:POS' > ${vcf_sites}
        bcftools query \
            -f "%CHROM:%POS\\n" \
            ${vcf_file} >> ${vcf_sites}
        """
}

process cat_freq {
    tag "cat_frq_${dataset}"
    label "bigmem"
    publishDir "${params.outDir}/freq/${dataset}/", mode: 'copy'

    input:
        tuple val(dataset), val(frqs), val(chrms)

    output:
        tuple val(dataset), file(frq)

    script:
        frq = "${dataset}.frq"

        """
        head -n 1 ${frqs[0]} > ${frq}
        tail -q -n +2 ${frqs.join(' ')} >> ${frq}
        """
        // head -n 1 ${frqs[0]} > ${frq}
        // tail -q -n +2 ${frqs.join(' ')} | \
        // awk -F"\\t" '{OFS="\\t"; n = split( \$1, a, ";" ); \$1=a[1]; print \$0}' >> ${frq}
}

process cat_sites {
    tag "cat_frq_${dataset}"
    label "bigmem"
    publishDir "${params.outDir}/freq/${dataset}/", mode: 'copy'

    input:
        tuple val(dataset), val(sites), val(chrms)

    output:
        tuple val(dataset), file(site)

    script:
        site = "${dataset}.sites"
        """
        echo -e '#CHRM:POS\\t${dataset}' > ${site}
        tail -q -n +2 ${sites.join(' ')} | awk '{print \$1"\\t1"}' >> ${site}
        """
}

process fisher_test_plink {
    tag "fisher_test_plink_${dataset}_${chrm}"
    label "bigmem"
    publishDir "${params.outDir}/fisher-test", mode:'copy', overwrite: true

    input:
        tuple val(dataset), file(plink_bed), file(plink_bim), file(plink_fam), file(pheno), val(chrm), val(prefix)
    output:
        tuple val(dataset), file("${out}.assoc.fisher"), val(chrm)
    script:
        out = "${dataset}_${chrm}_${prefix}"
        """
        plink \
            --bfile ${plink_bed.baseName} \
            --pheno ${pheno} \
            --allow-no-sex \
            --assoc fisher \
            --out ${out}
        """
}

process fisher_test_plink_pops {
    tag "fisher_test_plink_${dataset}_${chrm}"
    label "bigmem"
    publishDir "${params.outDir}/fisher-test", mode:'copy', overwrite: true

    input:
        tuple val(dataset), file(plink_bed), file(plink_bim), file(plink_fam), file(pheno), val(chrm), val(prefix)
    output:
        tuple val(dataset), file(out), val(chrm)
    script:
        base = "${dataset}_${chrm}_${prefix}"
        out = "${base}.assoc.fisher.csv"
        """
        plink \
            --bfile ${plink_bed.baseName} \
            --pheno ${pheno} \
            --allow-no-sex \
            --assoc fisher \
            --out ${base}
        awk 'NR==1{print "POPS\\t"\$0}' ${base}.assoc.fisher > ${out}
        awk 'NR>1{print "${dataset}\\t"\$0}' ${base}.assoc.fisher >> ${out}
        """
}

process filter_fisher_test{
    tag "filter_fisher_${dataset}"
    label "bigmem"

    input:
      tuple val(dataset), file(fisher), val(chrm)

    output:
      tuple val(dataset), file(out), val(chrm)

    script:
        out = "${fisher.baseName}.assoc.fisher.filtered.csv"
        template "filter_fisher_test.py"
}

process add_pops_fisher_test {
    tag "fisher_test_plink_${dataset}_${chrm}"
    label "bigmem"
    publishDir "${params.outDir}/fisher-test", mode:'copy', overwrite: true

    input:
        tuple val(dataset), file(fisher_file), val(chrm)
    output:
        tuple val(dataset), file(out), val(chrm)
    script:
        out = "${fisher_file}.csv"
        """
        awk 'NR==1{print "POPS\\t"\$0}' ${fisher_file} > ${out}
        awk 'NR>1{print "${dataset}\\t"\$0}' ${fisher_file} >> ${out}
        """
}

process get_map_vcf {
    tag "get_map_${dataset}"
    label "bigmem"

    input:
        tuple val(dataset), file(dataset_vcf), file(dataset_vcf_index)
    output:
        tuple val(dataset), file(dataset_vcf), file(dataset_vcf_index), file(dataset_map)
    script:
        base = file(dataset_vcf.baseName).baseName
        dataset_map = "${base}.map"
        """
        bcftools query -f '%CHROM\\n' ${dataset_vcf} | sort -n | uniq > ${dataset_map}
        """
}

process get_map {
    tag "get_map_${dataset}"
    label "bigmem"

    input:
        tuple val(dataset), file(dataset_vcf), file(vcf_idx), file(dataset_sample)
    output:
        tuple val(dataset), file(dataset_vcf), file(vcf_idx), file(dataset_sample), file(dataset_map)
    script:
        base = file(dataset_vcf.baseName).baseName
        dataset_map = "${base}.map"
        """
        bcftools query -f '%CHROM\\t%POS\\n' ${dataset_vcf} > ${dataset_map}
        """
}

process add_clinvar_to_freq {
    tag "add_clinvar_to_freq_${dataset}"
    publishDir "${params.outDir}/freq/pgx_only/${dataset}", mode:'copy', overwrite: true
    label "bigmem"

    input:
        tuple val(dataset), file(freq_file), file(clinvar_file), file(pharmgkb_file), file(gwas_file)

    output:
        tuple val(dataset), file(freq_out)

    script:
        freq_out = "${freq_file.getSimpleName()}_clinvar.csv"
        template "add_clinvar_to_freq.py"
}

/*
 * Annotate per-population MAFs into the VCF.
 * Computes AF for each population using sample subsets,
 * then adds {POP}_MAF and AIBST_MAF INFO fields to the VCF.
 */
process annotate_pop_mafs {
    tag "pop_mafs_${dataset}"
    label "bigmem"

    input:
        tuple val(dataset), file(vcf_file), file(sample_file)

    output:
        tuple val(dataset), file("${base}_pop_mafs.vcf.gz"), file(sample_file)

    script:
        base = file(vcf_file.baseName).baseName
        """
        set -euo pipefail

        # Index input if needed
        [ -f ${vcf_file}.tbi ] || tabix -f ${vcf_file}

        # Step 1: Add overall AIBST_MAF (= AF from all samples)
        bcftools +fill-tags ${vcf_file} -- -t AF,AC,AN | \
            bcftools annotate --rename-chrs /dev/null -Oz -o step0.vcf.gz 2>/dev/null || \
            bcftools +fill-tags ${vcf_file} -- -t AF,AC,AN -Oz -o step0.vcf.gz
        tabix -f step0.vcf.gz

        # Rename AF -> AIBST_MAF using bcftools annotate
        echo '##INFO=<ID=AIBST_MAF,Number=A,Type=Float,Description="AIBST overall allele frequency">' > aibst_hdr.txt
        bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AF\\n' step0.vcf.gz | bgzip > aibst_af.tsv.gz
        tabix -s1 -b2 -e2 aibst_af.tsv.gz
        bcftools annotate -a aibst_af.tsv.gz -h aibst_hdr.txt \
            -c CHROM,POS,REF,ALT,AIBST_MAF step0.vcf.gz -Oz -o step1.vcf.gz
        tabix -f step1.vcf.gz

        # Step 2: For each population, compute AF and annotate
        POPS=\$(awk -F'\\t' 'NR>1{print \$2}' ${sample_file} | sort -u)
        cp step1.vcf.gz current.vcf.gz
        cp step1.vcf.gz.tbi current.vcf.gz.tbi

        for POP in \$POPS; do
            # Extract sample list for this population
            awk -F'\\t' -v pop="\$POP" 'NR>1 && \$2==pop {print \$1}' ${sample_file} > \${POP}_samples.txt

            # Skip if no samples
            [ -s \${POP}_samples.txt ] || continue

            # Compute per-pop AF
            bcftools view -S \${POP}_samples.txt current.vcf.gz | \
                bcftools +fill-tags -- -t AF | \
                bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AF\\n' | \
                bgzip > \${POP}_af.tsv.gz
            tabix -s1 -b2 -e2 \${POP}_af.tsv.gz

            # Create header for this pop
            printf '##INFO=<ID=%s_MAF,Number=A,Type=Float,Description="%s population allele frequency">\n' \${POP} \${POP} > \${POP}_hdr.txt

            # Annotate
            bcftools annotate -a \${POP}_af.tsv.gz -h \${POP}_hdr.txt \
                -c CHROM,POS,REF,ALT,\${POP}_MAF current.vcf.gz -Oz -o next.vcf.gz
            tabix -f next.vcf.gz
            mv next.vcf.gz current.vcf.gz
            mv next.vcf.gz.tbi current.vcf.gz.tbi
        done

        mv current.vcf.gz ${base}_pop_mafs.vcf.gz
        tabix -f ${base}_pop_mafs.vcf.gz
        """
}
