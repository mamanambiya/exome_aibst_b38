cd /scratch3/users/mamana/exome_aibst/work/52/69abe9880c0e1888a4829bae54040a


# Get missense variant - genes
zcat AIBST.vcf.gz | SnpSift filter "(ANN[*].EFFECT has 'missense_variant')" | bgzip -c > AIBST.missense.vcf.gz
# Count missense genes
zcat AIBST.missense.vcf.gz | SnpSift extractFields - -e "." ANN[0].GENE | uniq | wc -l 
zcat AIBST.missense.vcf.gz | SnpSift extractFields - -e "." ID CHROM POS AF AC AN ANN[0].GENE ANN[0].EFFECT CLNDISDB CLNDN CLNSIG KNL_MAF KNK_MAF KNM_MAF KNP_MAF NGI_MAF NGH_MAF NGY_MAF TZA_MAF TZB_MAF SAV_MAF ZWD_MAF ZWS_MAF ACB_MAF ASW_MAF CEU_MAF ESN_MAF GBR_MAF GWD_MAF LWK_MAF MSL_MAF YRI_MAF Zulu_MAF > AIBST.missense.gene.csv

while read GENE; do
    echo $GENE:$(zcat AIBST.missense.vcf.gz | grep $GENE | wc -l)
done < AIBST.missense.gene.csv > AIBST.missense.gene.count.csv


# In VIP
while read GENE; do
    echo $(cat ~/exome_aibst/PGX_DATA/2021/VIP_gene.csv | grep $GENE)
done < AIBST.missense.gene.csv | uniq > AIBST.missense.VIP.count.csv




## LOF
zcat AIBST.vcf.gz | SnpSift filter "(LOF[*].PERC >= 0.5)" | bgzip -c > AIBST.lof.vcf.gz
zcat AIBST.lof.vcf.gz | SnpSift extractFields - -e "." ID CHROM POS AF AC AN ANN[0].GENE ANN[0].EFFECT CLNDISDB CLNDN CLNSIG KNL_MAF KNK_MAF KNM_MAF KNP_MAF NGI_MAF NGH_MAF NGY_MAF TZA_MAF TZB_MAF SAV_MAF ZWD_MAF ZWS_MAF ACB_MAF ASW_MAF CEU_MAF ESN_MAF GBR_MAF GWD_MAF LWK_MAF MSL_MAF YRI_MAF Zulu_MAF > AIBST.lof.gene.csv
# Per gene
zcat AIBST.vcf.gz | SnpSift filter "(ANN[*].GENE= '$G') && (LOF[*].PERC >= 0)" | grep -v "^#" | wc -l
zcat AIBST.vcf.gz | SnpSift filter "(ANN[*].GENE= '$G') && (ANN[*].EFFECT has 'missense_variant')" | grep -v "^#" | wc -l


for G in HLA-B CYP2D6  SLCO1B1  TPMT NUDT15  DPYD  HLA-A  CYP2C19  CYP2B6  CYP2C9  NAT  UGT1A1  CYP4F2  VKORC1
do
    echo $G
    zcat AIBST.vcf.gz | SnpSift filter "(ANN[*].GENE= '$G')" | grep -v "^#" | wc -l
    zcat AIBST.vcf.gz | SnpSift filter "(ANN[*].GENE= '$G') && (LOF[*].PERC >= 0)" | grep -v "^#" | wc -l
    zcat AIBST.vcf.gz | SnpSift filter "(ANN[*].GENE= '$G') && (ANN[*].EFFECT has 'missense_variant')" | grep -v "^#" | wc -l
done


CYP2D6

zcat AIBST.vcf.gz | SnpSift filter "(AC>3)" | SnpSift extractFields - -e "." ID CHROM POS AF AC AN ANN[0].GENE ANN[0].EFFECT CLNDISDB CLNDN CLNSIG KNL_MAF KNK_MAF KNM_MAF KNP_MAF NGI_MAF NGH_MAF NGY_MAF TZA_MAF TZB_MAF SAV_MAF ZWD_MAF ZWS_MAF ACB_MAF ASW_MAF CEU_MAF ESN_MAF GBR_MAF GWD_MAF LWK_MAF MSL_MAF YRI_MAF Zulu_MAF > Novel_ac_3.csv

zcat AIBST.vcf.gz | SnpSift filter "(LOF[*].PERC >= 0)"
zcat AIBST.vcf.gz | SnpSift filter "(ANN[*].EFFECT has 'missense_variant')" | grep -v "#" | wc -l
zcat AIBST.vcf.gz | SnpSift filter "(AC>3)" | grep -v "#" | wc -l



Rifampicin, Efavirenz, Nevirapine, Nelfinavir, Atazanavir, Ritonavir, Tenofovir, Isoniazid, Rifampin, Pyrazinamide, Ethambutol, Hydroxychloroquine, Tafenoquine, Primaquine, Amodiaquine, Lumefantrine




