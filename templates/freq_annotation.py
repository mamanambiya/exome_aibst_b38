#!/usr/bin/env python3.7
'''

'''
import argparse,sys,time
import pandas as pd
import numpy as np
import csv

parser = argparse.ArgumentParser()
parser.add_argument("--inTSV", default="${inTSV}", help="One or many TSV annotation files, if many use ';' as separator")
parser.add_argument("--outAnnot", default="${outAnnot}", help="TSV annotation file")
parser.add_argument("--outHdr", default="${outHdr}", help="TSV annotation header file")
parser.add_argument("--sites", default="${sites}", help="Base annotation file")
parser.add_argument("--annot", default="${annot}", action='store_true', help="Annotated TSV file")
args = parser.parse_args()

datasets = {
             'KNL_AF': 'AIBST',
             'KNP_AF': 'AIBST',
             'ZWS_AF': 'AIBST',
             'NGI_AF': 'AIBST',
             'NGH_AF': 'AIBST',
             'KNK_AF': 'AIBST',
             'KNM_AF': 'AIBST',
             'SAV_AF': 'AIBST',
             'TZB_AF': 'AIBST',
             'TZA_AF': 'AIBST',
             'NGY_AF': 'AIBST',
             'ZWD_AF': 'AIBST',
             'TOPMED': 'TOPMED',
             'AGVP_AF': 'AGVP',
             'SAHGP_AF': 'SAHGP',
             'SAN_AF': 'SAN',
             'TRYPANOGEN_AF': 'TRYPANOGEN',
             'KG_AC': 'KG', 'KG_AF': 'KG', 'KG_AFR_AC': 'KG', 'KG_AFR_AF': 'KG', 'KG_EUR_AC': 'KG', 'KG_EUR_AF': 'KG', 'KG_AMR_AC': 'KG', 'KG_AMR_AF': 'KG', 'KG_EAS_AC': 'KG', 'KG_EAS_AF': 'KG', 'KG_SAS_AC': 'KG', 'KG_SAS_AF': 'KG',
             'TWINSUK_AC': 'TWINSUK', 'TWINSUK_AF': 'TWINSUK',
             'ALSPAC_AC': 'ALSPAC', 'ALSPAC_AF': 'ALSPAC',
             'ESP6500_AA_AC': 'ESP6500', 'ESP6500_AA_AF': 'ESP6500', 'ESP6500_EA_AC': 'ESP6500', 'ESP6500_EA_AF': 'ESP6500',
             'Wolayta_AF': 'AGVP', 'Amhara_AF': 'AGVP' ,'Baganda_AF': 'AGVP' ,'Gumuz_AF': 'AGVP' ,'Oromo_AF': 'AGVP' ,'Somali_AF': 'AGVP' ,'Wolayta_AF': 'AGVP' ,'Zulu_AF': 'AGVP',
             'ExAC_AC': 'ExAC', 'ExAC_AF': 'ExAC', 'ExAC_Adj_AC': 'ExAC', 'ExAC_Adj_AF': 'ExAC', 'ExAC_AFR_AC': 'ExAC', 'ExAC_AFR_AF': 'ExAC', 'ExAC_AMR_AC': 'ExAC', 'ExAC_AMR_AF': 'ExAC', 'ExAC_EAS_AC': 'ExAC', 'ExAC_EAS_AF': 'ExAC', 'ExAC_FIN_AC': 'ExAC', 'ExAC_FIN_AF': 'ExAC', 'ExAC_NFE_AC': 'ExAC', 'ExAC_NFE_AF': 'ExAC', 'ExAC_SAS_AC': 'ExAC', 'ExAC_SAS_AF': 'ExAC', 'ExAC_nonTCGA_AC': 'ExAC', 'ExAC_nonTCGA_AF': 'ExAC', 'ExAC_nonTCGA_Adj_AC': 'ExAC', 'ExAC_nonTCGA_Adj_AF': 'ExAC', 'ExAC_nonTCGA_AFR_AC': 'ExAC', 'ExAC_nonTCGA_AFR_AF': 'ExAC', 'ExAC_nonTCGA_AMR_AC': 'ExAC', 'ExAC_nonTCGA_AMR_AF': 'ExAC', 'ExAC_nonTCGA_EAS_AC': 'ExAC', 'ExAC_nonTCGA_EAS_AF': 'ExAC', 'ExAC_nonTCGA_FIN_AC': 'ExAC', 'ExAC_nonTCGA_FIN_AF': 'ExAC', 'ExAC_nonTCGA_NFE_AC': 'ExAC', 'ExAC_nonTCGA_NFE_AF': 'ExAC', 'ExAC_nonTCGA_SAS_AC': 'ExAC', 'ExAC_nonTCGA_SAS_AF': 'ExAC', 'ExAC_nonpsych_AC': 'ExAC', 'ExAC_nonpsych_AF': 'ExAC', 'ExAC_nonpsych_Adj_AC': 'ExAC', 'ExAC_nonpsych_Adj_AF': 'ExAC', 'ExAC_nonpsych_AFR_AC': 'ExAC', 'ExAC_nonpsych_AFR_AF': 'ExAC', 'ExAC_nonpsych_AMR_AC': 'ExAC', 'ExAC_nonpsych_AMR_AF': 'ExAC', 'ExAC_nonpsych_EAS_AC': 'ExAC', 'ExAC_nonpsych_EAS_AF': 'ExAC', 'ExAC_nonpsych_FIN_AC': 'ExAC', 'ExAC_nonpsych_FIN_AF': 'ExAC', 'ExAC_nonpsych_NFE_AC': 'ExAC', 'ExAC_nonpsych_NFE_AF': 'ExAC', 'ExAC_nonpsych_SAS_AC': 'ExAC', 'ExAC_nonpsych_SAS_AF': 'ExAC',
             'gnomAD_exomes_AC':'gnomAD', 'gnomAD_exomes_AN':'gnomAD', 'gnomAD_exomes_AF':'gnomAD', 'gnomAD_exomes_AFR_AC':'gnomAD', 'gnomAD_exomes_AFR_AN':'gnomAD', 'gnomAD_exomes_AFR_AF':'gnomAD', 'gnomAD_exomes_AMR_AC':'gnomAD', 'gnomAD_exomes_AMR_AN':'gnomAD', 'gnomAD_exomes_AMR_AF':'gnomAD', 'gnomAD_exomes_ASJ_AC':'gnomAD', 'gnomAD_exomes_ASJ_AN':'gnomAD', 'gnomAD_exomes_ASJ_AF':'gnomAD', 'gnomAD_exomes_EAS_AC':'gnomAD', 'gnomAD_exomes_EAS_AN':'gnomAD', 'gnomAD_exomes_EAS_AF':'gnomAD', 'gnomAD_exomes_FIN_AC':'gnomAD', 'gnomAD_exomes_FIN_AN':'gnomAD', 'gnomAD_exomes_FIN_AF':'gnomAD', 'gnomAD_exomes_NFE_AC':'gnomAD', 'gnomAD_exomes_NFE_AN':'gnomAD', 'gnomAD_exomes_NFE_AF':'gnomAD', 'gnomAD_exomes_SAS_AC':'gnomAD', 'gnomAD_exomes_SAS_AN':'gnomAD', 'gnomAD_exomes_SAS_AF':'gnomAD', 'gnomAD_exomes_OTH_AC':'gnomAD', 'gnomAD_exomes_OTH_AN':'gnomAD', 'gnomAD_exomes_OTH_AF':'gnomAD', 'gnomAD_AC':'gnomAD', 'gnomAD_AN':'gnomAD', 'gnomAD_AF':'gnomAD', 'gnomAD_AFR_AC':'gnomAD', 'gnomAD_AFR_AN':'gnomAD', 'gnomAD_AFR_AF':'gnomAD', 'gnomAD_AMR_AC':'gnomAD', 'gnomAD_AMR_AN':'gnomAD', 'gnomAD_AMR_AF':'gnomAD', 'gnomAD_ASJ_AC':'gnomAD', 'gnomAD_ASJ_AN':'gnomAD', 'gnomAD_ASJ_AF':'gnomAD', 'gnomAD_EAS_AC':'gnomAD', 'gnomAD_EAS_AN':'gnomAD', 'gnomAD_EAS_AF':'gnomAD', 'gnomAD_FIN_AC':'gnomAD', 'gnomAD_FIN_AN':'gnomAD', 'gnomAD_FIN_AF':'gnomAD', 'gnomAD_NFE_AC':'gnomAD', 'gnomAD_NFE_AN':'gnomAD', 'gnomAD_NFE_AF':'gnomAD', 'gnomAD_OTH_AC':'gnomAD', 'gnomAD_OTH_AN':'gnomAD', 'gnomAD_OTH_AF':'gnomAD'}

annots = {
          "AIBST":"Allele Frequency from AIBST Project (AIBST)",
          "KG":"Allele Frequency from Thousand Genomes Project (1KG)",
          "gnomAD":"Allele Frequency from the Genome Aggregation Database (gnomAD)",
          "ExAC":"Allele Frequency from the Exome Aggregation Consortium (ExAC)",
          "ESP6500":"Allele Frequency from the NHLBI GO Exome Sequencing Project (ESP)",
          "ALSPAC":"Allele Frequency from the Avon Longitudinal Study of Parents and Children (ALSPAC)",
          "TWINSUK":"Allele Frequency from the TwinsUK project",
          "AGVP": "Allele Frequency among the African Genome Variation Project (AGVP)",
          "SAHGP": "Allele Frequency from the Southern African Human Genome Project (SAHGP)",
          "TRYPANOGEN": "Allele Frequency from H3Africa Trypanogen Project",
          "TOPMED": "Allele Frequency from TOPMED Project",
          "SAN": "Allele frequency of SAN population from Namibia "}

def readTSV(inTSV_list, sites='', outAnnot='', outHdr=''):
    '''
    :param inBed:
    :param outBed:
    :return:
    '''
    inTSV_list = inTSV_list.split(';')

    base = pd.read_csv(sites, sep='\\\\s+', engine='python', error_bad_lines=False)

    for annot in inTSV_list:
        data = pd.read_csv(annot, sep='\\\\s+', engine='python', error_bad_lines=False)
        base = pd.merge(base, data, on='#CHRM:POS', how='left')
    base = base.replace(np.nan, '0.0', regex=True)
    base = base.replace('.', '0.0')

    datas = base.to_dict()  # Transform to dict
    for col in datas:
        if '#CHRM:POS' not in col:
            for idx in datas[col]:
                frq = str(datas[col][idx]).strip().split(',')[0]
                if frq == '.':
                    frq = '0.0'
                try:
                    frq = float(frq)
                except:
                    print(frq)
                if frq > 0.5:
                    frq = 1 - frq
                datas[col][idx] = format(frq, '.5f')

    base = pd.DataFrame.from_dict(datas)

    chrom = base['#CHRM:POS'].str.split(':').str[0]
    pos = base['#CHRM:POS'].str.split(':').str[1]
    base = base.drop(columns=['#CHRM:POS'])
    base.insert(0, 'POS', pos, allow_duplicates=True)
    base.insert(0, 'POS', pos, allow_duplicates=True)
    base.insert(0, 'CHROM', chrom)

    # Writing to file
    base.to_csv(outAnnot, sep='\\t', index=False)

    info = pd.DataFrame([], columns=['info'])
    for col in base.columns[3:]:
        if col not in datasets:
            info = info.append({'info': "##INFO=<ID={},Number=1,Type=Float,Description=\"Allele frequency of {}\">".format(col, col)}, ignore_index=True)
        else:
            info = info.append({'info': "##INFO=<ID={},Number=1,Type=Float,Description=\"{}\">".format(col, annots[datasets[col]])}, ignore_index=True)

    #  Writing header
    info.to_csv(outHdr, sep=';', index=False, header=False, quoting=csv.QUOTE_NONE, quotechar='\\0')
    
if args.annot:
    readTSV(args.inTSV, args.sites, args.outAnnot, args.outHdr)


