#!/usr/bin/env python3

import argparse
import pandas as pd
import re

parser = argparse.ArgumentParser()
parser.add_argument("--freq_file", default="${freq_file}", help="")
parser.add_argument("--clinvar_file", default="${clinvar_file}", help="")
parser.add_argument("--pharmgkb_file", default="${pharmgkb_file}", help="")
parser.add_argument("--gwas_file", default="${gwas_file}", help="")
parser.add_argument("--freq_out", default="${freq_out}", help="")
args = parser.parse_args()


def parse_vcf_info(info_str, key):
    """Extract value from VCF INFO field"""
    match = re.search(f'{key}=([^;]+)', info_str)
    return match.group(1) if match else None


def add_clinvar_to_freq(freq_file, clinvar_file, pharmgkb_file, gwas_file, freq_out):
    freq = pd.read_csv(freq_file, sep='\\s+', engine='python')

    # Read ClinVar VCF (skip header lines starting with ##)
    clinvar = pd.read_csv(
        clinvar_file,
        sep='\\t',
        comment='#',  # Skip ## lines
        compression='gzip',
        names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'],
        dtype={'CHROM': str, 'POS': str, 'ID': str}
    )

    # Parse INFO field to extract CLNSIG and CLNDN
    clinvar['CLNSIG'] = clinvar['INFO'].apply(lambda x: parse_vcf_info(x, 'CLNSIG'))
    clinvar['CLNDN'] = clinvar['INFO'].apply(lambda x: parse_vcf_info(x, 'CLNDN'))

    # Keep only necessary columns
    clinvar1 = clinvar[['CHROM', 'POS', 'CLNSIG', 'CLNDN']].copy()

    # Process frequency data
    freq['CHROM'] = freq['CHROM'].astype(str)
    freq['POS'] = freq['POS'].astype(str)
    freq = freq.assign(
        var1=freq['#CHRM:POS'].str.split(';')).explode('var1')
    freq['#CHRM:POS'] = freq['var1']
    freq = freq.drop(['var1'], axis=1)

    # Merge with ClinVar
    clinvar_freq = pd.merge(
        freq, clinvar1, on=['CHROM', 'POS'], how='left', suffixes=('', '_y'))

    # Add PharmGKB annotations
    pharmgkb = pd.read_csv(pharmgkb_file, sep='\\t', engine='python')
    pharmgkb['#CHRM:POS'] = pharmgkb['Variant/Haplotypes'].astype(str)
    pharmgkb = pharmgkb[['#CHRM:POS', 'Level of Evidence',
                        'Phenotype Category', 'Drug(s)', 'Phenotype(s)']]

    clinvar_freq = pd.merge(
        clinvar_freq, pharmgkb, on=['#CHRM:POS'], how='left'
    ).drop_duplicates().reset_index(drop=True).sort_values(
        by=['CHROM', 'POS'], ascending=True
    )

    # Add GWAS annotations
    gwas = pd.read_csv(gwas_file, sep='\\t')
    gwas['#CHRM:POS'] = gwas['SNPS'].astype(str)
    gwas['GWAS_STUDY'] = gwas['STUDY'].astype(str)
    gwas['GWAS_DISEASE'] = gwas['DISEASE/TRAIT'].astype(str)
    gwas['GWAS_P-VALUE'] = gwas['P-VALUE'].astype(str)
    gwas['GWAS_RISK_ALLELE_FREQUENCY'] = gwas['RISK ALLELE FREQUENCY'].astype(str)
    gwas = gwas[['#CHRM:POS', 'GWAS_STUDY', 'GWAS_DISEASE',
                'GWAS_P-VALUE', 'GWAS_RISK_ALLELE_FREQUENCY']]
    gwas = gwas.assign(
        var1=gwas['#CHRM:POS'].str.split(';')).explode('var1')
    gwas['#CHRM:POS'] = gwas['var1']
    gwas = gwas.drop(['var1'], axis=1)

    clinvar_freq = pd.merge(
        clinvar_freq, gwas, on=['#CHRM:POS'], how='left'
    ).drop_duplicates().reset_index(drop=True).sort_values(
        by=['CHROM', 'POS'], ascending=True
    )

    # Final cleanup
    clinvar_freq = clinvar_freq.assign(
        var1=clinvar_freq['#CHRM:POS'].str.split(';')).explode('var1')
    clinvar_freq['#CHRM:POS'] = clinvar_freq['var1']
    clinvar_freq = clinvar_freq.drop(['var1'], axis=1)
    clinvar_freq['CHROM'] = clinvar_freq['CHROM'].str.replace('chr', '')

    clinvar_freq.to_csv(f'{freq_out}', index=False, header=True, sep='\\t')


if __name__ == '__main__':
    add_clinvar_to_freq(
        args.freq_file, args.clinvar_file,
        args.pharmgkb_file, args.gwas_file, args.freq_out
    )
