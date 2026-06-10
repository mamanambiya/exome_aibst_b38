#!/usr/bin/env python3.7
'''

'''
import argparse
import sys
import time
import pandas as pd
import csv


parser = argparse.ArgumentParser()
parser.add_argument(
    "--vcf_sites_bed", default="${vcf_sites_bed}", help="vcf")
parser.add_argument("--novel_sites_file",
                    default="${novel_sites_file}", help="Novel sites output file")

args = parser.parse_args()

def get_novel(vcf_sites_bed, novel_sites_file):
    """
    :param vcf_sites_bed:
    :param novel_sites_file:
    :return:
    """
    sites_out = open(novel_sites_file, 'w')
    a = pd.read_csv(vcf_sites_bed, delimiter='\\\\s+',
                    quotechar='\\"', engine='python')
    datas = a.to_dict()
    header = datas.keys()
    ids = datas['#CHRM:POS'].keys()
    for idx in ids:
        freqs = []
        chrm = str(datas['CHROM'][idx])
        pos = str(datas['POS'][idx])
        for col in header:
            if '_AF' in col or 'TOPMED' in col:
                allele_frq = datas[col][idx]
                freqs.append(float(allele_frq))
        if all(frq == 0 for frq in [float(frq) for frq in freqs]):
            data = '\\t'.join([chrm, pos, pos])+"\\n"
            sites_out.writelines(data)
    sites_out.close()

if __name__ == '__main__':
    get_novel(args.vcf_sites_bed, args.novel_sites_file)
