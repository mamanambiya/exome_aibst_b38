#!/usr/bin/env python3

import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--hdv", default="${hdv}", help="")
parser.add_argument("--annot", default="${annot}", help="")
parser.add_argument("--hdv_out", default="${hdv_out}", help="")
args = parser.parse_args()


def annot_hdv(hdv, annot, hdv_out=None):
    """_summary_

    Args:
        hdv (_type_): _description_
        annot (_type_): _description_
        hdv_out (_type_, optional): _description_. Defaults to None.
    """
    annot = pd.read_csv(annot, sep='\\t', engine='python')
    hdv = pd.read_csv(hdv, sep='\\t', engine='python')
    
    annot['#CHRM:POS'] = annot['#CHRM:POS'].astype(str)
    annot['CHROM'] = annot['CHROM'].astype(str)
    annot['POS'] = annot['POS'].astype(str)
    annot['REF'] = annot['REF'].astype(str)
    annot['ALT'] = annot['ALT'].astype(str)
    
    hdv['#CHRM:POS'] = hdv['#CHRM:POS'].astype(str)
    hdv['CHROM'] = hdv['CHROM'].astype(str)
    hdv['POS'] = hdv['POS'].astype(str)
    hdv['REF'] = hdv['REF'].astype(str)
    hdv['ALT'] = hdv['ALT'].astype(str)

    hdv_annot = pd.merge(hdv, annot, on=[
                         '#CHRM:POS', 'CHROM', 'POS', 'REF', 'ALT'], how='left').drop_duplicates(subset=['CHROM', 'POS', 'REF', 'ALT'], keep='last').reset_index(drop=True).sort_values(by=['CHROM', 'POS'], ascending=True)
    hdv_annot.to_csv(f'{hdv_out}', index=False, header=True, sep='\\t')


if __name__ == '__main__':
    annot_hdv(args.hdv, args.annot, args.hdv_out)
