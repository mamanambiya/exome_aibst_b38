#!/usr/bin/env python3

import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--annot", default="${annot}", help="")
parser.add_argument("--outAnnot", default="${outAnnot}", help="")
parser.add_argument("--hdr", default="${hdr}", help="")
parser.add_argument("--outHdr", default="${outHdr}", help="")
args = parser.parse_args()


def fix_mafs(annot, outAnnot, hdr, outHdr):
    '''
    :param fts_2by2_input:
    :return fst_matrix_output:
    '''
    data = pd.read_csv(annot, delimiter='\\\\s+', quotechar='\\"',
                       engine='python')
    data.columns = data.columns.str.replace('_x', '')
    data.columns = data.columns.str.replace('_y', '')
    data.columns = data.columns.str.replace('POS.1', '')
    data.to_csv(outAnnot, sep='\\t', index=False)

    ## Headers
    hd = set([i.replace('_x', '').replace('_y', '')
                for i in open(hdr).readlines()])
    outHdr = open(outHdr, 'w')
    outHdr.writelines(''.join(hd))
    outHdr.close()


fix_mafs(args.annot, args.outAnnot, args.hdr, args.outHdr)
