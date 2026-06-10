#!/usr/bin/env python3.7
'''

'''
import argparse, sys, time
import pandas as pd
import numpy as np
import csv

parser = argparse.ArgumentParser()
parser.add_argument("--maf", default="", help="One or many TSV annotation files, if many use ';' as separator")
parser.add_argument("--out", default="", help="TSV annotation file")
args = parser.parse_args()

def split_maf_to_chromosome(maf, out):
    """[summary]
    
    Arguments:
        maf {[type]} -- [description]
        out {[type]} -- [description]
    """
    outs = {}
    chrms = {}
    datas = {}
    mafs = pd.DataFrame.to_dict(pd.read_csv(maf, sep='\\s+', engine='python'))
    for col in mafs:
        if col == '#CHRM:POS':
            for idx in mafs[col]:
                if idx not in datas:
                    datas[idx] = []
                chrm = mafs['#CHRM:POS'][idx].strip().split(':')[0]
                pos = mafs['#CHRM:POS'][idx].strip().split(':')[1]
                chrm_pos = mafs['#CHRM:POS'][idx]
                if chrm not in outs:
                    outs[chrm] = open("{}.chr{}.mafs".format(out, chrm), 'w')
                    outs[chrm].writelines('\t'.join(mafs.keys()) + '\n')
                chrms[idx] = chrm
                datas[idx].append(mafs[col][idx])
        else:
            for idx in mafs[col]:
                datas[idx].append(mafs[col][idx])
    print("Writing ouput files ...")
    for idx in datas:
        outs[chrms[idx]].writelines('\t'.join([str(it) for it in datas[idx]])+'\n')


split_maf_to_chromosome(args.maf, args.out)
