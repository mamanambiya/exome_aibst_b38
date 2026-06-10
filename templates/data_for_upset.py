#!/usr/bin/env python3.7
'''

'''
import argparse
import pandas as pd
import re
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument(
    "--sites", default="${sites}", help="TSV annotation file")
parser.add_argument(
    "--csv_out", default="${csv_out}", help="Annotation output file")
args = parser.parse_args()


def data_for_upset(sites, outFile):
    """
    :param frq_file:
    :return:
    """

    sites = [it.strip() for it in sites.split(',')]
    
    datas = pd.read_csv(sites[0], delimiter='\\\\s+',
                        quotechar='\\"', engine='python')
    for site in sites[1:]:
        b = pd.read_csv(site, delimiter='\\\\s+', quotechar='\\"', engine='python')
        datas = pd.merge(datas, b, how='outer', on='#CHRM:POS')
    datas = datas.replace(np.nan, '0', regex=True)
    # datas.drop(['#CHRM:POS'], axis='columns', inplace=True)
    datas.to_csv(outFile, sep=',', index=False)


if __name__ == '__main__':
    data_for_upset(args.sites, args.csv_out)
