#!/usr/bin/env python3.7
'''

'''
import argparse
import sys
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument("--counts", default="${counts}", help="")
parser.add_argument("--dataset", default="${dataset}", help="output file")
parser.add_argument("--out", default="${count}", help="output file")
parser.add_argument("--header", default="${header}", help="output file")
args = parser.parse_args()


def combine_counts(counts, out, dataset, header='true'):
    """
    :param counts:
    :param out:
    :return:
    """
    counts = counts.split(' ')
    
    if header == 'true':
        datas = pd.read_csv(counts[0], delimiter=';', quotechar='\\"', skipinitialspace = True,
                            engine='python', header=0, )
        datas.columns = datas.columns.str.replace(' ', '')
    else:
        datas = pd.read_csv(counts[0], delimiter=';', quotechar='\\"', skipinitialspace = True,
                            engine='python')
        datas.columns = datas.columns.str.replace(' ', '')

    for count in counts[1:]:
        if header == 'true':
            dat = pd.read_csv(count, delimiter=';', quotechar='\\"', skipinitialspace = True,
                              engine='python', header=0, )
            dat.columns = dat.columns.str.replace(' ', '')
        else:
            dat = pd.read_csv(count, delimiter=';', quotechar='\\"', skipinitialspace = True,
                              engine='python')
            dat.columns = dat.columns.str.replace(' ', '')
        datas = pd.merge(datas, dat, how='inner', on=['Group'])
        
    datas1 = datas[["Group"]]
    datas1[f'Min_{dataset}'] = datas.min(axis=1).round()
    datas1[f'Max_{dataset}'] = datas.max(axis=1).round()
    datas1[f'Mean_{dataset}'] = datas.mean(axis=1).round()
    datas1[f'Total_{dataset}'] = datas.sum(axis=1)
    datas2 = datas[["Group"]]
    datas2[f'Total_{dataset}'] = datas.sum(axis=1)
    datas.to_csv(out+".all.csv", sep=';', index=False)
    datas1.to_csv(out+".summary.csv", sep=';', index=False)
    datas2.to_csv(out+".total.csv", sep=';', index=False)

if __name__ == '__main__':
    combine_counts(args.counts, args.out, args.dataset, args.header)
