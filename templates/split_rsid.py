#!/usr/bin/env python3.7
'''

'''
import argparse
import sys
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument("--frq_file", default="${frq_file}", help="")
parser.add_argument("--out", default="${out}", help="output file")
args = parser.parse_args()


def split_rsid(frq, out):
    """
    :param counts:
    :param out:
    :return:
    """
    freq = pd.read_csv(frq, delimiter='\\t')

    # Handle empty input files (header only, no data rows)
    if freq.empty:
        # Write empty output with header
        freq.to_csv(out, index=False, header=True, sep='\\t')
        return

    freqs = freq.T.to_dict()

    header = list(freq.columns)
    new = []
    for key in freqs:
        for col in header:
            if '_MAF' in col or '_AF' in col:
                maf = freqs[key][col]
                freqs[key][col] = round(float(maf), 5)
        rsID = str(freqs[key]['#CHRM:POS'])
        try:
            if ';' in rsID:
                rsID = rsID.split(';')
                for rs in rsID:
                    if 'rs' in rs:
                        freqs[key]['#CHRM:POS'] = rs
                        new.append(list(freqs[key].values()))
            else:
                new.append(list(freqs[key].values()))
        except:
            print(list(freqs[key].values()))

    new_ = pd.DataFrame(new, columns=header)
    new_.to_csv(out, index=False, header=True, sep='\\t')


if __name__ == '__main__':
    split_rsid(args.frq_file, args.out)
