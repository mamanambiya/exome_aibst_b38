#!/usr/bin/env python3

import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--fisher", help="", default="${fisher}")
parser.add_argument("--out", help="", default="${out}")
args = parser.parse_args()

def filter_fisher_test(fisher, out):
    """_summary_
    """
    data = pd.read_csv(fisher, delimiter='\\s+')
    data = data.astype({'CHR': 'Int64', 'BP': 'Int64'})
    ## filter for pv < 0.05
    data = data[(data['P'] <= 0.05)
                ].drop_duplicates().reset_index(drop=True)
    data.to_csv(out, index=False, header=True, sep='\\t')


if __name__ == '__main__':
    filter_fisher_test(args.fisher, args.out)
