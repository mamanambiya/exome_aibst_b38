#!/usr/bin/env python3

import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--fst_weir", help="", default="${fst_weir}")
parser.add_argument("--mean_fst_out", help="", default="${mean_fst_out}")
parser.add_argument("--fst_weir_out", help="", default="${fst_weir_out}")
parser.add_argument("--pops", help="", default="${pops}")
args = parser.parse_args()


def process_fst(fst_weir, fst_weir_out, mean_fst_out, pops=''):
    """_summary_
    """
    fst = pd.read_csv(fst_weir, sep='\\s+', na_filter=True)
    fst.rename(columns={'WEIR_AND_COCKERHAM_FST': 'FST'}, inplace=True)
    fst = fst.dropna(subset=['FST'])
    # fst = fst[(fst['FST'] > 0)]
    mean = fst['FST'].mean()
    pops = pops.split('__')
    fst.insert(0, 'POPS', '_'.join(pops))
    # pops = fst['POPS'].to_list()[0].split('_')
    fst.to_csv(f'{fst_weir_out}', index=False, header=True, sep='\\t')
    # mean_fst_header = f'POS\\tPOP\\t\\n'
    dat = f'{pops[0]}\\t{pops[1]}\\t{mean}\\n'
    out = open(f'{mean_fst_out}', 'w')
    out.writelines(dat)
    out.close()


if __name__ == '__main__':
    process_fst(args.fst_weir, args.fst_weir_out, args.mean_fst_out, args.pops)
