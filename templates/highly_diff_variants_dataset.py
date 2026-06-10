#!/usr/bin/env python3.7
'''

'''
import argparse
import sys
import time
import gzip
import csv
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--annot", default="${annot}", help="TSV annotation file")
parser.add_argument("--pvalues", default="${pvalues}", help="P-value file")
parser.add_argument("--fst_file", default="${fst_file}", help="FST file")
parser.add_argument("--fst_cutoff", default="${fst_cutoff}", help="FST cutoff")
parser.add_argument(
    "--step", default="${fold}", help="number of folds to use")
parser.add_argument(
    "--ac", default="${ac}", help="")
parser.add_argument(
    "--group_treshold", default="${group_treshold}", help="number of populations to use")
parser.add_argument(
    "--test_pops", default="${test_pops}", help="")
parser.add_argument(
    "--base_pops", default="${base_pops}", help="")
parser.add_argument(
    "--outFile", default="${outFile}", help="Annotation output file")
args = parser.parse_args()

# Get pvalues fron fisher test statistic from Plink


def get_pvalues(fisher):
    """_summary_
    """
    pvalues = pd.read_csv(fisher, delimiter='\\s+')
    pvalues = pvalues.astype({'CHR': 'str', 'BP': 'Int64'})
    # Normalize chr prefix: ensure CHR has 'chr' prefix to match FST/freq files
    pvalues['CHR'] = pvalues['CHR'].astype(str).apply(
        lambda x: x if x.startswith('chr') else f'chr{x}')
    # filter for pv < 0.05
    pvalues = pvalues[(pvalues['P'] <= 0.05)
                      ].drop_duplicates().reset_index(drop=True)
    pvs = pvalues.T.to_dict()
    pvalue_snps = {}
    for idx in pvs:
        pv = pvs[idx]
        id = f'{pv["CHR"]}:{pv["BP"]}'
        if id not in pvalue_snps:
            pvalue_snps[id] = {}
        # Handle both multi-pop (POPS column) and single-pop Fisher files
        pop_key = pv.get('POPS', 'ALL')
        pvalue_snps[id][pop_key] = pv['P']

    return pvalue_snps


def get_fst(fst_file, fst_cutoff='', include_snps=''):
    """_summary_
    """
    # data = pd.read_csv(fst_file, delimiter='\\\\s+')
    fst_snps = {}

    for fst in pd.read_csv(fst_file, sep='\\s+', na_filter=True, iterator=True, chunksize=1000000, names=['POPS', 'CHROM', 'POS', 'FST']):
        fst = fst.dropna(subset=['FST'])
        fst['FST'] = fst['FST'].astype(str)
        fst = fst[fst['FST'].str.contains('FST', regex=True) == False]
        fst['FST'] = fst['FST'].astype(float)
        if fst_cutoff != '':
            fst_cutoff = float(fst_cutoff)
            fst = fst[(fst['FST'] > fst_cutoff)]

        fst = fst.T.to_dict()

        for idx in fst:
            data = fst[idx]
            id = f'{data["CHROM"]}:{data["POS"]}'
            if include_snps != '' and len(include_snps) > 0:
                if id in include_snps:
                    if id not in fst_snps:
                        fst_snps[id] = {}
                    fst_snps[id][data['POPS']] = data['FST']
            else:
                if id not in fst_snps:
                    fst_snps[id] = {}
                fst_snps[id][data['POPS']] = data['FST']

    return fst_snps


# Get base dataset name from list
def get_base_dataset_name(header):
    """_summary_
    """
    base_dataset_names = []
    # Get base dataset and AF
    for data in header:
        if data.endswith('_AF_') or data.endswith('_MAF_'):
            if data.endswith('_AF_'):
                dataset = data.split('_AF_')[0]
            elif data.endswith('_MAF_'):
                dataset = data.split('_MAF_')[0]
            # af = float(datas[idx][f'{dataset}_AF_'])
            # maf = float(datas[idx][f'{dataset}_MAF_'])
            if dataset not in base_dataset_names:
                base_dataset_names.append(dataset)
    if (len(base_dataset_names)) == 1:
        dataset = base_dataset_names[0]
    else:
        print(f'Multiple base datasets')
        sys.exit(1)

    return dataset


# TODO add dataset or pop from the process
def highly_diff(frq_file, pvalue_file, fst_file, outFile, test_pops='', base_pops='', fold=3, ac_cutoff=2, group_treshold=2, fst_cutoff='0.25'):
    """
    :param frq_file header #CHRM:POS CHROM   POS     GENE    EFFECT  AC      MAF POP_AF ...:
    :return:
    """

    print(f'Reading frequency file {frq_file} ...')
    datas = pd.read_csv(frq_file, delimiter='\\s+',
                        quotechar='\\"', engine='python')
    header = datas.iloc[0, ].to_dict().keys()
    datas = list(datas.T.to_dict().values())

    ac_cutoff = int(ac_cutoff)
    test_pops = test_pops.strip().split(',')
    base_pops = base_pops.strip().split(',')
    # group_treshold = len(base_pops)*float(group_treshold)
    group_treshold = int(group_treshold)

    freqs = {}
    freqs_pops = {}
    records = {}
    records_zeros = {}
    records_non_zeros = {}
    records_all = {}
    hdv_non_zeros = {}
    records_zeros_list = {}

    af = 0
    maf = 0

    # Get pvalues
    pvalues = get_pvalues(pvalue_file)  # TODO add dataset to filter

    # Get fst
    # TODO add dataset to filter
    fsts = get_fst(fst_file, fst_cutoff, pvalues)

    base = get_base_dataset_name(header)  # TODO add dataset and compare

    if base in base_pops or base in test_pops:
        for record in datas:
            af_zeros = 0
            # Get base allele frequencies
            base_af = float(record[f'{base}_MAF_'])
            base_ac = float(record[f'{base}_AC_'])
            if base_af > 0 and base_af <= 0.5 and base_ac > ac_cutoff and base in base_pops:
                id = f'{record["CHROM"]}:{record["POS"]}'
                records_zeros[id] = 0
                records_non_zeros[id] = 0
                records_zeros_list[id] = []
                if id not in records_all:
                    records_all[id] = 1
                freqs[id] = []
                freqs_pops[id] = []
                records[id] = record

                af_fold_up = base_af/int(fold)
                af_fold_down = base_af*int(fold)

                # Check if variant has FST scores
                if id in fsts:
                    # Get node allele frequencies
                    for data in record:
                        if data.endswith('_AF') and data[:-3] in base_pops and f'{base}_AF_' not in data:
                            node = data[:-3]
                            pops = f'{base}_{node}'
                            pops1 = f'{node}_{base}'
                            node_af = float(record[f'{node}_AF'])

                            if node_af != 0.0:
                                records_all[id] += 1
                            if float(node_af) == float(0):
                                records_zeros[id] += 1
                                records_zeros_list[id].append(node)
                            if pops in fsts[id] or pops1 in fsts[id]:
                                if node_af > 0 and node_af < 0.5:
                                    if node in base_pops:
                                        if af_fold_up >= node_af or af_fold_down <= node_af:
                                            freqs[id].append(float(node_af))
                                            freqs_pops[id].append(node)
                                            records_non_zeros[id] += 1

    hdvs = {}
    for id in freqs_pops:
        if (len(freqs_pops[id]) >= group_treshold) or (records_zeros[id] > (group_treshold - 1) and records_non_zeros[id] >= (group_treshold - 1)):
            if id not in hdvs:
                hdvs[id] = {}
            hdvs[id]['#CHRM:POS'] = records[id]['#CHRM:POS']
            hdvs[id]['CHROM'] = records[id]['CHROM']
            hdvs[id]['POS'] = records[id]['POS']
            hdvs[id]['REF'] = records[id]['REF']
            hdvs[id]['ALT'] = records[id]['ALT']
    hdvs = pd.DataFrame(hdvs)
    if len(hdvs) > 0:
        hdvs = hdvs.T.drop_duplicates(subset=[
            'CHROM', 'POS', 'REF', 'ALT'], keep='last').sort_values(['CHROM', 'POS']).reset_index(drop=True)
    hdvs.to_csv(f'{outFile}_base.tsv', index=False, header=True, sep='\\t')
    hdvs_all = {}
    for id in freqs:
        if (len(freqs_pops[id]) >= group_treshold) or (records_zeros[id] > (group_treshold - 1) and records_non_zeros[id] >= (group_treshold - 1)):
            if id not in hdvs_all:
                hdvs_all[id] = {}
            hdvs_all[id]['#CHRM:POS'] = records[id]['#CHRM:POS']
            hdvs_all[id]['CHROM'] = records[id]['CHROM']
            hdvs_all[id]['POS'] = records[id]['POS']
            hdvs_all[id]['REF'] = records[id]['REF']
            hdvs_all[id]['ALT'] = records[id]['ALT']
    hdvs_all = pd.DataFrame(hdvs_all)
    if len(hdvs_all) > 0:
        hdvs_all = hdvs_all.T.drop_duplicates(subset=['CHROM', 'POS', 'REF', 'ALT'], keep='last').sort_values([
            'CHROM', 'POS']).reset_index(drop=True)
    hdvs_all.to_csv(f'{outFile}_all.tsv', index=False, header=True, sep='\\t')


if args.annot:
    highly_diff(args.annot, args.pvalues, args.fst_file,
                args.outFile, args.test_pops, args.base_pops, args.step, args.ac, args.group_treshold, args.fst_cutoff)
