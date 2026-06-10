#!/usr/bin/env python3.7

import argparse
'''

'''

parser = argparse.ArgumentParser()
parser.add_argument("--sites", default="${sites}", help="TSV annotation file")
parser.add_argument(
    "--csv_out", default="${csv_out}", help="Annotation output file")
args = parser.parse_args()


def make_circos(sites, csv_out):
    '''
    '''
    lists = sites.split(';')
    datas = {}
    for fil in lists:
        fil = fil.split('==')
        dataset = fil[0]
        site = fil[1]
        if dataset not in datas:
            datas[dataset] = []
        datas[dataset] += [it.strip() for it in open(site).readlines() if '#CHRM:POS' not in it]
    datas1 = []
    store = []
    my_datas = sorted(datas, reverse=True)
    for dataset in my_datas:
        for it in datas[dataset]:
            for dataset1 in my_datas:
                # if dataset != dataset1:
                if it in datas[dataset1]:
                    if dataset+'_'+dataset1+'_'+it not in store and dataset1+'_'+dataset+'_'+it not in store:
                        datas1.append([dataset, dataset1, it])
                        store.append(dataset+'_'+dataset1+'_'+it)
                        store.append(dataset1+'_'+dataset+'_'+it)
    out = open(csv_out, 'w')
    out.writelines(' '.join(['from', 'to', 'value'])+'\\n')
    datas1 = sorted(datas1, key=(lambda x: (x[0], x[1], x[2])), reverse=True)
    for data in datas1:
        out.writelines(' '.join(data)+'\\n')
    out.close()


if __name__ == '__main__':
    make_circos(args.sites, args.csv_out)
