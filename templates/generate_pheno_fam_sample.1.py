#!/usr/bin/env python3

import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--inSample", help="", default="${inSample}")
parser.add_argument("--inFam", help="", default="${inFam}")
parser.add_argument("--pheno", help="", default="${pheno}")
args = parser.parse_args()


def update_sample(inFam, inSample, pheno):
    """
    """
    out = open(pheno, 'w')
    samples = {}
    groups = []
    for line in open(inSample):
        line = line.strip().split()
        pop = line[1]
        sample = line[0]
        if pop not in groups:
            groups.append(pop)
        samples[sample] = str(groups.index(pop)+1)
        # if group in cases:
        #     samples[line[0]] = '2'
        # elif group in controls:
        #     samples[line[0]] = '1'
        # else:
        #     print("Wrong phenotype!!!!")
        #     sys.exit(1)
    for line in open(inFam):
        line = line.strip().split()
        id = line[1]
        out.writelines('\\t'.join(['0', id, samples[id]+"\\n"]))
    out.close()


if __name__ == '__main__':
    update_sample(args.inFam, args.inSample, args.pheno)
