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
parser.add_argument(
    "--outFile", default="${outFile}", help="Annotation output file")
args = parser.parse_args()


def highly_diff(frq_file, outFile):
    """
    :param frq_file header #CHRM:POS CHROM   POS     GENE    EFFECT  AC      MAF POP_AF ...:
    :return:
    """
    fold = 2
    results_up = {}
    results_down = {}
    df = pd.read_csv(frq_file, delimiter='\\\\s+',
                        quotechar='\\"', engine='python')

    # Handle empty input files (header only, no data rows)
    if df.empty:
        print(f'WARNING: Frequency file {frq_file} is empty (no data rows). Generating empty output files.')
        # Create empty output files
        empty_df = pd.DataFrame()
        empty_df.to_csv(f'{outFile}.sites.tsv', index=False, header=True, sep='\\t')
        empty_df.to_csv(f'{outFile}.all.tsv', index=False, header=True, sep='\\t')
        empty_df.to_csv(f'{outFile}.genes.tsv', index=False, header=True, sep='\\t')
        empty_df.to_csv(outFile, index=False, header=True, sep='\\t')
        return

    datas = df.T.to_dict()
    header = datas[0].keys()
    ids = datas.keys()
    results_up_top = []
    snp_freqs = {}

    for idx in ids:
        freqs = []
        freqs1 = []
        dataset = ''
        for data in header:
            if '_MAF_' in data:
                maf = float(datas[idx][data])
                maf_fold_up = maf/fold
                maf_fold_down = maf*fold
                dataset = data.split('_MAF_')[0]
                datas[idx][data] = round(float(maf), 2)
            if '_AF' in data or '_MAF' in data or 'TOPMED' in data:
                allele_frq = datas[idx][data]
                datas[idx][data] = round(float(allele_frq), 2)
            if ('_AF' in data or '_MAF' in data or 'TOPMED' in data) and '_MAF_' not in data and dataset+"_MAF" not in data:
                allele_frq = datas[idx][data]
                if data in ['AGVP_AF', 'KG_AF', 'KG_AFR_AF', 'gnomAD_AF', 'gnomAD_AFR_AF', 'TOPMED', 'GBR_MAF', 'CEU_MAF']:
                    freqs.append(float(allele_frq))
                    if len(data.split('_')) <= 2:
                        freqs1.append(float(allele_frq))

        # Upstream
        if len(freqs) > 0 and maf_fold_up > 0 and not all(frq == 0 for frq in freqs) and len([resp for resp in [maf_fold_up >= frq for frq in freqs if frq > 0.0] if resp]) >= 3:
            effect = datas[idx]['EFFECT']
            try:
                effect = (effect.split('&')[0]).split('_gene_variant')
            except:
                print(effect)
            datas[idx]['EFFECT'] = ' '.join(effect)
            diffs = []
            n = 0
            for frq in freqs1:
                if frq != 0 and maf_fold_up >= frq:
                    n += 1
            if n >= int(len(freqs)/2):
                results_up_top.append(idx)
            for freq in freqs:
                try:
                    diffs.append(float(maf)/float(freq))
                except:
                    if int(freq) == 0:
                        continue
            diff = int(sum(diffs)/len(diffs))
            group = diff - (diff % 10) + 10
            if not group in results_up:
                results_up[group] = []
            if idx in results_up[group]:
                results_up[group].remove(idx)
            results_up[group].insert(0, idx)
            snp = datas[idx]['#CHRM:POS']
            if len([resp for resp in [maf_fold_up >= frq for frq in freqs if frq > 0.0] if resp]) >= 2:
                snp_freqs[snp] = freqs

        # Downsteam
        if len(freqs) > 0 and maf_fold_down > 0 and not all(frq == 0 for frq in freqs) and len([resp for resp in [maf_fold_down <= frq for frq in freqs if frq > 0.0] if resp]) >= 3:
            # TODO in at least 2
            effect = datas[idx]['EFFECT']
            effect = (effect.split('&')[0]).split('_gene_variant')
            datas[idx]['EFFECT'] = ' '.join(effect)
            diffs = []
            n = 0
            for frq in freqs1:
                if frq != 0 and maf_fold_up <= frq:
                    n += 1
            if n >= int(len(freqs)/2):
                results_up_top.append(idx)
            for freq in freqs:
                try:
                    diffs.append(float(maf)/float(freq))
                except:
                    if int(freq) == 0:
                        continue
            diff = int(sum(diffs)/len(diffs))
            group = diff - (diff % 10) + 10
            if not group in results_down:
                results_down[group] = []
            if idx in results_down[group]:
                results_down[group].remove(idx)
            results_down[group].insert(0, idx)
            snp = datas[idx]['#CHRM:POS']
            if len([resp for resp in [maf_fold_down*2 <= frq for frq in freqs if frq > 0.0] if resp]) >= 2:
                snp_freqs[snp] = freqs

    out = open(outFile, 'w')  # Top 20
    out_CLNDN = open(outFile+'.CLNDN.tsv', 'w')  # All
    out_all = open(outFile+'.all.tsv', 'w')  # All
    sites = open(outFile+'.sites.tsv', 'w')  # Sites
    sites.writelines('\\t'.join(['#CHRM:POS']) + '\\n')
    out.writelines('\\t'.join(header) + '\\n')
    out_all.writelines('\\t'.join(header) + '\\n')
    genes = {}
    results = {}
    n = 0

    # Upstream
    n = 0
    for group in sorted(results_up, reverse=True):
        for idx in results_up[group]:
            data = []
            for col in header:
                data.append(str(datas[idx][col]))
                if '_AC_' in col:
                    ac = float(datas[idx][col])
            if ';' in data[0]:
                snps = data[0].split(';')
                for snp in snps:
                    if ac >= 2:
                        out_all.writelines('\\t'.join(
                            [snp] + data[1:]) + '\\n')
            else:
                if ac >= 3:
                    out_all.writelines('\\t'.join(data) + '\\n')
            snp = datas[idx]['#CHRM:POS']
            gene = datas[idx]['GENE']

            CLNDN = datas[idx]['CLNDN']
            if CLNDN != '' and CLNDN != '.':
                if ';' in data[0]:
                    snps = data[0].split(';')
                    for snp in snps:
                        out_CLNDN.writelines(
                            '\\t'.join([snp] + data[1:])+'\\n')
                else:
                    out_CLNDN.writelines('\\t'.join(data)+'\\n')
                if idx in results_up_top and ac >= 2:
                    if ';' in data[0]:
                        snps = data[0].split(';')
                        for snp in snps:
                            results[snp] = [snp] + data[1:]
                    else:
                        results[snp] = [snp] + data[1:]

                if gene not in genes:
                    genes[gene] = []
                if idx in results_up_top and ac >= 2 and maf > 0:
                    genes[gene].append(snp)

                sites.writelines('\\t'.join([snp])+'\\n')

    # Downstream
    n = 0
    for group in sorted(results_down, reverse=True):
        for idx in results_down[group]:
            data = []
            for col in header:
                data.append(str(datas[idx][col]))
                if '_AC_' in col:
                    ac = float(datas[idx][col])
            if ';' in data[0]:
                snps = data[0].split(';')
                for snp in snps:
                    if ac >= 3:
                        out_all.writelines('\\t'.join(
                            [snp] + data[1:]) + '\\n')
            else:
                if ac >= 2:
                    out_all.writelines('\\t'.join(data) + '\\n')
            snp = datas[idx]['#CHRM:POS']
            gene = datas[idx]['GENE']

            CLNDN = datas[idx]['CLNDN']
            if CLNDN != '' and CLNDN != '.':
                if ';' in data[0]:
                    snps = data[0].split(';')
                    for snp in snps:
                        out_CLNDN.writelines(
                            '\\t'.join([snp] + data[1:])+'\\n')
                else:
                    out_CLNDN.writelines('\\t'.join(data)+'\\n')
                if idx in results_up_top and ac >= 2:
                    if ';' in data[0]:
                        snps = data[0].split(';')
                        for snp in snps:
                            results[snp] = [snp] + data[1:]
                    else:
                        results[snp] = [snp] + data[1:]
                n += 1

                if gene not in genes:
                    genes[gene] = []
                if idx in results_up_top and ac >= 2 and maf > 0:
                    genes[gene].append(snp)

                sites.writelines('\\t'.join([snp])+'\\n')

    out_genes = open(outFile + '.genes.tsv', 'w')
    out_genes.writelines('\\t'.join(['GENE', 'SNPS'])+'\\n')
    for gene in genes:
        out_genes.writelines('\\t'.join([gene, str(len(genes[gene]))])+'\\n')

    g = sorted(genes, key=lambda k: len(genes[k]), reverse=True)
    for gene in g:
        if len(genes[gene]) > 0:
            #         if gene in [ 'ABCB1','ABCC2','CYP2A6','CYP2B6','CYP2D6','CYP2C19','CYP2C8','CYP3A5','CYP3A4','NAT2','SLCO1B1','UGT2B7','UGT1A1','NF-κB1','TNF-α','NR1I2','AGBL4','NR1I3','ABCC4','SLC28A2']:
            #             print(g.index(gene), gene, len(genes[gene]))
            if gene in g[:20]:
                n = 0
                snps = list(genes[gene])
                for snp1 in snps:
                    if snp1 in results and n == 0:
                        if snp1 in snp_freqs:
                            if snp1.startswith('rs') and (float(results[snp1][6]) >= 0.1 or float(results[snp1][9])-float(results[snp1][6]) >= 0.05):
                                results[snp1][3] = f"{gene} ({len(genes[gene])})"
                                n += 1
                                out.writelines('\\t'.join(results[snp1])+'\\n')
                    else:
                        continue

    out.close()
    sites.close()
    out_all.close()
    out_genes.close()

if args.annot:
    highly_diff(args.annot, args.outFile)
