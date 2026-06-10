#!/usr/bin/env python3.7

import argparse
import sys
import time as t
import numpy as np
from pandas import DataFrame
import seaborn as sns
# %matplotlib inline

parser = argparse.ArgumentParser()
parser.add_argument("--fst_matrix", help="", default="${fst_matrix}")
parser.add_argument("--fst_matrix_plot", help="", default="${fst_matrix_plot}")
args = parser.parse_args()


def plot_fst_matrix(fst_matrix, fst_matrix_plot):
    """_summary_

    Args:
        fst_matrix (_type_): _description_
        fst_matrix_plot (_type_): _description_
    """

    # Generate a mask for the upper triangle
    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(230, 50, as_cmap=True)
    mask = np.triu(np.ones_like(fst_matrix, dtype=bool))
    plt.figure(figsize=(20, 15))
    sns.heatmap(fst_matrix, cmap=cmap, annot=True, center=0, vmin=0.03,
                square=True, linewidths=.5, cbar_kws={"shrink": .3}, mask=mask)
    plt.savefig(fst_matrix_plot)


if __name__ == "__main__":
    plot_fst_matrix(args.fst_matrix, args.fst_matrix_plot)
