# many to many tests of relationships

import os
import sys
import itertools
import networkx as nx
import pandas as pd
from utils import read_ersa, read_king
import seaborn as sns
from matplotlib import pyplot as plt

from utils import get_kinship
from utils import read_pedigree
from utils import read_pipeline_output


def compare(total, correct, plot_name=None, cutoff=1):
    ds = []
    cs = []
    for i in sorted(total):
        if i >= cutoff:
            c = correct.get(i, 0)
            cs.append((c / total[i]) * 100)
            ds.append(i)
    df = pd.DataFrame(columns=['Degree of Relatedness', '% predicted correct +/- 1 degree relationship'])
    df.iloc[:, 0] = ds
    df.iloc[:, 1] = cs
    ax = sns.barplot(x="Degree of Relatedness", y="% predicted correct +/- 1 degree relationship", data=df,
                     palette='muted')
    for _, row in df.iterrows():
        ax.text(row.name, row.iloc[1] + 1, round(row.iloc[1], 1), color='black', ha="center")
    ax.set(ylim=(0, 100))
    if not plot_name:
        plt.show()
    else:
        plt.savefig(plot_name)
        plt.close()


def evaluate(result, fam, plot_name, only_client=False):
    """If only_client=True, all pairwise relatives between client are evaluated"""
    iids, pedigree = read_pedigree(fn=fam)
    _, kinship = get_kinship(pedigree)
    inferred, clients = read_pipeline_output(result, only_client)
    #print(len(inferred), inferred.edges)
    #print('*'*100)
    #print(len(kinship), kinship.edges)
    total = {}
    correct = {}
    if only_client:
        iterator = itertools.combinations(clients, 2)
    else:
        iterator = itertools.product(clients, iids - clients)
    for i, j in iterator:
        # check if name is this format 'fid_iid'
        if '_' in i:
            i = "_".join(i.split('_')[1:])
        if '_' in j:
            j = "_".join(j.split('_')[1:])
        if kinship.has_edge(i, j):
            degree = kinship[i][j]['degree']
            total[degree] = total.get(degree, 0) + 1
            if inferred.has_edge(i, j) and inferred[i][j]['ersa'] != 'NA':
                ersa_d = int(inferred[i][j]['ersa'])
                print(i, j, degree, ersa_d, inferred[i][j]['king'])
                if (degree <= ersa_d <= degree + 2):
                    correct[degree] = correct.get(degree, 0) + 1
            else:
                # print(i, j, degree, "NA", "NA")
                pass
    if not total:
        return
    compare(total, correct, plot_name)


if __name__ == '__main__':
    evaluate(snakemake.input['rel'], snakemake.input['fam'], snakemake.output[0], only_client=True)
