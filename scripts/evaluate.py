# many to many tests of relationships

import os
import sys
import itertools
import networkx as nx
import pandas as pd
from utils import read_ersa, read_king
#import seaborn as sns
#from matplotlib import pyplot as plt
import gzip

from utils import get_kinship
from utils import read_pedigree
from utils import read_pipeline_output
from utils import line_generator


def read_segments_graph(matchfiles):
    i1, i2, l = 1, 3, 10
    segments_graph = nx.Graph()
    for line in line_generator(matchfiles):
        items = line.split()
        id1 = items[i1]
        id2 = items[i2]
        length = float(items[l])

        if segments_graph.has_edge(id1, id2):
            segments_graph[id1][id2]['length'] += length
            segments_graph[id1][id2]['N'] += 1
        else:
            segments_graph.add_edge(id1, id2, N=1, length=length)

    return segments_graph


def ersa_king(rel_ersa, rel_king, segments_graph, outname):
    """Note FID2 and IID2 are subjects in the background data
    If all=True, give all pairwise relatives between clients,
    instead of that between client and background"""
    for i in segments_graph.nodes:
        if i in rel_king:
            for n in rel_king.neighbors(i):
                # not consider relationships in the group
                degree = rel_king[i][n]['degree']

                if rel_ersa.has_edge(i, n):
                    rel_ersa[i][n]["King_est"] = degree
                else:
                    rel_ersa.add_edge(i, n, King_est=degree)

    w = open(outname, 'w')
    ids_no_relatives = []
    w.write(
        "Target_FID1\tTarget_IID1\tBackground_FID2\tBackground_IID2\tRel_est1\tRel_est2\tDegree\tKing_est\tTotal_seg_num\tTotal_seg_len\n")
    for sample_id in segments_graph:
        fid1, iid1 = sample_id.split('_')
        if sample_id in rel_ersa:
            for n in rel_ersa.neighbors(sample_id):
                fid2, iid2 = n.split('_')
                segs = segments_graph.edges.get((sample_id, n), {})
                # for duplicate values both in target and background
                fid2 = fid2.replace('2:', '')
                dc = rel_ersa[sample_id][n]
                dst = dc.get('d_est', 'NA')
                dst_k = dc.get('King_est', "NA")
                if dst != 'NA':
                    if dst_k == 'PO':
                        dst = 1
                    elif dst_k == 'FS':
                        dst = 2
                    elif dst_k == '2nd':
                        # KING '2nd': half-sibs, avuncular pairs and grandparentgrandchild pairs
                        # so here, dst can be either 2 or 3 by our definition
                        dst = 3
                    else:
                        # KING will override ERSA for degree 1,2 and 3
                        dst = int(dst) + 1
                w.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}\n".format(
                    fid1, iid1, fid2, iid2, dc.get('Rel_est1', "NA"),
                    dc.get('Rel_est2', "NA"), str(dst), dst, segs.get('N', 0), segs.get('length', 0)
                ))
        else:
            ids_no_relatives.append((fid1, iid1))

    '''
    for fid, iid in ids_no_relatives:
        # also output clients with no relatives inferred
        w.write("{}\t{}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n".format(fid, iid))
    '''
    w.close()
    return


def compare(total, correct, outname=None, cutoff=1):
    df = pd.DataFrame(columns=['Degree of Relatedness', '% predicted correct +/- 1 degree relationship'])
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
    if not outname:
        plt.show()
    else:
        plt.savefig(outname + '.png')
        plt.close()


def evaluate(result, fam, output, only_client=False):
    """If only_client=True, all pairwise relatives between client are evaluated"""
    iids, pedigree = read_pedigree(fn=fam)
    _, kinship = get_kinship(pedigree)
    inferred, clients = read_pipeline_output(result, only_client)
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
                # print(i, j, degree, ersa_d, inferred[i][j]['king'])
                if (degree - 1 <= ersa_d <= degree + 1):
                    correct[degree] = correct.get(degree, 0) + 1
            else:
                # print(i, j, degree, "NA", "NA")
                pass
    if not total:
        return
    compare(total, correct, output)


if __name__ == '__main__':
    rel_king = read_king(snakemake.input['king'][0])
    segments_graph = read_segments_graph(snakemake.input['germline'])
    rel_ersa = read_ersa(snakemake.input['ersa'][0])
    ersa_king(rel_ersa, rel_king, segments_graph, snakemake.output[0])
