# many to many tests of relationships

import networkx as nx
from utils.relatives import read_ersa, read_king
from utils.relatives import line_generator


def read_segments_graph(matchfile):
    i1, i2, l = 1, 3, 10
    segments_graph = nx.Graph()
    for line in line_generator([matchfile]):
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
    """Note FID2 and IID2 are subjects in the background data"""
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
        # print('sample_id 1 is: ', sample_id)
        fid1, iid1 = sample_id.split('_')
        if sample_id in rel_ersa:
            for i, n in enumerate(rel_ersa.neighbors(sample_id)):
                # print('rel-ersa is ', n)
                fid2, iid2 = n.split('_')
                if fid1 == fid2 and iid1 == iid2:
                    # we do not write relative degree of samples with themselves
                    continue
                segs = segments_graph.edges.get((sample_id, n), {})
                # for duplicate values both in target and background
                fid2 = fid2.replace('2:', '')
                dc = rel_ersa[sample_id][n]
                dst = dc.get('d_est', 'NA')
                dst_k = dc.get('King_est', "NA")
                # print(iid1, iid2, dst, dst_k)
                if dst_k != 'NA':
                    if dst_k == 'PO':
                        dst = 1
                    elif dst_k == 'FS':
                        dst = 2
                    elif dst_k == '2nd':
                        # KING '2nd': half-sibs, avuncular pairs and grandparentgrandchild pairs
                        # so here, dst can be either 2 or 3 by our definition
                        dst = 2
                    elif dst_k == '3rd':
                        dst = 3
                    #elif dst_k == '4th':
                    #    dst = 4
                    else:
                        pass
                        # KING will override ERSA for degree 1,2 and 3
                        #dst = int(dst) + 1
                w.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}\n".format(
                    fid1, iid1, fid2, iid2, dc.get('Rel_est1', "NA"),
                    dc.get('Rel_est2', "NA"), str(dst), dst_k, segs.get('N', 0), segs.get('length', 0)
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


if __name__ == '__main__':
    rel_king = read_king(snakemake.input['king'][0])
    # print(rel_king.nodes)
    # print()
    # print(rel_king.edges)
    degrees = {}
    for (u, v, deg) in rel_king.edges.data('degree', default='NA'):
        degrees[deg] = degrees.get(deg, 0) + 1
    print('king degrees are: ')
    print(degrees)
    segments_graph = read_segments_graph(snakemake.input['germline'])
    rel_ersa = read_ersa(snakemake.input['ersa'][0])
    ersa_king(rel_ersa, rel_king, segments_graph, snakemake.output[0])
    print('merged king and ersa output')

