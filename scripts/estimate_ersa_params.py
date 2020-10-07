import networkx as nx
from utils.relatives import line_generator, read_king


def write_l_th(l, th, t, file_path):
    with open(file_path, 'w') as file:
        file.write(f'ERSA_L={l:.3f}\n')
        file.write(f'ERSA_TH={th:.3f}\n')
        file.write(f'ERSA_T={t:.3f}')


def get_l_th(matchfiles, rel_king, ERSA_t=2.5):
    """
        Infer l and th from pairs not in King for ERSA run
        l is the average number of IBD segments in population
        th is the average length of IBD segment
        ERSA_t is min length of segment to be considered in segment aggregation
    """

    # Indices of first sample id, second sample id and segment length columns
    i1, i2, l = 1, 3, 10
    # Graph with edges not in close kinship as detected by King
    unrelated = nx.Graph()

    for line in line_generator(matchfiles):
        items = line.split()
        id1 = items[i1]
        id2 = items[i2]
        length = float(items[l])
        if (items[11] == 'cM') and (not rel_king.has_edge(id1, id2)) and (length > ERSA_t) and (
                length < 10):
            # < 10 is only for background
            # TODO: why magic < 10 number
            if unrelated.has_edge(id1, id2):
                unrelated[id1][id2]['length'] += length
                unrelated[id1][id2]['N'] += 1
            else:
                unrelated.add_edge(id1, id2, N=1, length=length)

    total_len, total_N = 0, 0
    for i, j in unrelated.edges:
        total_len += unrelated[i][j]['length']
        total_N += unrelated[i][j]['N']

    th = total_len / total_N
    l = total_N / len(unrelated.edges)
    return l, th


if __name__ == '__main__':
    rel_king = read_king(snakemake.input['king'][0])
    #l, th = get_l_th(snakemake.input['germline'], rel_king)

    # best params are
    #l = 1.33
    #th = 3.5

    l = 1.33
    th = 4.0
    write_l_th(l, th, 2.5, snakemake.output[0])
