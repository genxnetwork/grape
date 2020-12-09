import os
import glob

import shutil
import networkx as nx
import pandas

from itertools import combinations


def clean_dir(directory, recursive=False):
    if recursive:
        if os.path.isdir(directory):
            shutil.rmtree(directory)
    else:
        for i in os.listdir(directory):
            name = os.path.join(directory, i)
            if os.path.isfile(name):
                os.remove(name)


def remove(filename, one=False):
    """Remove all files with filename or specified prefix"""
    if one:
        try:
            os.remove(filename)
        except OSError:
            pass
    else:
        for i in glob.glob("{}.*".format(filename)):
            os.remove(i)
    return


def check_dir(dirname):
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
    return dirname


def check_exist(filename, directory=False):
    """Check is a file or directory exist"""
    if os.path.exists(filename):
        if directory:
            if os.path.isdir(filename):
                return True
        else:
            if os.path.isfile(filename):
                return True
    return False


def read_pedigree(fn):
    """Read .fam into a network object"""
    pedigree = nx.DiGraph()
    iids = set()
    for i in open(fn):
        fid, iid, father, mother = i.split()[:4]
        #iid = iid.split('_')[1]
        full_id = f'{fid}_{iid}' if '_' not in iid else iid
        iids.add(full_id)
        if father != '0':
            #f = f.split('_')[1]
            father_iid = f'{fid}_{father}' if '_' not in father else father
            pedigree.add_edge(father_iid, full_id)
        if mother != '0':
            #m = m.split('_')[1]
            mother_iid = f'{fid}_{mother}' if '_' not in mother else mother
            pedigree.add_edge(mother_iid, full_id)
    return iids, pedigree


def relatives_to_graph(path, only_client=False):
    # id1     id2     total_seg_len   seg_count       king_degree     ersa_degree     final_degree
    g = nx.Graph()

    relatives = pandas.read_table(path)
    clients = set(relatives.id1)
    if only_client:
        clients |= set(relatives.id2)
    edges = [(r['id1'], r['id2'], {'ersa': r['ersa_degree'], 'king': r['king_degree']}) for _, r in relatives.iterrows()]
    g.add_edges_from(edges)
    return g, clients


def read_pipeline_output(fn, only_client=False):
    """Given our pipeline output, return a networkx object"""
    clients = set()
    g = nx.Graph()
    with open(fn) as f:
        next(f)
        for line in f:
            items = line.strip().split(sep="\t")
            clients.add(items[0])
            if only_client:
                clients.add(items[1])
            if items[-1] != 'NA':
                g1, g2 = f'{items[0]}_{items[1]}', f'{items[2]}_{items[3]}'
                g.add_edge(g1, g2, ersa=items[-4], king=items[-3])
    return g, clients


def get_kinship(pedigree: nx.DiGraph) -> nx.Graph:
    """
    Given a networkx object, get all the pairwise relationships

    Args:
        pedigree (networkx.DiGraph): directed pedigree graph with only 1st degree relations
    Returns:
        networkx.Graph: undirected graph with all pairwise relationships
    """
    # get all the descendants (0 length means self)
    v_relatives_ = list(nx.shortest_path_length(pedigree))
    v_relatives = []
    for i, v in v_relatives_:
        temp = {}
        for j in v:
            if v[j] != 0:
                # remove self relationship
                temp[j] = v[j]
        v_relatives.append((i, temp))

    relatives = nx.Graph()
    for ancestor, descendants in v_relatives:
        for i, j in combinations(descendants, 2):
            # len(descendants) = 0/1 will return []
            degree = descendants[i] + descendants[j]
            if relatives.has_edge(i, j):
                if relatives[i][j]['degree'] > degree:
                    relatives.add_edge(i, j, common_ancestor=ancestor, degree=degree)
                elif relatives[i][j]['degree'] == degree:
                    pass
            else:
                relatives.add_edge(i, j, common_ancestor=ancestor, degree=degree)
        # add vertical relationship
        for m in descendants:
            relatives.add_edge(ancestor, m, common_ancestor=ancestor, degree=descendants[m])

    return relatives


def read_across_kin(filepath: str) -> nx.Graph:
    """
    Reads king output and returns networkx.Graph with all relations found by king as nodes in the graph

    Args:
        filepath (str): full path to king output file
    Returns:
        networkx.Graph: relations between samples found by king
    """
    relation = nx.Graph()
    with open(filepath) as file:

        start = next(file, None)
        if start is None:
            return relation

        for line in file:
            items = line.strip().split(sep="\t")
            if items[-1] == 'UN':
                continue
            else:
                name1 = '_'.join(items[:2])
                name2 = '_'.join(items[2:4])
                relation.add_edge(name1, name2, degree=items[-1])
    return relation


def read_king(prefix):
    if os.path.exists(prefix):
        relations = read_across_kin(prefix)
    else:
        relations = nx.Graph()

    return relations


def read_ersa(filepath: str) -> nx.Graph:
    """
    Reads ersa output and returns networkx.Graph with all relations found by ersa as nodes in the graph

    Args:
        filepath (str): full path to ersa output file
    Returns:
        networkx.Graph: relations between samples found by ersa
    """
    relations = nx.Graph()
    with open(filepath, 'r') as file:
        next(file)
        for line in file:
            items = [item.strip() for item in line.split('\t')]
            if items[2] != 'NA' or items[3] != 'NA':
                relations.add_edge(items[0], items[1], Rel_est1=items[2], Rel_est2=items[3], d_est=items[4])

    return relations


def line_generator(matchfiles):
    for i in matchfiles:
        iterator = open(i)
        for line in iterator:
            yield line
