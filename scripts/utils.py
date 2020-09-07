import os
import re
import sys
import glob
import argparse

import shutil
import random
import networkx as nx

from itertools import combinations

ATGC = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
atgc = {'A': 1, 'T': 1, 'G': 0, 'C': 0}


def is_ambiguous(a1, a2):
    """if SNP is ambiguous or not"""
    if atgc[a1] == atgc[a2]:
        return True
    return False


def match_ref(a1, a2, ref, alt):
    """Match a1, a2 with ref and alt
    None: match, 1: swap, 2: flip, 3: flip&swap, 4: exclude
    """
    if a1 == ref and a2 == alt:
        return
    elif a1 == alt and a2 == ref:
        return 1
    else:
        ambiguous = is_ambiguous(a1, a2)
        if ambiguous:
            return 4
        else:
            a1_f, a2_f = [ATGC[i] for i in [a1, a2]]
            if a1_f == ref and a2_f == alt:
                return 2
            elif a1_f == alt and a2_f == ref:
                return 3
            else:
                return 4


def check_duplicate_positions(bimfile):
    c = {}
    remove_list = set()
    for i in open(bimfile):
        items = i.split()
        k = (items[0], items[3])
        c[k] = c.get(k, 0) + 1
        if c[k] > 1:
            remove_list.add(k)
    return remove_list


def check_fam(famfile):
    """'_' are reserved for Eagle, space for Germline, ':' duplicate samples; '-' for ersa; '%' for RaPID"""
    ids = []
    for i in open(famfile):
        fid, iid = i.split()[:2]
        if re.search(r'[\s_:%-]', fid):
            raise (Exception("FID has illegal _/:/-/space for {}".format(famfile)))
        if re.search(r'[\s_:%-]', iid):
            raise (Exception("IID has illegal _/:/-/space for {}".format(famfile)))
        if (fid, iid) in ids:
            raise (Exception("Duplicate {} {} in {}".format(fid, iid, famfile)))
        ids.append((fid, iid))
    return ids


def check_rs_id(bimfile):
    """Check if bim file contains SNP rs ids"""
    snps = [i.split()[1] for i in open(bimfile)]
    snps_r = random.sample(snps, 10)
    for i in snps_r:
        # select random snp, to guess if there is rs id
        if re.match(r"^rs\d+$", i):
            return True
    # if not rs id, modified the bim file
    os.rename(bimfile, bimfile + '_')
    w = open(bimfile, 'w')
    for i in open(bimfile + '_'):
        items = i.split()
        items[1] = "{}:{}".format(items[0], items[3])
        w.write("\t".join(items) + "\n")
    w.close()


def get_parser():
    """ Default: perferm a relative matching
    -p: only perferm phasing
    -i: perferm phasing and imputation
    """

    description = "Perform a relative matching between client and background"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-c", "--client", help='Client input VCF file name')
    parser.add_argument("-o", "--output", required=True,
                        help='output directory for phased/imputed/inferred/evaluation file')
    parser.add_argument("-r", "--result", help='Pipeline result file, only needed when perform evaluation')
    parser.add_argument("-e", "--evaluate", default='', help="Given fam file, evaluate the inference result")
    parser.add_argument("-v", "--verbose", action='store_true', help="Whether to save temp data to output folder")

    group1 = parser.add_mutually_exclusive_group()
    group1.add_argument("-p", "--phase", action='store_true', help="Phase the client input")
    group1.add_argument("-i", "--impute", action='store_true', help="Phase and impute the client input")

    group2 = parser.add_mutually_exclusive_group()
    group2.add_argument("-b", "--background", help='Phased background input VCF file name')
    group2.add_argument("-a", "--all", action='store_true',
                        help='Perform a pairwise searching for all samples in clients')
    return parser


def check_params(parser):
    """Check parser parameters before running"""
    if not check_exist(os.path.dirname(os.path.abspath(parser.output)), directory=True):
        sys.exit(('ERROR: Incorrect path specified by --output'))
    elif not check_exist(parser.output, directory=True):
        os.mkdir(parser.output)
    if parser.result:
        # if only evaluate the pipeline output
        if not parser.evaluate:
            sys.exit(('ERROR: The --result argument requires the --evaluate'))
        if not check_exist(parser.result):
            sys.exit(('ERROR: The pipeline result file specified by --result does not exist!'))
        if not check_exist(parser.evaluate):
            sys.exit(('ERROR: The .fam file specified by --evaluate does not exist!'))
    else:
        if not check_exist(parser.client) and not check_exist(parser.client, directory=True):
            sys.exit(('ERROR: The client data specified by --client does not exist!'))
        if (parser.impute or parser.phase):
            if parser.background and (not check_exist(parser.background)):
                # if do inference for phased or imputed data
                sys.exit(('ERROR: The background VCF file specified by --background does not exist!'))
            return
        if not parser.all:
            if not parser.background:
                sys.exit(('ERROR: The --background argument are required without using "-a" option'))
            if not check_exist(parser.background):
                sys.exit(('ERROR: The background VCF file specified by --background does not exist!'))
    return


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
        _, iid, f, m = i.split()[:4]
        iid = iid.split('_')[1]

        iids.add(iid)
        if f != '0':
            f = f.split('_')[1]
            pedigree.add_edge(f, iid)
        if m != '0':
            m = m.split('_')[1]
            pedigree.add_edge(m, iid)
    return iids, pedigree


def read_pipeline_output(fn, only_client=False):
    """Given our pipeline output, return a networkx object"""
    clients = set()
    g = nx.Graph()
    with open(fn) as f:
        next(f)
        for line in f:
            items = line.strip().split(sep="\t")
            clients.add(items[1])
            if only_client:
                clients.add(items[3])
            if items[-1] != 'NA':
                #print('items: ', items)
                g1, g2 = items[1], items[3]
                g.add_edge(g1, g2, ersa=items[-4], king=items[-3])
    return g, clients


def get_kinship(pedigree):
    """Given a networkx object, get all the pairwise relationships"""
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
                    pass  # below will print not only the lowest common ancestry
                    # print("{} and {} has ancestor {} and {}".format(i, j, relatives[i][j]['common_ancestor'], ancestor))
                    # relatives[i][j]['common_ancestor2'] = ancestor
            else:
                relatives.add_edge(i, j, common_ancestor=ancestor, degree=degree)
        # add vertical relationship
        for m in descendants:
            relatives.add_edge(ancestor, m, common_ancestor=ancestor, degree=descendants[m])

    return v_relatives, relatives


def read_across_kin(name):
    relation = nx.Graph()
    with open(name) as f:
        next(f)
        for line in f:
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


def read_ersa(filename):
    relation = nx.Graph()
    with open(filename) as f:
        next(f)
        for line in f:
            items = line.strip().split(sep="\t")
            items = [item.strip() for item in items if item != '']
            if items[2] != 'NA' or items[3] != 'NA':
                relation.add_edge(items[0], items[1], Rel_est1=items[2], Rel_est2=items[3], d_est=items[4])
    return relation


def line_generator(matchfiles):
    for i in matchfiles:
        iterator = open(i)
        for line in iterator:
            yield line
