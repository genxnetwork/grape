import os
import subprocess
import networkx as nx

from itertools import combinations


BCFTOOLS = "bcftools"


def bcftools_index(infile):
    args = [BCFTOOLS, 'index', '--force', '-f', infile]
    pipes = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std_out, std_err = pipes.communicate()
    if pipes.returncode != 0:
        print(std_out.decode())
        raise Exception(std_err.decode())
    return


def bcftools_samples(infile):
    """Return individual ids of given vcf file"""
    #if not os.path.exists(infile + '.csi'):
    args = [BCFTOOLS, 'query', '-l', infile]
    pipes = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std_out, std_err = pipes.communicate()
    if pipes.returncode != 0:
        print(std_out.decode())
        raise Exception(std_err.decode())
    return std_out.decode().split()


def bcftools_reheader(infile, sample, outname=None):
    """Rename vcf samples given a sample file"""
    if not outname:
        outname = infile

    # workaround of https: // github.com / samtools / bcftools / issues / 1288
    vcf_uncompressed = infile[:-3]  # remove .gz
    view_args = [BCFTOOLS, 'view', infile, '-O', 'v', '-o', vcf_uncompressed]
    pipes = subprocess.Popen(view_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std_out, std_err = pipes.communicate()
    if pipes.returncode != 0:
        print(std_out.decode())
        raise Exception(std_err.decode())

    vcf_reheaded = vcf_uncompressed[:-4] + '.reheaded.vcf'
    args = [BCFTOOLS, 'reheader', vcf_uncompressed, '-s', sample, '-o', vcf_reheaded]
    pipes = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std_out, std_err = pipes.communicate()
    if pipes.returncode != 0:
        print(std_out.decode())
        raise Exception(std_err.decode())

    args = [BCFTOOLS, 'view', vcf_reheaded, '-Oz', '-o', outname]
    pipes = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std_out, std_err = pipes.communicate()
    if pipes.returncode != 0:
        print(std_out.decode())
        raise Exception(std_err.decode())

    os.remove(vcf_reheaded)
    os.remove(vcf_uncompressed)
    return outname


def read_pedigree(fn = './hapmap_sim.fam'):
    pedigree = nx.DiGraph()
    for i in open(fn):
        _, iid ,f ,m = i.split()[:4]
        if f != '0':
            pedigree.add_edge(f, iid)
        if m != '0':
            pedigree.add_edge(m, iid)
    return pedigree


def get_kinship(pedigree):
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
                    #print("{} and {} has ancestor {} and {}".format(i, j, relatives[i][j]['common_ancestor'], ancestor))
                    #relatives[i][j]['common_ancestor2'] = ancestor
            else:
                relatives.add_edge(i, j, common_ancestor=ancestor, degree=degree)
        # add vertical relationship
        for m in descendants:
            relatives.add_edge(ancestor, m, common_ancestor=ancestor, degree=descendants[m])

    return v_relatives, relatives


if __name__ == "__main__":
    print('input:')
    print(snakemake.input)
    print(snakemake.input.keys())
    fam_name = snakemake.input['fam']
    fam_output = snakemake.output['fam']
    with open(fam_name, 'r') as old_file, open(fam_output, 'w') as new_file:
        old_ids = old_file.readlines()
        for i, old_id in enumerate(old_ids):
            # make names legal for input
            #new_file.write(old_id.replace("-", "").replace("_", ""))
            new_file.write(old_id)

    print()

    _, kinship = get_kinship(read_pedigree(fam_name))
    kin_name = snakemake.output['kin']
    nx.write_edgelist(kinship, kin_name, data=True)

    vcf_name = snakemake.input['vcf']
    vcf_output = snakemake.output['vcf']
    sample_ids = bcftools_samples(vcf_name)
    with open(vcf_output + '.samples', 'w') as samples_file:
        for i, _id in enumerate(sample_ids):
            i1, i2 = _id.split(sep='_')
            samples_file.write("{}_{}\n".format(i1, i2))

    filename = bcftools_reheader(vcf_name, vcf_output + '.samples', vcf_output)
    bcftools_index(vcf_output)
