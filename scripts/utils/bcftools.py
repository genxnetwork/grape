import subprocess
import os


def bcftools_index(infile):
    args = ["bcftools", 'index', '-f', infile]
    pipes = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std_out, std_err = pipes.communicate()
    if pipes.returncode != 0:
        print(std_out.decode())
        raise Exception(std_err.decode())
    print("Indexing {}...".format(infile))
    return


def bcftools_query(infile, arg):
    """Return one column of snp information of a given vcf file
    arg = "%ID\n" for SNP ids
    arg = "%POS\n" for SNP position
    """
    print(infile)
    if not os.path.exists(infile + '.csi'):
        bcftools_index(infile)
    args = ["bcftools", "query", "-f", arg, infile]
    pipes = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std_out, std_err = pipes.communicate()
    if pipes.returncode != 0:
        print(std_out.decode())
        raise Exception(std_err.decode())
    return std_out.decode().split()
