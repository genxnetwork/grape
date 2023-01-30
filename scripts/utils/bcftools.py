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


def bcftools_query(infile, arg, save_output=False):
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

    if save_output:
        with open(save_output, 'w') as output:
            output.write(std_out)

    return std_out.decode().split()


def bcftools_view(infile, outfile, samples):
    """Return one column of snp information of a given vcf file
    """
    print(infile)
    if not os.path.exists(infile + '.csi'):
        bcftools_index(infile)
    args = ["bcftools", "view", "-S", samples, infile, "-O", "z", "-o", outfile, "--force-samples"]
    pipes = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std_out, std_err = pipes.communicate()
    if pipes.returncode != 0:
        print(std_out.decode())
        raise Exception(std_err.decode())


def bcftools_stats(infile, samples, save_output=False, save_stats=''):
    """Returns vcf file stats as raw string table
    """
    print(infile)
    if not os.path.exists(infile + '.csi'):
        bcftools_index(infile)
    args = ["bcftools", "stats", "-S", samples, infile,
            "|", "tee", save_stats, "|",
            "|", "grep", "'^PSC'"]
    pipes = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std_out, std_err = pipes.communicate()
    if pipes.returncode != 0:
        print(std_out.decode())
        raise Exception(std_err.decode())

    if save_output:
        with open(save_output, 'w') as output:
            output.write(std_out)

    return std_out.decode().split()
