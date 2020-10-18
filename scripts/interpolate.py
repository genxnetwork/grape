import pandas
import numpy
from utils.bcftools import bcftools_query


def interpolate(vcf_path, cm_path, output_path, chrom):

    genetic_map = pandas.read_csv(cm_path, sep=" ")

    pos = list(map(int, bcftools_query(vcf_path, "%POS\n")))
    # interpolate CM from BP
    #chro = genetic_map[genetic_map.chr == int(chrom)]
    # because for chr19 'Genetic_Map(cM)' is 'Genetic_Map.cM.'
    cms = numpy.interp(pos, genetic_map.position.values, genetic_map.iloc[:, -1].values)

    with open(output_path, 'w') as w_map:
        for i, cm in enumerate (cms):
            w_map.write("{}\t{}\n".format(i, cm))


if __name__ == '__main__':
    vcf_path = snakemake.input['vcf']
    cm_path = snakemake.input['cmmap']

    chrom = snakemake.wildcards['chrom']

    output_path = snakemake.output[0]
    print(vcf_path)
    interpolate(vcf_path, cm_path, output_path, chrom)