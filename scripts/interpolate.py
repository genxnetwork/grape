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

    eps = 1e-7
    new_cms = []
    for i in range(0, cms.shape[0]):
        if i == 0:
            new_cms.append(cms[i])
            continue

        if new_cms[-1] >= cms[i]:
            new_cms.append(new_cms[-1] + eps)
        else:
            new_cms.append(cms[i])

    print(f'cms of len {len(cms)} after:')

    print(new_cms[-10:])
    with open(output_path, 'w') as w_map:
        for i, cm in enumerate(new_cms):
            if i > 0 and cm <= new_cms[i-1]:
                raise ValueError(f'cm positions {i-1} and {i} are NOT strictly increasing: {new_cms[:i+1]}')
            w_map.write(f"{i}\t{cm:.8f}\n")


if __name__ == '__main__':
    vcf_path = snakemake.input['vcf']
    cm_path = snakemake.input['cmmap']

    chrom = snakemake.wildcards['chrom']

    output_path = snakemake.output[0]
    print(vcf_path)
    interpolate(vcf_path, cm_path, output_path, chrom)