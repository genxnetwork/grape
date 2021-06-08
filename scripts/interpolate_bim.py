import pandas
import numpy
from collections import namedtuple


def interpolate(bim_path, cm_path, output_path):

    genetic_map = pandas.read_csv(cm_path, sep=" ", compression='gzip')
    bim = pandas.read_table(bim_path, header=None, names=['chrom', '_id', 'cm', 'pos', 'ref', 'alt'])

    for chrom, group in genetic_map.groupby(by='chr'):
        cms = numpy.interp(bim.pos[bim.chrom == chrom].values, group.position.values, group.iloc[:, -1].values)
        bim.loc[bim.chrom == chrom, 'cm'] = cms

    bim.to_csv(output_path, sep='\t', header=None, index=False)


if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        test_dir = '/media_ssd/pipeline_data/TF-CEU-TRIBES-ibis-king-2'
        Snakemake = namedtuple('Snakemake', ['input', 'output'])
        snakemake = Snakemake(
            input={'bim': f'{test_dir}/preprocessed/data_unmapped.bim',
                   'cmmap': f'/media_ssd/new_ref/tables/genetic_map_hg19_withX.txt.gz'},

            output=[f'test_data/interpolated.bim']
        )

    bim_path = snakemake.input['bim']
    cm_path = snakemake.input['cmmap']

    output_path = snakemake.output[0]
    interpolate(bim_path, cm_path, output_path)