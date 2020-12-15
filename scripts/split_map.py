import pandas
import os


if __name__ == '__main__':
    path = snakemake.input['bim']
    map_dir = snakemake.params['cm_dir']
    bim = pandas.read_table(path, header=None, names=['chrom', '_id', 'cm_pos', 'bp_pos', 'ref', 'alt'])

    maps = bim.loc[:, ['chrom', '_id', 'cm_pos', 'bp_pos']].groupby(by='chrom')
    for chrom, cmmap in maps:
        map_path = os.path.join(map_dir, f'chr{chrom}.cm.map')

        frame = cmmap.reset_index()
        frame.to_csv(map_path, index=False, sep='\t', header=None)
        print(f'written map for chrom {chrom} to {map_path}')

    print(f'written {len(maps)} maps to {map_dir}')
