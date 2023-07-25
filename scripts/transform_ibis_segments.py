import pandas
import logging
import mmh3
from typing import List
import os
from collections import namedtuple
import subprocess


def process_chunk_with_hash(data: pandas.DataFrame, denominator: int, dest_dir: str):
    data.loc[:, 'id1'] = data.sample1.str.replace(':', '_') + ' ' + data.sample1.str.replace(':', '_')
    data.loc[:, 'id2'] = data.sample2.str.replace(':', '_') + ' ' + data.sample2.str.replace(':', '_')
    # awk '{{sub(":", "_", $1); sub(":", "_", $2); print $1, $1 "\t" $2, $2 "\t" $3 "\t" $4, $5 "\t" 0, 0 "\t" $10 "\t" $9 "\t" "cM" "\t" 0 "\t" 0 "\t" 0}};
    data.loc[:, 'phys_pos'] = [f'{start} {end}' for start, end in zip(data.phys_start_pos, data.phys_end_pos)]
    data.loc[:, 'zero1'] = '0 0'
    data.loc[:, 'zero2'] = '0'
    data.loc[:, 'cM'] = 'cM'
    data.loc[:, 'ibd21'] = [2 if ibd_type == 'IBD2' else 0 for ibd_type in data.IBD_type]
    data.loc[:, 'ibd22'] = data.ibd21
    data.loc[:, 'bucket_id'] = [mmh3.hash(id1) % denominator for id1 in data.id1]

    to_write = data.loc[:,
               ['bucket_id', 'id1', 'id2', 'chrom', 'phys_pos', 'zero1', 'marker_count', 'genetic_seg_length', 'cM', 'zero2',
                'ibd21', 'ibd22']]

    for bucket_id, group in to_write.groupby('bucket_id'):
        dest_file = os.path.join(dest_dir, f'{bucket_id}.tsv')
        group.drop('bucket_id', inplace=True, axis='columns')
        group.to_csv(dest_file, index=False, header=None, sep='\t', mode='a')


def split_by_id(input_ibd: str, samples_count: int, dest_dir: str):
    names = [
        'sample1',
        'sample2',
        'chrom',
        'phys_start_pos',
        'phys_end_pos',
        'IBD_type',
        'genetic_start_pos',
        'genetic_end_pos',
        'genetic_seg_length',
        'marker_count',
        'error_count',
        'error_density'
    ]
    read_chunksize = int(1e+6)
    samples_chunksize = 2000
    denominator = samples_count // samples_chunksize + 1

    for i, chunk in enumerate(pandas.read_csv(input_ibd, header=None, names=names, sep='\t', chunksize=read_chunksize)):
        if not chunk.empty:
            process_chunk_with_hash(chunk, denominator, dest_dir)
            logging.info(f'Chunk {i} of size {read_chunksize} was written to {dest_dir} and split into {denominator} buckets')
        else:
            logging.info(f'No IBD segments were found in {input_ibd}. Nothing was written to {dest_dir}')


if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        Snakemake = namedtuple('Snakemake', ['input', 'output', 'log'])
        snakemake = Snakemake(
            input={
                'ibd': '/media/data1/relatives/runs/run100k/ibis/merged_ibis.seg',
                'fam': '/media/data1/relatives/runs/run100k/preprocessed/data.fam'},
            output={'bucket_dir': '/media/data1/relatives/test100koutput.txt'},
            log=['/media/data1/relatives/test100koutput.log']
        )

    print(subprocess.check_output('which python', shell=True).decode('utf-8'))
    print(subprocess.check_output('conda info', shell=True).decode('utf-8'))
    print(subprocess.check_output('conda list', shell=True).decode('utf-8'))

    ibd = snakemake.input['ibd']
    bucket_dir = snakemake.output['bucket_dir']
    if not os.path.exists(bucket_dir):
        os.makedirs(bucket_dir)
    for file in os.listdir(bucket_dir):
        os.remove(os.path.join(bucket_dir, file))

    samples = pandas.read_table(snakemake.input['fam'])
    samples_count = samples.shape[0] + 1 # there is no header in the file
    print(snakemake.log)
    logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG, format='%(levelname)s:%(asctime)s %(message)s')

    split_by_id(ibd, samples_count, bucket_dir)
