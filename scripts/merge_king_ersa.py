import pandas
import numpy
import os
import logging
from utils.ibd import read_king_segments as rks, interpolate_all


def is_non_zero_file(fpath):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


def _sort_ids(data):
    unsorted_mask = data.id2 < data.id1
    data.loc[unsorted_mask, 'id1'], data.loc[unsorted_mask, 'id2'] = data.loc[unsorted_mask, 'id2'], data.loc[
        unsorted_mask, 'id1']

def read_germline(ibd_path):
    germline_names = [
        'fid_iid1',
        'fid_iid2',
        'chrom',
        'start_end_bp',
        'start_end_snp',
        'snp_len',
        'genetic_len',
        'len_units',
        'mismatches',
        'is_homozygous1',
        'is_homozygous2'
    ]

    data = pandas.read_table(ibd_path, header=None, names=germline_names)
    if data.shape[0] == 0:
        return pandas.DataFrame(columns=['id1', 'id2', 'seg_count_germline', 'total_seg_len_germline']).set_index(['id1', 'id2'])

    data.loc[:, 'id1'] = data.fid_iid1.str.split().str[1]
    data.loc[:, 'id2'] = data.fid_iid2.str.split().str[1]
    _sort_ids(data)
    segments_info = data.loc[:, ['id1', 'id2', 'chrom', 'genetic_len']].groupby(by=['id1', 'id2']).agg({'genetic_len': 'sum', 'chrom': 'count'})
    segments_info.rename({'genetic_len': 'total_seg_len_germline', 'chrom': 'seg_count_germline'}, axis=1, inplace=True)
    return segments_info


def map_king_degree(king_degree):
    degree_map = {
        'Dup/MZ': 0,
        'PO': 1,
        'FS': 2,
        '2nd': 2,
        '3rd': 3,
        '4th': 4,
        'NA': numpy.nan
    }
    return [degree_map[kd] for kd in king_degree]


def map_king_relation(king_degree):
    degree_map = {
        'Dup/MZ': 0,
        'PO': 'PO',
        'FS': 'FS',
        '2nd': 2,
        '3rd': 3,
        '4th': 4,
        'NA': numpy.nan
    }
    return [degree_map[kd] for kd in king_degree]


def read_king(king_path):
    # FID1    ID1     FID2    ID2     MaxIBD1 MaxIBD2 IBD1Seg IBD2Seg PropIBD InfType
    try:
        data = pandas.read_table(king_path)
        data.loc[:, 'id1'] = data.FID1.astype(str) + '_' + data.ID1.astype(str)
        data.loc[:, 'id2'] = data.FID2.astype(str) + '_' + data.ID2.astype(str)
        _sort_ids(data)

        data.loc[:, 'king_degree'] = map_king_degree(data.InfType)
        data.loc[:, 'king_relation'] = map_king_relation(data.InfType)
        data.loc[:, 'king_degree'] = data.king_degree.astype(float).astype(pandas.Int32Dtype())
        data.rename({'PropIBD': 'shared_genome_proportion'}, axis='columns', inplace=True)
        indexed = data.loc[:, ['id1', 'id2', 'king_degree', 'king_relation', 'shared_genome_proportion']].\
            set_index(['id1', 'id2'])
        return indexed

    except pandas.errors.EmptyDataError:
        return pandas.DataFrame(data=None, index=None,
                                columns=['id1',
                                         'id2',
                                         'king_degree',
                                         'king_relation',
                                         'shared_genome_proportion']).set_index(['id1', 'id2'])


def read_king_segments(king_segments_path, map_dir):
    try:
        segments = rks(king_segments_path)
        segments = interpolate_all(segments, map_dir)
        data = pandas.DataFrame(columns=['id1', 'id2', 'total_seg_len_king', 'seg_count_king'])
        logging.info(f'loaded and interpolated segments for {len(segments)} pairs')
        for key, segs in segments.items():
            row = {'id1': key[0],
                   'id2': key[1],
                   'total_seg_len_king': sum([s.cm_len for s in segs]),
                   'seg_count_king': len(segs)}
            data = data.append(row, ignore_index=True)

        _sort_ids(data)
        return data.set_index(['id1', 'id2'])

    except pandas.errors.EmptyDataError:
        return pandas.DataFrame(columns=['id1', 'id2', 'total_seg_len_king', 'seg_count_king']).set_index(['id1', 'id2'])


def _read_kinship_data(kinship_path):
    try:
        data = pandas.read_table(kinship_path)
        print(kinship_path, data.shape)
        # If no relations were found, king creates file with only header
        if data.shape[0] != 0:
            if 'FID' in data.columns:
                data.loc[:, 'id1'] = data.FID.astype(str) + '_' + data.ID1.astype(str)
                data.loc[:, 'id2'] = data.FID.astype(str) + '_' + data.ID2.astype(str)
            else:
                data.loc[:, 'id1'] = data.FID1.astype(str) + '_' + data.ID1.astype(str)
                data.loc[:, 'id2'] = data.FID2.astype(str) + '_' + data.ID2.astype(str)

            _sort_ids(data)
            data.rename({'Kinship': 'kinship'}, axis=1, inplace=True)
            data = data.loc[:, ['id1', 'id2', 'kinship']].set_index(['id1', 'id2'])
        else:
            data = pandas.DataFrame(columns=['id1', 'id2', 'kinship']).set_index(['id1', 'id2'])
        return data
    except pandas.errors.EmptyDataError:
        return pandas.DataFrame(columns=['id1', 'id2', 'kinship']).set_index(['id1', 'id2'])


def read_kinship(kinship_path, kinship0_path):
    # parse within families
    # FID     ID1     ID2     N_SNP   Z0      Phi     HetHet  IBS0    Kinship Error
    within = _read_kinship_data(kinship_path)
    logging.info(f'loaded {within.shape[0]} pairs from within-families kinship estimation results')

    # FID1    ID1     FID2    ID2     N_SNP   HetHet  IBS0    Kinship
    across = _read_kinship_data(kinship0_path)
    logging.info(f'loaded {across.shape[0]} pairs from across-families kinship estimation results')
    return pandas.concat([within, across], axis=0)


def read_ersa(ersa_path):
    # Indv_1     Indv_2      Rel_est1      Rel_est2      d_est     N_seg     Tot_cM
    data = pandas.read_table(ersa_path, header=0,
                             names=['id1', 'id2', 'rel_est1', 'rel_est2', 'ersa_degree', 'seg_count', 'total_seg_len'],
                             dtype={'ersa_degree': str})

    data = data.loc[(data.rel_est1 != 'NA') | (data.rel_est2 != 'NA'), :]
    data.loc[:, 'id1'] = data.id1.str.strip()
    data.loc[:, 'id2'] = data.id2.str.strip()
    _sort_ids(data)
    data.loc[:, 'ersa_degree'] = pandas.to_numeric(data.ersa_degree.str.strip(), errors='coerce').astype(pandas.Int32Dtype())

    logging.info(f'read {data.shape[0]} pairs from ersa output')
    logging.info(f'ersa after reading has {pandas.notna(data.ersa_degree).sum()}')
    return data.loc[data.id1 != data.id2, ['id1', 'id2', 'ersa_degree']].set_index(['id1', 'id2'])


if __name__ == '__main__':

    '''
    ibd_path = 'test_data/merge_king_ersa/merged_ibd.tsv'
    king_path = 'test_data/merge_king_ersa/merged_imputed_king.seg'
    king_segments_path = 'test_data/merge_king_ersa/merged_imputed_king.segments.gz'
    # within families
    kinship_path = 'test_data/merge_king_ersa/merged_imputed_kinship.kin'
    # across families
    kinship0_path = 'test_data/merge_king_ersa/merged_imputed_kinship.kin0'
    ersa_path = 'test_data/merge_king_ersa/relatives.tsv'
    map_dir = '/media/pipeline_data/sim-vcf-to-ped/cm'
    output_path = 'test_data/relatives_merge_king_ersa.tsv'
    '''
    logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG, format='%(levelname)s:%(asctime)s %(message)s')

    ibd_path = snakemake.input['ibd']
    king_path = snakemake.input['king']
    king_segments_path = snakemake.input['king_segments']
    # within families
    kinship_path = snakemake.input['kinship']
    # across families
    kinship0_path = snakemake.input['kinship0']
    ersa_path = snakemake.input['ersa']
    map_dir = snakemake.params['cm_dir']
    output_path = snakemake.output[0]

    ibd = read_germline(ibd_path)
    king = read_king(king_path)
    kinship = read_kinship(kinship_path, kinship0_path)
    king_segments = read_king_segments(king_segments_path, map_dir)
    ersa = read_ersa(ersa_path)

    logging.info(f'ibd shape: {ibd.shape[0]}, ersa shape: {ersa.shape[0]}')
    # print('ibd test:',  ibd[('GRC12118091', 'GRC12118096')])
    # print('ibd test2:', ibd[('GRC12118096', 'GRC12118091')])

    relatives = ibd.merge(king, how='outer', left_index=True, right_index=True).\
                    merge(kinship, how='outer', left_index=True, right_index=True).\
                    merge(ersa, how='outer', left_index=True, right_index=True).\
                    merge(king_segments, how='outer', left_index=True, right_index=True)

    prefer_ersa_mask = pandas.isnull(relatives.king_degree) | (relatives.king_degree > 3)
    relatives.loc[:, 'final_degree'] = relatives.king_degree
    # if king is unsure or king degree > 3 then we use ersa distant relatives estimation
    relatives.loc[prefer_ersa_mask, 'final_degree'] = relatives.ersa_degree
    logging.info(f'king is null or more than 3: {prefer_ersa_mask.sum()}')
    logging.info(f'ersa is not null: {pandas.notna(relatives.ersa_degree).sum()}')

    if 'total_seg_len_king' in relatives.columns:
        relatives.loc[:, 'total_seg_len'] = relatives.total_seg_len_king
        relatives.loc[:, 'seg_count'] = relatives.seg_count_king

    relatives.loc[prefer_ersa_mask, 'total_seg_len'] = relatives.total_seg_len_germline
    relatives.loc[prefer_ersa_mask, 'seg_count'] = relatives.seg_count_germline

    # approximate calculations, IBD share is really small in this case
    relatives.loc[prefer_ersa_mask, 'shared_genome_proportion'] = 0.5*relatives.loc[prefer_ersa_mask, 'total_seg_len'] / 3580
    relatives.drop(['total_seg_len_king', 'seg_count_king', 'total_seg_len_germline', 'seg_count_germline'],
                   axis='columns', inplace=True)

    logging.info(f'final degree not null: {pandas.notna(relatives.final_degree).sum()}')

    relatives.loc[pandas.notna(relatives.final_degree), :].to_csv(output_path, sep='\t')
