import pandas
import numpy
import os
import sys
from utils.ibd import read_king_segments as rks, interpolate_all


def is_non_zero_file(fpath):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


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
    data.loc[:, 'id1'] = data.fid_iid1.str.split().str[1]
    data.loc[:, 'id2'] = data.fid_iid2.str.split().str[1]

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
    data = pandas.read_table(king_path)
    data.loc[:, 'id1'] = data.FID1.astype(str) + '_' + data.ID1.astype(str)
    data.loc[:, 'id2'] = data.FID2.astype(str) + '_' + data.ID2.astype(str)

    data.loc[:, 'king_degree'] = map_king_degree(data.InfType)
    data.loc[:, 'king_relation'] = map_king_relation(data.InfType)
    data.loc[:, 'king_degree'] = data.king_degree.astype(float).astype(pandas.Int32Dtype())
    data.rename({'PropIBD': 'shared_genome_proportion'}, axis='columns', inplace=True)
    indexed = data.loc[:, ['id1', 'id2', 'king_degree', 'king_relation', 'shared_genome_proportion']].\
        set_index(['id1', 'id2'])

    return indexed


def read_king_segments(king_segments_path, map_dir):
    segments = rks(king_segments_path)
    segments = interpolate_all(segments, map_dir)
    data = pandas.DataFrame(columns=['id1', 'id2', 'total_seg_len_king', 'seg_count_king'])
    print(f'loaded and interpolated segments for {len(segments)} pairs')
    for key, segs in segments.items():
        row = {'id1': key[0],
               'id2': key[1],
               'total_seg_len_king': sum([s.cm_len for s in segs]),
               'seg_count_king': len(segs)}
        data = data.append(row, ignore_index=True)

    return data.set_index(['id1', 'id2'])


def read_kinship(kinship_path, kinship0_path):
    # parse within families
    # FID     ID1     ID2     N_SNP   Z0      Phi     HetHet  IBS0    Kinship Error
    within, across = None, None
    if is_non_zero_file(kinship_path):
        within = pandas.read_table(kinship_path)
        within.loc[:, 'id1'] = within.FID.astype(str) + '_' + within.ID1.astype(str)
        within.loc[:, 'id2'] = within.FID.astype(str) + '_' + within.ID2.astype(str)
        within.rename({'Kinship': 'kinship'}, axis=1, inplace=True)
        within = within.loc[:, ['id1', 'id2', 'kinship']].set_index(['id1', 'id2'])
        print(f'loaded {within.shape[0]} pairs from within-families kinship estimation results')

    # FID1    ID1     FID2    ID2     N_SNP   HetHet  IBS0    Kinship
    if is_non_zero_file(kinship0_path):
        across = pandas.read_table(kinship0_path)
        across.loc[:, 'id1'] = across.FID1.astype(str) + '_' + across.ID1.astype(str)
        across.loc[:, 'id2'] = across.FID2.astype(str) + '_' + across.ID2.astype(str)
        across.rename({'Kinship': 'kinship'}, axis=1, inplace=True)
        across = across.loc[:, ['id1', 'id2', 'kinship']].set_index(['id1', 'id2'])
        print(f'loaded {across.shape[0]} pairs from across-families kinship estimation results')

    if within is None and across is None:
        return None
    elif within is None and across is not None:
        return across
    elif within is not None and across is None:
        return within
    else:
        return pandas.concat([within, across], axis=0)


def read_ersa(ersa_path):
    # Indv_1     Indv_2      Rel_est1      Rel_est2      d_est     N_seg     Tot_cM
    data = pandas.read_table(ersa_path, header=0,
                             names=['id1', 'id2', 'rel_est1', 'rel_est2', 'ersa_degree', 'seg_count', 'total_seg_len'])

    data = data.loc[(data.rel_est1 != 'NA') | (data.rel_est2 != 'NA'), :]
    data.loc[:, 'id1'] = data.id1.str.strip()
    data.loc[:, 'id2'] = data.id2.str.strip()
    data.loc[:, 'ersa_degree'] = pandas.to_numeric(data.ersa_degree.str.strip(), errors='coerce').astype(pandas.Int32Dtype())

    print(f'read {data.shape[0]} pairs from ersa output')

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

    if kinship is not None:
        print(f'kinship is not none')
        print(kinship.columns)
        relatives = ibd.merge(king, how='outer', left_index=True, right_index=True).\
            merge(kinship, how='outer', left_index=True, right_index=True).\
            merge(ersa, how='outer', left_index=True, right_index=True).\
            merge(king_segments, how='outer', left_index=True, right_index=True)
    else:
        print('kinship is none')
        relatives = ibd.merge(king, how='outer', left_index=True, right_index=True).\
            merge(ersa, how='outer', left_index=True, right_index=True). \
            merge(king_segments, how='outer', left_index=True, right_index=True)

    prefer_ersa_mask = pandas.isnull(relatives.king_degree) | (relatives.king_degree > 3)
    relatives.loc[:, 'final_degree'] = relatives.king_degree
    # if king is unsure or king degree > 3 then we use ersa distant relatives estimation
    relatives.loc[prefer_ersa_mask, 'final_degree'] = relatives.ersa_degree

    if 'total_seg_len_king' in relatives.columns:
        relatives.loc[:, 'total_seg_len'] = relatives.total_seg_len_king
        relatives.loc[:, 'seg_count'] = relatives.seg_count_king

    relatives.loc[prefer_ersa_mask, 'total_seg_len'] = relatives.total_seg_len_germline
    relatives.loc[prefer_ersa_mask, 'seg_count'] = relatives.seg_count_germline
    relatives.drop(['total_seg_len_king', 'seg_count_king', 'total_seg_len_germline', 'seg_count_germline'],
                   axis='columns', inplace=True)

    relatives.loc[pandas.notna(relatives.final_degree), :].to_csv(output_path, sep='\t')
