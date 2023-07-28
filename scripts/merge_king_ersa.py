import pandas
import numpy
import os
import logging
from collections import Counter, namedtuple
from utils.ibd import read_king_segments as rks, interpolate_all, Segment


def is_non_zero_file(fpath):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


def _sort_ids(data):
    unsorted_mask = data.id2 < data.id1
    data.loc[unsorted_mask, 'id1'], data.loc[unsorted_mask, 'id2'] = data.loc[unsorted_mask, 'id2'], data.loc[
        unsorted_mask, 'id1']


def read_bucket_dir(bucket_dir):

    total = None
    for file in os.listdir(bucket_dir):
        if not file.endswith('tsv'):
            continue
        path = os.path.join(bucket_dir, file)
        bucket = read_germline(path)
        if total is None:
            total = bucket
        else:
            total = total.append(bucket)
    return total


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
        data.rename({'PropIBD': 'shared_genome_proportion', 'IBD1Seg': 'king_ibd1_prop', 'IBD2Seg': 'king_ibd2_prop'}, axis='columns', inplace=True)
        indexed = data.loc[:, ['id1', 'id2', 'king_degree', 'king_relation', 'king_ibd1_prop', 'king_ibd2_prop', 'shared_genome_proportion']].\
            set_index(['id1', 'id2'])
        return indexed

    except pandas.errors.EmptyDataError:
        return pandas.DataFrame(data=None, index=None,
                                columns=['id1',
                                         'id2',
                                         'king_degree',
                                         'king_relation',
                                         'king_ibd1_prop',
                                         'king_ibd2_prop',
                                         'shared_genome_proportion']).set_index(['id1', 'id2'])


def read_king_segments_chunked(king_segments_path, map_dir):
    total = None
    try:
        for i, chunk in enumerate(pandas.read_table(
                                                    king_segments_path, 
                                                    compression='gzip', 
                                                    dtype={'FID1': str, 'ID1': str, 'FID2': str, 'ID2': str}, 
                                                    chunksize=1e+6)):

            segments = {}
            for i, row in chunk.iterrows():
                id1 = row['FID1'] + '_' + row['ID1']
                id2 = row['FID2'] + '_' + row['ID2']
                seg = Segment(id1, id2, row['Chr'],
                            bp_start=row['StartMB']*1e+6, bp_end=row['StopMB']*1e+6)
                key = tuple(sorted((seg.id1, seg.id2)))
                if key not in segments:
                    segments[key] = [seg]
                else:
                    segments[key].append(seg)

            segments = interpolate_all(segments, map_dir)
            data = pandas.DataFrame(columns=['id1', 'id2', 'total_seg_len_king', 'seg_count_king'])
            logging.info(f'loaded and interpolated segments from chunk {i} for {len(segments)} pairs')
            rows = []
            for key, segs in segments.items():
                row = {'id1': key[0],
                    'id2': key[1],
                    'total_seg_len_king': sum([s.cm_len for s in segs]),
                    'seg_count_king': len(segs)}
                rows.append(row)
            
            rows_frame = pandas.DataFrame.from_records(rows, columns=['id1', 'id2', 'total_seg_len_king', 'seg_count_king'])
            data = pandas.concat([data, rows_frame], ignore_index=True)

            _sort_ids(data)
            data = data.set_index(['id1', 'id2'])
            total = total.append(data) if total is not None else data

        return total

    except pandas.errors.EmptyDataError:
        return pandas.DataFrame(columns=['id1', 'id2', 'total_seg_len_king', 'seg_count_king']).set_index(['id1', 'id2'])


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


def kinship_to_degree(kinship):
    degrees = []
    # intervals are from KING manual
    # >0.354, [0.177, 0.354], [0.0884, 0.177] and [0.0442, 0.0884]
    for k in kinship:
        if k > 0.354:
            degrees.append(0)
        elif 0.177 < k <= 0.354:
            degrees.append(1)
        elif 0.0884 < k <= 0.177:
            degrees.append(2)
        elif 0.0442 < k <= 0.0884:
            degrees.append(3)
        else:
            degrees.append(None)
    return degrees


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
            data.loc[:, 'kinship_degree'] = kinship_to_degree(data.kinship)
            data = data.loc[:, ['id1', 'id2', 'kinship', 'kinship_degree']].set_index(['id1', 'id2'])
        else:
            data = pandas.DataFrame(columns=['id1', 'id2', 'kinship', 'kinship_degree']).set_index(['id1', 'id2'])
        return data
    except pandas.errors.EmptyDataError:
        return pandas.DataFrame(columns=['id1', 'id2', 'kinship', 'kinship_degree']).set_index(['id1', 'id2'])


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
    # Indv_1     Indv_2      Rel_est1      Rel_est2      d_est     lower_d  upper_d     N_seg     Tot_cM
    data = pandas.read_table(ersa_path, header=0,
                             names=['id1', 'id2', 'rel_est1', 'rel_est2',
                                    'ersa_degree', 'ersa_lower_bound', 'ersa_upper_bound', 'seg_count_ersa', 'total_seg_len_ersa'],
                             dtype={'ersa_degree': str, 'seg_count_ersa': str, 'total_seg_len_ersa': str})

    data = data.loc[(data.rel_est1 != 'NA') | (data.rel_est2 != 'NA'), :]
    data.loc[:, 'id1'] = data.id1.str.strip()
    data.loc[:, 'id2'] = data.id2.str.strip()
    _sort_ids(data)
    data.loc[:, 'ersa_degree'] = pandas.to_numeric(data.ersa_degree.str.strip(), errors='coerce').astype(
        pandas.Int32Dtype())
    data.loc[:, 'ersa_lower_bound'] = pandas.to_numeric(data.ersa_lower_bound.str.strip(), errors='coerce').astype(
        pandas.Int32Dtype())
    data.loc[:, 'ersa_upper_bound'] = pandas.to_numeric(data.ersa_upper_bound.str.strip(), errors='coerce').astype(
        pandas.Int32Dtype())
    print('dtypes are: ', data.dtypes)
    data.loc[:, 'seg_count_ersa'] = pandas.to_numeric(data.seg_count_ersa.str.strip(), errors='coerce').astype(
        pandas.Int32Dtype())
    data.loc[:, 'total_seg_len_ersa'] = pandas.to_numeric(data.total_seg_len_ersa.str.strip(), errors='coerce')

    data.loc[data.ersa_lower_bound != 1, 'ersa_lower_bound'] = data.ersa_lower_bound - 1
    data.loc[data.ersa_upper_bound != 1, 'ersa_upper_bound'] = data.ersa_upper_bound - 1

    logging.info(f'read {data.shape[0]} pairs from ersa output')
    logging.info(f'ersa after reading has {pandas.notna(data.ersa_degree).sum()}')
    data.loc[:, 'is_niece_aunt'] = [True if 'Niece' in d else False for d in data.rel_est1]
    logging.info(f'we have total of {data.is_niece_aunt.sum()} possible Niece/Aunt relationships')

    return data.loc[data.id1 != data.id2, ['id1', 'id2', 'ersa_degree', 'ersa_lower_bound', 'ersa_upper_bound', 'is_niece_aunt', 'seg_count_ersa', 'total_seg_len_ersa']].\
        set_index(['id1', 'id2'])


def read_ersa2(ersa_path):
    # individual_1    individual_2    est_number_of_shared_ancestors  est_degree_of_relatedness       0.95 CI_2p_lower:4
    # 2p_upper:5     1p_lower:6 1p_upper:7        0p_lower:8        0p_upper:9
    # maxlnl_relatedness      maxlnl_unrelatednes
    data = pandas.read_table(ersa_path, comment='#')

    data = data.loc[data.est_degree_of_relatedness != 'no_sig_rel', :]
    data.rename({'individual_1': 'id1', 'individual_2': 'id2', 'est_number_of_shared_ancestors': 'shared_ancestors'}, axis=1, inplace=True)
    data.loc[:, 'ersa_degree'] = pandas.to_numeric(data['est_degree_of_relatedness'], errors='coerce').\
        astype(pandas.Int32Dtype())

    lower, upper = [], []
    for i, ancestors in enumerate(data.shared_ancestors):
        if ancestors == 0:
            lower.append(data.iloc[i, 8])
            upper.append(data.iloc[i, 9])
        if ancestors == 1:
            lower.append(data.iloc[i, 6])
            upper.append(data.iloc[i, 7])
        if ancestors == 2:
            lower.append(data.iloc[i, 4])
            upper.append(data.iloc[i, 5])

    data.loc[:, 'ersa_lower_bound'] = lower
    data.loc[:, 'ersa_upper_bound'] = upper

    cols = ['id1', 'id2', 'ersa_degree', 'ersa_lower_bound', 'ersa_upper_bound', 'shared_ancestors']
    return data.loc[data.id1 != data.id2, cols].set_index(['id1', 'id2'])


if __name__ == '__main__':

    try:
        snakemake
    except NameError:
        test_dir = '/media_ssd/pipeline_data/TF-CEU-TRIBES-ibis-king-2'
        Snakemake = namedtuple('Snakemake', ['input', 'output', 'params', 'log'])
        snakemake = Snakemake(
            input={'king': f'{test_dir}/king/data.seg',
                   'king_segments': f'{test_dir}/king/data.segments.gz',
                   'kinship': f'{test_dir}/king/data.kin',
                   'kinship0': f'{test_dir}/king/data.kin0',
                   'ersa': f'{test_dir}/ersa/relatives.tsv'},

            output=[f'test_data/merge_king_ersa.tsv'],
            params={'cm_dir': f'{test_dir}/cm'},
            log=['test_data/merge_king_ersa.log']
        )

    logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG, format='%(levelname)s:%(asctime)s %(message)s')

    king_path = snakemake.input['king']
    king_segments_path = snakemake.input['king_segments']
    # within families
    kinship_path = snakemake.input['kinship']
    # across families
    kinship0_path = snakemake.input['kinship0']
    ersa_path = snakemake.input['ersa']
    map_dir = snakemake.params['cm_dir']
    output_path = snakemake.output[0]

    king = read_king(king_path)
    kinship = read_kinship(kinship_path, kinship0_path)
    king_segments = read_king_segments_chunked(king_segments_path, map_dir)
    ersa = read_ersa(ersa_path)

    relatives = king.merge(kinship, how='outer', left_index=True, right_index=True).\
                     merge(ersa, how='outer', left_index=True, right_index=True).\
                     merge(king_segments, how='outer', left_index=True, right_index=True)

    prefer_ersa_mask = pandas.isnull(relatives.king_degree) | (relatives.king_degree > 3)
    relatives.loc[:, 'final_degree'] = relatives.king_degree
    # if king is unsure or king degree > 3 then we use ersa distant relatives estimation
    relatives.loc[prefer_ersa_mask, 'final_degree'] = relatives.ersa_degree
    logging.info(f'king is null or more than 3: {prefer_ersa_mask.sum()}')
    logging.info(f'ersa is not null: {pandas.notna(relatives.ersa_degree).sum()}')
    print(f'relatives columns are {relatives.columns}')
    if 'total_seg_len_king' in relatives.columns:
        relatives.loc[:, 'total_seg_len'] = relatives.total_seg_len_king*relatives.king_ibd1_prop
        relatives.loc[:, 'total_seg_len_ibd2'] = relatives.total_seg_len_king*relatives.king_ibd2_prop
        relatives.loc[:, 'seg_count'] = relatives.seg_count_king

    relatives.loc[prefer_ersa_mask, 'total_seg_len'] = relatives.total_seg_len_ersa
    relatives.loc[prefer_ersa_mask, 'seg_count'] = relatives.seg_count_ersa
    print('is na: ', pandas.isna(relatives.loc[prefer_ersa_mask, 'total_seg_len_ersa']).sum())
    # approximate calculations, IBD share is really small in this case
    relatives.loc[prefer_ersa_mask, 'shared_genome_proportion'] = 0.5*relatives.loc[prefer_ersa_mask, 'total_seg_len'].values / 3440
    relatives.drop(['total_seg_len_king', 'seg_count_king', 'total_seg_len_ersa', 'seg_count_ersa', 'king_ibd1_prop', 'king_ibd2_prop'],
                   axis='columns', inplace=True)

    logging.info(f'final degree not null: {pandas.notna(relatives.final_degree).sum()}')

    relatives.loc[pandas.notna(relatives.final_degree), :].to_csv(output_path, sep='\t')
