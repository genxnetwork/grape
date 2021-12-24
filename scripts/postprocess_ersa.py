import pandas
import numpy
import os
import logging
from collections import Counter
from utils.ibd import read_king_segments as rks, interpolate_all


def is_non_zero_file(fpath):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


def _sort_ids(data):
    unsorted_mask = data.id2 < data.id1
    data.loc[unsorted_mask, 'id1'], data.loc[unsorted_mask, 'id2'] = data.loc[unsorted_mask, 'id2'], data.loc[
        unsorted_mask, 'id1']


def read_ibis(ibd_path):
    # sample1 sample2 chrom phys_start_pos phys_end_pos IBD_type
    # genetic_start_pos genetic_end_pos genetic_seg_length marker_count error_count error_density
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
    data = pandas.read_table(ibd_path, header=None, names=names)
    data.loc[:, 'id1'] = data.sample1.str.replace(':', '_')
    data.loc[:, 'id2'] = data.sample2.str.replace(':', '_')
    _sort_ids(data)
    ibd2_info = data.loc[data.IBD_type == 'IBD2', ['id1', 'id2', 'chrom', 'genetic_seg_length']].groupby(by=['id1', 'id2']).agg(
        {'genetic_seg_length': 'sum', 'chrom': 'count'})
    ibd2_info.rename({'genetic_seg_length': 'total_seg_len_ibd2', 'chrom': 'seg_count_ibd2'}, axis=1, inplace=True)
    print(ibd2_info)
    return ibd2_info


def read_ersa(ersa_path):
    # Indv_1     Indv_2      Rel_est1      Rel_est2      d_est     lower_d  upper_d     N_seg     Tot_cM
    data = pandas.read_table(ersa_path, header=0,
                             names=['id1', 'id2', 'rel_est1', 'rel_est2',
                                    'ersa_degree', 'ersa_lower_bound', 'ersa_upper_bound', 'seg_count', 'total_seg_len'],
                             dtype={'ersa_degree': str, 'ersa_lower_bound': str, 'ersa_upper_bound': str,
                                    'seg_count': str, 'total_seg_len': str})

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
    data.loc[:, 'seg_count'] = pandas.to_numeric(data.seg_count.str.strip(), errors='coerce').astype(
        pandas.Int32Dtype())
    data.loc[:, 'total_seg_len'] = pandas.to_numeric(data.total_seg_len.str.replace(',', '').str.strip(),
                                                     errors='coerce').astype(float)
    print(f'parsed total_seg_len, found {(data.total_seg_len > 1000).sum()} long entries')
    data.loc[data.ersa_lower_bound != 1, 'ersa_lower_bound'] = data.ersa_lower_bound - 1
    data.loc[data.ersa_upper_bound != 1, 'ersa_upper_bound'] = data.ersa_upper_bound - 1

    logging.info(f'read {data.shape[0]} pairs from ersa output')
    logging.info(f'ersa after reading has {pandas.notna(data.ersa_degree).sum()}')
    return data.loc[data.id1 != data.id2, ['id1', 'id2', 'ersa_degree', 'ersa_lower_bound', 'ersa_upper_bound',
                                           'seg_count', 'total_seg_len']].\
        set_index(['id1', 'id2'])


if __name__ == '__main__':

    logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG, format='%(levelname)s:%(asctime)s %(message)s')

    ibd_path = snakemake.input['ibd']
    # within families
    # across families
    ersa_path = snakemake.input['ersa']
    output_path = snakemake.output[0]

    ibd = read_ibis(ibd_path)
    ersa = read_ersa(ersa_path)

    logging.info(f'ibd shape: {ibd.shape[0]}, ersa shape: {ersa.shape[0]}')

    relatives = ibd.merge(ersa, how='outer', left_index=True, right_index=True)

    # It is impossible to have more than 50% of ibd2 segments unless you are monozygotic twins or duplicates.
    GENOME_CM_LEN = 3440
    DUPLICATES_THRESHOLD = 1750
    FS_TOTAL_THRESHOLD = 2100
    FS_IBD2_THRESHOLD = 450
    dup_mask = relatives.total_seg_len_ibd2 > DUPLICATES_THRESHOLD
    fs_mask = (relatives.total_seg_len > FS_TOTAL_THRESHOLD) & (relatives.total_seg_len_ibd2 > FS_IBD2_THRESHOLD) & (~dup_mask)
    po_mask = (relatives.total_seg_len / GENOME_CM_LEN > 0.8) & (~fs_mask) & (~dup_mask)
    print(f'We have {sum(dup_mask)} duplicates, {sum(fs_mask)} full siblings and {sum(po_mask)} parent-offspring relationships')

    relatives.loc[:, 'final_degree'] = relatives.ersa_degree
    relatives.loc[dup_mask, 'final_degree'] = 0
    relatives.loc[po_mask, 'final_degree'] = 1

    relatives.loc[(~po_mask) & (relatives.final_degree == 1), 'final_degree'] = 2  # to eliminate some ersa false positives
    relatives.loc[fs_mask, 'final_degree'] = 2
    relatives.loc[:, 'relation'] = relatives.ersa_degree.astype(str)
    relatives.loc[po_mask, 'relation'] = 'PO'
    relatives.loc[fs_mask, 'relation'] = 'FS'
    relatives.loc[dup_mask, 'relation'] = 'MZ/Dup'

    # approximate calculations, IBD share is really small in this case
    relatives.loc[pandas.isna(relatives.total_seg_len_ibd2), 'total_seg_len_ibd2'] = 0
    relatives.loc[pandas.isna(relatives.seg_count_ibd2), 'seg_count_ibd2'] = 0
    '''
        IBIS and ERSA do not distinguish IBD1 and IBD2 segments
        shared_genome_proportion is a proportion of identical alleles
        i.e. shared_genome_proportion of [(0/0), (0/1)] and [(0/0), (1/1)] should be 3/4
        GENOME_SEG_LEN is length of centimorgans of the haplotype
        ibd2 segments span two haplotypes, therefore their length should be doubled
        total_seg_len includes both ibd1 and ibd2 segments length
    '''
    relatives.loc[:, 'shared_genome_proportion'] = (relatives.total_seg_len - relatives.total_seg_len_ibd2) / (2*GENOME_CM_LEN) + relatives.total_seg_len_ibd2*2 / (2*GENOME_CM_LEN)
    relatives.loc[relatives.shared_genome_proportion > 1.0, 'shared_genome_proportion'] = 1.0
    logging.info(f'final degree not null: {pandas.notna(relatives.final_degree).sum()}')
    relatives.loc[pandas.notna(relatives.final_degree), :].to_csv(output_path, sep='\t')
