import pandas
import numpy
import os


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
    segments_info.rename({'genetic_len': 'total_seg_len', 'chrom': 'seg_count'}, axis=1, inplace=True)
    return segments_info


def map_king_degree(king_degree):
    degree_map = {
        'PO': 1,
        'FS': 2,
        '2nd': 2,
        '3rd': 3,
        '4th': 4,
        'NA': numpy.nan
    }
    return [degree_map[kd] for kd in king_degree]


def read_king(king_path):
    # FID1    ID1     FID2    ID2     MaxIBD1 MaxIBD2 IBD1Seg IBD2Seg PropIBD InfType
    data = pandas.read_table(king_path)
    data.loc[:, 'id1'] = data.FID1 + '_' + data.ID1
    data.loc[:, 'id2'] = data.FID2 + '_' + data.ID2

    data.loc[:, 'king_degree'] = map_king_degree(data.InfType)
    data.loc[:, 'king_degree'] = data.king_degree.astype(float).astype(pandas.Int32Dtype())
    indexed = data.loc[:, ['id1', 'id2', 'king_degree']].set_index(['id1', 'id2'])
    return indexed


def read_kinship(kinship_path, kinship0_path):
    # parse within families
    # FID     ID1     ID2     N_SNP   Z0      Phi     HetHet  IBS0    Kinship Error
    within, across = None, None
    if is_non_zero_file(kinship_path):
        within = pandas.read_table(kinship_path)
        within.loc[:, 'id1'] = within.FID + '_' + within.ID1
        within.loc[:, 'id2'] = within.FID + '_' + within.ID2
        within.rename({'Kinship': 'kinship'}, axis=1, inplace=True)
        within = within.loc[:, ['id1', 'id2', 'kinship']].set_index(['id1', 'id2'])

    # FID1    ID1     FID2    ID2     N_SNP   HetHet  IBS0    Kinship
    if is_non_zero_file(kinship0_path):
        across = pandas.read_table(kinship0_path)
        across.loc[:, 'id1'] = across.FID1 + '_' + across.ID1
        across.loc[:, 'id2'] = across.FID2 + '_' + across.ID2
        across.rename({'Kinship': 'kinship'}, axis=1, inplace=True)
        across = across.loc[:, ['id1', 'id2', 'kinship']].set_index(['id1', 'id2'])

    if within is None and across is None:
        return None
    elif within is None and across is not None:
        return across
    elif within is not None and across is None:
        return within
    else:
        pandas.concat([within, across], axis=0)


def read_ersa(ersa_path):
    # Indv_1     Indv_2      Rel_est1      Rel_est2      d_est     N_seg     Tot_cM
    data = pandas.read_table(ersa_path, header=0,
                             names=['id1', 'id2', 'rel_est1', 'rel_est2', 'ersa_degree', 'seg_count', 'total_seg_len'])

    data = data.loc[(data.rel_est1 != 'NA') | (data.rel_est2 != 'NA'), :]
    data.loc[:, 'id1'] = data.id1.str.strip()
    data.loc[:, 'id2'] = data.id2.str.strip()
    data.loc[:, 'ersa_degree'] = pandas.to_numeric(data.ersa_degree.str.strip(), errors='coerce').astype(pandas.Int32Dtype())


    print(f'read {data.shape[0]} pairs from ersa output')
    #print(data.iloc[0, :])

    print(len(numpy.unique(data.id1)), len(numpy.unique(data.id2)))

    return data.loc[:, ['id1', 'id2', 'ersa_degree']].set_index(['id1', 'id2'])


if __name__ == '__main__':

    '''
    ibd_path = 'test_data/merge_king_ersa/merged_ibd.tsv'
    king_path = 'test_data/merge_king_ersa/merged_imputed_king.seg'
    # within families
    kinship_path = 'test_data/merge_king_ersa/merged_imputed_kinship.kin'
    # across families
    kinship0_path = 'test_data/merge_king_ersa/merged_imputed_kinship.kin0'
    ersa_path = 'test_data/merge_king_ersa/relatives.tsv'
    '''

    ibd_path = snakemake.input['ibd']
    king_path = snakemake.input['king']
    # within families
    kinship_path = snakemake.input['kinship']
    # across families
    kinship0_path = snakemake.input['kinship0']
    ersa_path = snakemake.input['ersa']

    ibd = read_germline(ibd_path)
    king = read_king(king_path)
    kinship = read_kinship(kinship_path, kinship0_path)
    ersa = read_ersa(ersa_path)

    if kinship is not None:
        relatives = ibd.merge(king, how='outer', left_index=True, right_index=True).\
            merge(kinship, how='outer', left_index=True, right_index=True).\
            merge(ersa, how='outer', left_index=True, right_index=True)
    else:
        relatives = ibd.merge(king, how='outer', left_index=True, right_index=True).\
            merge(ersa, how='outer', left_index=True, right_index=True)

    prefer_king_mask = pandas.isnull(relatives.king_degree) | (relatives.king_degree > 3)
    relatives.loc[:, 'final_degree'] = relatives.king_degree
    # if king is unsure or king degree > 3 then we use ersa distant relatives estimation
    relatives.loc[prefer_king_mask, 'final_degree'] = relatives.ersa_degree

    output_path = snakemake.output[0]

    relatives.loc[pandas.notna(relatives.final_degree), :].to_csv(output_path, sep='\t')
