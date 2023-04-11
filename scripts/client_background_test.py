import pandas


# find rows in one dataframe that are not in the other
def rows_in_one_dataframe(first: pandas.DataFrame, second: pandas.DataFrame):
    second = second.copy()
    second.loc[:, 'indicator'] = 'second'
    left = first.merge(second, on=['id1', 'id2'], how='left', suffixes=('_first', '_second'))
    return left[pandas.isnull(left.indicator)]


def reorder_ids(matches: pandas.DataFrame):
    matches.loc[matches.id1 < matches.id2, ['id1', 'id2']] = matches.loc[matches.id1 < matches.id2, ['id2', 'id1']].values
    
# load client-background matches and compare it with the local matches
def compare_cb_local_matches(cb_matches_path: str, local_matches_path: str):
    cb_matches = pandas.read_table(cb_matches_path)
    matches = pandas.read_table(local_matches_path)
    
    cb_matches.loc[:, 'id2'] = cb_matches.id2.str.split('@').str[2]
    reorder_ids(cb_matches)
    reorder_ids(matches)
    common = cb_matches.merge(matches, on=['id1', 'id2'], how='inner', suffixes=('_cb', '_local'))
    print(f'common: {common.shape[0]}\tcb: {cb_matches.shape[0]}\tlocal: {matches.shape[0]}')
    
    left_only = rows_in_one_dataframe(cb_matches, matches)
    right_only = rows_in_one_dataframe(matches, cb_matches)
    
    left_only.to_csv('workdir/left_only.tsv', sep='\t', index=False)
    right_only.to_csv('workdir/right_only.tsv', sep='\t', index=False)
    
    assert common.shape[0] == cb_matches.shape[0]
    assert common.shape[0] == matches.shape[0]
    
    assert common[common.final_degree_cb == common.final_degree_local].all()


# extract client-background matches from the local matches
def extract_cb_matches_from_local(input150_path: str, input18_path: str, local_matches_path: str, output_path: str):
    input150_samples = pandas.read_csv(input150_path, sep='\t', header=None)
    input18_samples = pandas.read_csv(input18_path, sep='\t', header=None)
    
    relatives = pandas.read_table(local_matches_path)
    
    matches18_150 = relatives[relatives.id1.isin(input18_samples[0]) & relatives.id2.isin(input150_samples[0])]
    matches150_18 = relatives[relatives.id1.isin(input150_samples[0]) & relatives.id2.isin(input18_samples[0])]
    
    matches = pandas.concat([matches18_150, matches150_18])
    matches.to_csv(output_path, sep='\t', index=False)
    

if __name__ == '__main__':
    input150_samples_path = 'workdir/input150.samples'    
    input18_samples_path = 'workdir/input18.samples'
    local_matches_path = 'workdir/relatives_ibis168.tsv'
    cb_matches_path = 'workdir/relatives_cb_input18_150.tsv'
    
    local_cb_matches_path = 'workdir/relatives_18_150_2.tsv'
    
    extract_cb_matches_from_local(input150_samples_path, input18_samples_path, local_matches_path, local_cb_matches_path)
    compare_cb_local_matches(cb_matches_path, local_cb_matches_path)
    