import pandas


if __name__ == '__main__':
    bim_path = 'test_data/sim-22-ibis-small/merged_ibis.bim'
    output_path = 'test_data/sim-22-ibis-small/bad_regions.tsv'
    bim = pandas.read_table(bim_path, names=['chr', 'rs', 'zero', 'pos', 'ref', 'alt'])
    bad_regions = []
    window_size = int(1e+6)  # 1 MB as in Gusev et al 2012
    print(bim.columns)
    for chr, group in bim.groupby(by='chr'):
        for left in range(0, group.pos.max(), window_size):
            if ((group.pos >= left) & (group.pos <= left + window_size)).sum() < 100:
                region = (chr, left, left + window_size)
                if bad_regions and bad_regions[-1][2] == left:
                    bad_regions[-1] = ((chr, bad_regions[-1][1], left + window_size))
                else:
                    bad_regions.append(region)

    data = pandas.DataFrame.from_records(bad_regions, columns=['chr', 'start_mb', 'end_mb'])
    data.to_csv(output_path, index=False, sep='\t')
