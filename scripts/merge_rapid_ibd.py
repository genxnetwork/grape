from utils.ibd import merge_rapid_ibd_data


if __name__ == '__main__':
    files = list(snakemake.input)
    output_file = snakemake.output['ibd']
    segments_count = merge_rapid_ibd_data(files, output_file)
    print(f'Merged {segments_count} ibd segments from all chromosomes')
