from utils.ibd import merge_germline_segments, segments_to_germline, interpolate_all, read_rapid_segments, read_pedsim_segments


if __name__ == '__main__':

    ibd_path = snakemake.input['ibd']
    map_dir = snakemake.params['cm_dir']
    gap = float(snakemake.params['merge_gap'])
    pedsim_path = snakemake.input['true_ibd']

    use_true_ibd = bool(snakemake.params['use_true_ibd'])
    if not use_true_ibd:
        segments = read_rapid_segments(ibd_path)
        segments = interpolate_all(segments, map_dir)
        segments = merge_germline_segments(segments, gap=gap)
        segments_to_germline(segments, snakemake.output['ibd'])
    else:
        true_segments = read_pedsim_segments(pedsim_path)
        segments_to_germline(true_segments, snakemake.output['ibd'])
