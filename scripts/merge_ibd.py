from utils.ibd import merge_germline_segments, segments_to_germline, interpolate_all, read_germline_segments, read_pedsim_segments


if __name__ == '__main__':

    germline_path = snakemake.input['germline']
    map_dir = snakemake.params['cm_dir']
    gap = float(snakemake.params['merge_gap'])
    pedsim_path = snakemake.input['true_ibd']

    use_true_ibd = bool(snakemake.params['use_true_ibd'])
    if not use_true_ibd:
        segments = read_germline_segments(germline_path)
        segments = interpolate_all(segments, map_dir)
        segments = merge_germline_segments(segments, gap=gap)
        segments_to_germline(segments, snakemake.output['ibd'])
    else:
        true_segments = read_pedsim_segments(pedsim_path)
        segments_to_germline(true_segments, snakemake.output['ibd'])
