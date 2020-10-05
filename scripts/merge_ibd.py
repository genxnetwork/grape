from ibd import merge_germline_segments, segments_to_germline, interpolate_all
import pandas


if __name__ == '__main__':
    germline_path = snakemake.input['germline']
    map_dir = snakemake.params['cm_dir']
    gap = float(snakemake.params['merge_gap'])

    segments = read_germline_segments(germline_path)
    segments = interpolate_all(segments, map_dir)
    segments = merge_germline_segments(segments, gap=gap)

    frame = pandas.DataFrame.from_dict(segments, orient='index')
    frame.to_csv(snakemake.output['overlap'], sep='\t')