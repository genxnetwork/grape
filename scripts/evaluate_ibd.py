from ibd import read_pedsim_segments, read_germline_segments, interpolate_all, total_overlap
import pandas


if __name__ == '__main__':
    pedsim_path = snakemake.input['pedsim']
    germline_path = snakemake.input['germline']
    map_dir = snakemake.params['cm_dir']
    true_segments = read_pedsim_segments(pedsim_path)
    germline_segments = read_germline_segments(germline_path)

    germline_segments = interpolate_all(germline_segments, map_dir)

    overlaps = total_overlap(true_segments, germline_segments)

    frame = pandas.DataFrame.from_dict(overlaps, orient='index')

    frame.to_csv(snakemake.output['overlap'], sep='\t')
