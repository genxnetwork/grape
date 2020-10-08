from utils.ibd import read_pedsim_segments, read_germline_segments, interpolate_all, total_overlap
import pandas
import matplotlib.pyplot as plt
import seaborn as sns


def plot_segment_accuracy(true_segments, found_segments, plot_name):

    accuracy = {'<10': [], '<100': [], '>100': []}

    for key, true_segs in true_segments.items():
        overlap = 0.0
        true_len = 0.0

        if key in found_segments:
            found_segs = found_segments[key]
            for ts in true_segs:
                for fs in found_segs:
                    overlap += ts.overlap(fs)

        true_len += sum([seg.cm_len for seg in true_segs])

        if true_len <= 10:
            accuracy['<10'].append(overlap/true_len)
        elif true_len <= 100:
            accuracy['<100'].append(overlap / true_len)
        else:
            accuracy['>100'].append(overlap / true_len)

    df = pandas.DataFrame.from_dict(accuracy).mean().transpose()
    df.columns = ['Accuracy']
    sns.barplot(y='Accuracy', data=df, palette='muted')
    if not plot_name:
        plt.show()
    else:
        print('plot saved to ', plot_name)
        plt.savefig(plot_name)
        plt.close()


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

    #plot_segment_accuracy(true_segments, germline_segments)
