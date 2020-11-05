from utils.ibd import read_pedsim_segments, read_germline_segments, interpolate_all, total_overlap, read_rapid_segments
import pandas
import matplotlib.pyplot as plt
import seaborn as sns
import numpy


def plot_segment_accuracy(true_segments, found_segments, plot_name):

    accuracy = {'<2.5': [], '<10': [], '<100': [], '>100': []}
    # key is (fid_iid1, fid_iid2), sorted
    for key, true_segs in true_segments.items():
        if key[0] == key[1]:
            continue
        if key in found_segments:
            found_segs = found_segments[key]
            for ts in true_segs:
                overlap = 0.0
                for fs in found_segs:
                    overlap += ts.overlap(fs)

                if ts.cm_len <= 2.5:
                    accuracy['<2.5'].append(overlap / ts.cm_len)
                elif ts.cm_len <= 10:
                    accuracy['<10'].append(overlap / ts.cm_len)
                elif ts.cm_len <= 100:
                    accuracy['<100'].append(overlap / ts.cm_len)
                else:
                    accuracy['>100'].append(overlap / ts.cm_len)

    records = {key: numpy.mean(value) for key, value in accuracy.items()}
    df = pandas.DataFrame.from_dict(records, orient='index', columns=['Accuracy'])
    df.reset_index(level=0, inplace=True, )
    sns.barplot(x='index', y='Accuracy', data=df, palette='muted')
    if not plot_name:
        plt.show()
    else:
        print('plot saved to ', plot_name)
        plt.savefig(plot_name)
        plt.close()


def plot_ibd_len_accuracy(true_segments, found_segments, plot_name):

    accuracy = {'<10': [], '<100': [], '<250': [], '>250': []}
    # key is (fid_iid1, fid_iid2), sorted
    for key, true_segs in true_segments.items():

        if key in found_segments:
            found_segs = found_segments[key]
            found_len = sum([seg.cm_len for seg in found_segs])
        else:
            found_len = 0.0

        true_len = sum([seg.cm_len for seg in true_segs])

        if true_len <= 1.0:
            continue
        if true_len <= 10:
            accuracy['<10'].append(found_len/true_len)
        elif true_len <= 100:
            accuracy['<100'].append(found_len/true_len)
        elif true_len <= 250:
            accuracy['<250'].append(found_len/true_len)
        else:
            accuracy['>250'].append(found_len/true_len)

    print('counts of data')
    print({key: len(value) for key, value in accuracy.items()})
    records = {key: numpy.mean(value) for key, value in accuracy.items()}
    df = pandas.DataFrame.from_dict(records, orient='index', columns=['Accuracy'])
    df.reset_index(level=0, inplace=True, )
    sns.barplot(x='index', y='Accuracy', data=df, palette='muted')
    if not plot_name:
        plt.show()
    else:
        print('plot saved to ', plot_name)
        plt.savefig(plot_name)
        plt.close()


if __name__ == '__main__':
    pedsim_path = snakemake.input['pedsim']
    ibd_path = snakemake.input['ibd']
    map_dir = snakemake.params['cm_dir']
    is_rapid_ibd = bool(snakemake.params['is_rapid_ibd'])

    true_segments = read_pedsim_segments(pedsim_path)
    if not is_rapid_ibd:
        found_segments = read_germline_segments(ibd_path)
    else:
        found_segments = read_rapid_segments(ibd_path)

    found_segments = interpolate_all(found_segments, map_dir)

    overlaps = total_overlap(true_segments, found_segments)

    frame = pandas.DataFrame.from_dict(overlaps, orient='index')

    frame.to_csv(snakemake.output['overlap'], sep='\t')

    plot_name = snakemake.output['seg_accuracy']
    total_acc_name = snakemake.output['total_len_accuracy']
    plot_segment_accuracy(true_segments, found_segments, plot_name)
    plot_ibd_len_accuracy(true_segments, found_segments, total_acc_name)
