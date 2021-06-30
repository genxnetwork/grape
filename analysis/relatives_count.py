import matplotlib.pyplot as plt
import pandas
from collections import Counter
import numpy


def get_data(relatives_path):
    relatives = pandas.read_table(relatives_path, index_col=['id1', 'id2'])
    counter = Counter(relatives.final_degree.values)

    data = []
    for i in range(1, 11):
        data.append(counter[i])

    return data


if __name__ == '__main__':
    relatives_hard = 'data/alpha001_zsl6_zsc05/relatives.tsv'
    # relatives_soft = 'data/alpha005_zsl5_zsc01/relatives.tsv'
    relatives_soft = '/media_ssd/pipeline_data/big_ibis_alpha005/results/relatives.tsv'
    image_path = 'data/hard_soft_params.png'

    plt.figure(figsize=(15, 10))
    plt.xlabel('Predicted degree')
    plt.ylabel('Number of predicted relatives')
    plt.grid(True)

    hard = get_data(relatives_hard)
    soft = get_data(relatives_soft)

    plt.bar(numpy.arange(1, 11) - 0.22, hard, width=0.35, alpha=0.8, label='alpha=0.01, zero_seg_len=6, zero_seg_count=0.5')
    plt.bar(numpy.arange(1, 11) + 0.22, soft, width=0.35, alpha=0.8, label='alpha=0.05, zero_seg_len=5, zero_seg_count=0.1')
    print(hard)
    print(soft)
    plt.xticks(numpy.arange(1, 11), fontsize=18)
    plt.legend(fontsize=18)
    plt.savefig(image_path)
