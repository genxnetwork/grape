# many to many tests of relationships

import itertools
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from utils.relatives import get_kinship, read_pedigree, read_pipeline_output


def precision_recall(total, confusion_matrix, plot_name=None):

    data = {'True Degree': [], 'Precision': [], 'Recall': []}

    max_degree = max(c[0] for c in confusion_matrix.keys())

    for true_degree in range(1, max_degree + 1):

        # precision: how many of the selected items are relevant
        # recall: how many of the relevant items are selected
        true_positives, false_negatives, false_positives = 0, 0, 0
        for key, value in confusion_matrix.items():
            predicted_degree = key[1]
            if true_degree == key[0]:
                if true_degree - 1 <= predicted_degree <= true_degree + 1:
                    true_positives += value
                #else:
                #    false_negatives += value
            if key[0] == -1 and true_degree == predicted_degree:
                false_positives += value

        true_total = total[true_degree]
        false_negatives = true_total - true_positives

        print(f'td: {true_degree}\t tp: {true_positives}\t fn: {false_negatives}\t fp: {false_positives}')
        if true_positives + false_positives == 0:
            data['Precision'].append(0.0)
        else:
            data['Precision'].append(true_positives / (true_positives + false_positives))
        if true_positives + false_negatives == 0:
            data['Recall'].append(0.0)
        else:
            data['Recall'].append(true_positives / (true_positives + false_negatives))
        data['True Degree'].append(true_degree)

    df = pd.DataFrame.from_dict(data)
    df.set_index('True Degree').plot.bar()

    if not plot_name:
        plt.show()
    else:
        print('plot saved to ', plot_name)
        plt.savefig(plot_name)
        plt.close()
    print('precision and recall plotting finished')


def compare(total, correct, plot_name=None, cutoff=1):
    ds = []
    cs = []
    for i in sorted(total):
        if i >= cutoff:
            c = correct.get(i, 0)
            cs.append((c / total[i]) * 100)
            ds.append(i)
    df = pd.DataFrame(columns=['Degree of Relatedness', '% predicted correct +/- 1 degree relationship'])
    df.iloc[:, 0] = ds
    df.iloc[:, 1] = cs
    ax = sns.barplot(x="Degree of Relatedness", y="% predicted correct +/- 1 degree relationship", data=df,
                     palette='muted')
    for _, row in df.iterrows():
        ax.text(row.name, row.iloc[1] + 1, round(row.iloc[1], 1), color='black', ha="center")
    ax.set(ylim=(0, 100))
    if not plot_name:
        plt.show()
    else:
        print('plot saved to ', plot_name)
        plt.savefig(plot_name)
        plt.close()
    print('comparison finished')


def evaluate(result, fam, plot_name, pr_plot_name, only_client=False):
    """If only_client=True, all pairwise relatives between client are evaluated"""
    iids, pedigree = read_pedigree(fn=fam)
    print('pedigree:')
    print(list(pedigree.edges)[:10])
    kinship = get_kinship(pedigree)
    print('results is: ', result)
    inferred, clients = read_pipeline_output(result, only_client)
    print('inferred:')
    print(len(inferred), list(inferred.edges)[:10])
    print('*'*100)
    print('kinship:')
    print(len(kinship), list(kinship.edges)[:10])
    confusion_matrix = {}
    total = {}
    correct = {}
    if only_client:
        iterator = itertools.combinations(clients, 2)
    else:
        iterator = itertools.product(clients, iids - clients)
    for i, j in iterator:
        # check if name is this format 'fid_iid'
        #if '_' in i:
        #    i = "_".join(i.split('_')[1:])
        #if '_' in j:
        #    j = "_".join(j.split('_')[1:])
        if kinship.has_edge(i, j):
            degree = kinship[i][j]['degree']
            total[degree] = total.get(degree, 0) + 1
            if inferred.has_edge(i, j) and inferred[i][j]['ersa'] != 'NA':
                ersa_d = int(inferred[i][j]['ersa'])
                #print(i, j, degree, ersa_d, inferred[i][j]['king'])
                key = (degree, ersa_d)
                confusion_matrix[key] = confusion_matrix.get(key, 0) + 1
                if degree - 1 <= ersa_d <= degree + 1:
                    correct[degree] = correct.get(degree, 0) + 1

        else:
            if inferred.has_edge(i, j) and inferred[i][j]['ersa'] != 'NA':
                ersa_d = int(inferred[i][j]['ersa'])
                #print(i, j, degree, ersa_d, inferred[i][j]['king'])
                degree = -1
                key = (degree, ersa_d)
                confusion_matrix[key] = confusion_matrix.get(key, 0) + 1

    if not total:
        print('total is not total')
        return
    print('correct: ', correct)
    print('total: ', total)
    keys = sorted(list(confusion_matrix.keys()))
    for key in keys:
        print(f'{key}\t{confusion_matrix[key]}')

    compare(total, correct, plot_name)
    precision_recall(total, confusion_matrix, pr_plot_name)


if __name__ == '__main__':
    evaluate(snakemake.input['rel'], snakemake.input['fam'], snakemake.output['accuracy'], snakemake.output['pr'], only_client=True)
