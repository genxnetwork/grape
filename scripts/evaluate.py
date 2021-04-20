# many to many tests of relationships

import itertools
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import logging
import math
import networkx as nx
import pydot
from networkx.drawing.nx_pydot import graphviz_layout

from utils.relatives import get_kinship, read_pedigree, relatives_to_graph


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
        logging.info(f'td: {true_degree}\t tp: {true_positives}\t fn: {false_negatives}\t fp: {false_positives}')

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
        logging.info(f'plot saved to {plot_name}')
        plt.savefig(plot_name)
        plt.close()


def compare(total, correct, plot_name=None, cutoff=1):
    ds = []
    cs = []
    for i in sorted(total):
        if i >= cutoff:
            c = correct.get(i, 0)
            cs.append((c / total[i]) * 100)
            ds.append(i)
    df = pd.DataFrame(columns=['Degree of Relatedness', '% of true degree in predicted interval'])
    df.iloc[:, 0] = ds
    df.iloc[:, 1] = cs
    ax = sns.barplot(x="Degree of Relatedness", y="% of true degree in predicted interval", data=df, palette='muted')
    for _, row in df.iterrows():
        ax.text(row.name, row.iloc[1] + 1, round(row.iloc[1], 1), color='black', ha="center")
    ax.set(ylim=(0, 100))
    if not plot_name:
        plt.show()
    else:
        logging.info(f'plot saved to {plot_name}')
        plt.savefig(plot_name)
        plt.close()


def plot_confusion_matrix(confusion_matrix, plot_name):
    x = []
    y = []
    sizes = []
    norel_index = None
    index = 0
    max_degree = max(max([c[0] for c in confusion_matrix.keys()]), max([c[1] for c in confusion_matrix.keys()]))
    for (true_degree, pred_degree), count in confusion_matrix.items():

        if true_degree == -1:
            true_degree = max_degree + 1
        if pred_degree == -1:
            pred_degree = max_degree + 1

        x.append(true_degree)
        y.append(pred_degree)
        if true_degree == max_degree + 1 and pred_degree == max_degree + 1:
            norel_index = index
            sizes.append(0)
        else:
            sizes.append(count)
        index += 1

    normed_sizes = []
    max_count = max(sizes)
    for i, size in enumerate(sizes):
        if i == norel_index:
            normed_sizes.append(2500)
        else:
            normed_sizes.append(2000 * size / max_count)
    sizes[norel_index] = confusion_matrix[-1, -1]
    plt.figure(figsize=(15, 15))
    plt.ylim(0, max_degree + 2)
    plt.xlim(0, max_degree + 2)
    plt.xticks(ticks=list(range(max_degree + 2)),
               labels=[str(d) for d in range(max_degree + 1)] + ['No rel'], fontsize=20)
    plt.yticks(ticks=list(range(max_degree + 2)),
               labels=[str(d) for d in range(max_degree + 1)] + ['No rel'], fontsize=20)
    plt.grid(zorder=0)
    plt.scatter(x, y, s=normed_sizes, zorder=3)
    for i, size in enumerate(sizes):
        plt.annotate(size, (x[i], y[i]), xytext=(x[i]+0.17, y[i]+0.29), fontsize=14)
    plt.xlabel('Ground truth', fontsize=24)
    plt.ylabel('Predicted degree', fontsize=24)
    plt.savefig(plot_name)
    plt.close()


def get_pred_degree(edge, source):
    king = edge['king']
    ersa = edge['ersa']
    final = edge['final']
    degree = {
        'both': final,
        'king': king,
        'ersa': ersa
    }[source]
    return int(degree) if not math.isnan(degree) else -1


def kinship_to_dataframe(kinship: nx.Graph):
    records = []
    for id1, id2, degree in kinship.edges.data("degree", default=-1):
        record = (id1, id2, degree) if id1 < id2 else (id2, id1, degree)
        records.append(record)
    data = pd.DataFrame.from_records(records, columns=['id1', 'id2', 'true_degree'])
    return data.set_index(['id1', 'id2'])


def draw_pedigree(pedigree: nx.DiGraph, pedigree_plot_name: str):

    plt.figure(figsize=(20, 20))
    first1 = pedigree.subgraph([n for n in pedigree.nodes if n.startswith('first1')])
    pos = graphviz_layout(first1, prog="dot")
    nx.draw(first1, pos, with_labels=True)
    plt.savefig(pedigree_plot_name)
    plt.close()


def interval_precision_recall(kinship, inferred, clients, source, plot_name):
    iterator = itertools.combinations(clients, 2)

    true_positives = {degree: 0 for degree in range(1, 15)}
    false_negatives = {degree: 0 for degree in range(1, 15)}
    false_positives = {degree: 0 for degree in range(1, 15)}

    for i, j in iterator:
        if kinship.has_edge(i, j):
            degree = kinship[i][j]['degree']
        else:
            degree = -1

        if inferred.has_edge(i, j):
            pred_degree = get_pred_degree(inferred[i][j], source)
        else:
            pred_degree = -1

        if inferred.has_edge(i, j):
            lower, upper = inferred[i][j].get('ersa_lower_bound', pred_degree), inferred[i][j].get('ersa_upper_bound', pred_degree)
        else:
            lower, upper = -1, -1

        if degree != -1 and pred_degree != -1:
            if degree <= 1 and pred_degree == degree:
                true_positives[degree] += 1
            elif 1 < degree <= 4 and degree - 1 <= pred_degree <= degree + 1:
                true_positives[degree] += 1
            elif degree > 4 and lower <= degree <= upper:
                true_positives[degree] += 1
            else:
                false_negatives[degree] += 1
                false_positives[pred_degree] += 1

        elif degree != -1 and pred_degree == -1:
            false_negatives[degree] += 1

    data = {'True Degree': [], 'Precision': [], 'Recall': []}
    for degree in range(1, 14):
        tp, fp, fn = true_positives[degree], false_positives[degree], false_negatives[degree]
        if tp + fp > 0:
            precision = tp / (tp + fp)
        else:
            precision = 0.0
        if tp + fn > 0:
            recall = tp / (tp + fn)
        else:
            recall = 0.0
        data['True Degree'].append(degree)
        data['Precision'].append(precision)
        data['Recall'].append(recall)

    df = pd.DataFrame.from_dict(data)
    df.set_index('True Degree').plot.bar()
    logging.info('Precision-Recall DataFrame:')
    logging.info(df)
    if not plot_name:
        plt.show()
    else:
        logging.info(f'plot saved to {plot_name}')
        plt.savefig(plot_name)
        plt.close()


def interval_evaluate(result, fam, plot_name, pr_plot_name, conf_matrix_plot_name, output_path, only_client=False, source='both', pedigree_plot_name=None, dist_plot_name=None):
    iids, pedigree = read_pedigree(fn=fam)
    if pedigree_plot_name is not None:
        draw_pedigree(pedigree, pedigree_plot_name)
    kinship = get_kinship(pedigree)
    inferred, clients = relatives_to_graph(result, only_client)
    confusion_matrix = {}
    total = {}
    correct = {}
    if only_client:
        iterator = itertools.combinations(clients, 2)
    else:
        iterator = itertools.product(clients, iids - clients)

    for i, j in iterator:
        if kinship.has_edge(i, j):
            degree = kinship[i][j]['degree']
            total[degree] = total.get(degree, 0) + 1
        else:
            degree = -1

        if inferred.has_edge(i, j):
            pred_degree = get_pred_degree(inferred[i][j], source)
        else:
            pred_degree = -1

        key = (degree, pred_degree)
        confusion_matrix[key] = confusion_matrix.get(key, 0) + 1
        if inferred.has_edge(i, j):
            lower, upper = inferred[i][j].get('ersa_lower_bound', pred_degree), inferred[i][j].get('ersa_upper_bound', pred_degree)
        else:
            lower, upper = -1, -1
        if degree != -1 and pred_degree != -1:
            if degree <= 1 and pred_degree == degree:
                correct[degree] = correct.get(degree, 0) + 1
            elif 1 < degree <= 4 and degree - 1 <= pred_degree <= degree + 1:
                correct[degree] = correct.get(degree, 0) + 1
            elif degree > 4 and lower <= degree <= upper:
                correct[degree] = correct.get(degree, 0) + 1

    print(f'Evaluation results for {source}:')
    logging.info(f'Evaluation results for {source}')
    if not total:
        print('total is not total')
        return
    print('correct: ', correct)
    print('total: ', total)
    keys = sorted(list(confusion_matrix.keys()))
    for key in keys:
        print(f'{key}\t{confusion_matrix[key]}')
        logging.info(f'{key}\t{confusion_matrix[key]}')
    print()

    compare(total, correct, plot_name)
    interval_precision_recall(kinship, inferred, clients, source, pr_plot_name)

    plot_confusion_matrix(confusion_matrix, conf_matrix_plot_name)

    kinship_frame = kinship_to_dataframe(kinship)
    predictions = pd.read_table(result, index_col=['id1', 'id2'])

    merged = predictions.merge(kinship_frame, left_index=True, right_index=True, how='outer')
    if dist_plot_name is not None:
        draw_distribution_plot(merged, dist_plot_name)
    merged.to_csv(output_path, sep='\t')


def evaluate(result, fam, plot_name, pr_plot_name, conf_matrix_plot_name, output_path, only_client=False, source='both', pedigree_plot_name=None, dist_plot_name=None):
    """If only_client=True, all pairwise relatives between client are evaluated"""
    iids, pedigree = read_pedigree(fn=fam)
    if pedigree_plot_name is not None:
        draw_pedigree(pedigree, pedigree_plot_name)
    kinship = get_kinship(pedigree)
    inferred, clients = relatives_to_graph(result, only_client)
    confusion_matrix = {}
    total = {}
    correct = {}
    if only_client:
        iterator = itertools.combinations(clients, 2)
    else:
        iterator = itertools.product(clients, iids - clients)

    for i, j in iterator:
        if kinship.has_edge(i, j):
            degree = kinship[i][j]['degree']
            total[degree] = total.get(degree, 0) + 1
        else:
            degree = -1

        if inferred.has_edge(i, j):
            pred_degree = get_pred_degree(inferred[i][j], source)
        else:
            pred_degree = -1

        key = (degree, pred_degree)
        confusion_matrix[key] = confusion_matrix.get(key, 0) + 1
        if degree != -1 and pred_degree != -1 and degree - 1 <= pred_degree <= degree + 1:
            correct[degree] = correct.get(degree, 0) + 1

    if not total:
        print('total is not total')
        return
    print('correct: ', correct)
    print('total: ', total)
    keys = sorted(list(confusion_matrix.keys()))
    for key in keys:
        print(f'{key}\t{confusion_matrix[key]}')
        logging.info(f'{key}\t{confusion_matrix[key]}')

    compare(total, correct, plot_name)
    precision_recall(total, confusion_matrix, pr_plot_name)
    plot_confusion_matrix(confusion_matrix, conf_matrix_plot_name)

    kinship_frame = kinship_to_dataframe(kinship)
    predictions = pd.read_table(result, index_col=['id1', 'id2'])

    merged = predictions.merge(kinship_frame, left_index=True, right_index=True, how='outer')
    if dist_plot_name is not None:
        draw_distribution_plot(merged, dist_plot_name)
    merged.to_csv(output_path, sep='\t')


def is_aunt(id1, id2, degree):
    if degree != 3:
        return False

    iid1 = id1.split('_')[1]
    iid2 = id2.split('_')[1]

    g1, b1, _ = iid1.split('-')
    g2, b2, _ = iid2.split('-')
    # print(g1, b1, g2, b2)
    if g1 == 'g2' and g2 == 'g3' and b1 != b2:
        return True
    return False


def is_grandmother(id1, id2, degree):
    if degree != 2:
        return False
    iid1 = id1.split('_')[1]
    iid2 = id2.split('_')[1]

    g1, b1, _ = iid1.split('-')
    g2, b2, _ = iid2.split('-')
    # print(g1, b1, g2, b2)
    if g1 != g2:
        return True
    return False


def draw_distribution_plot(merged, dist_plot_name):
    m = merged.reset_index()
    aunt_mask = [is_aunt(id1, id2, degree) for id1, id2, degree in zip(m.id1, m.id2, m.true_degree)]

    aunt_seg_len = m.loc[aunt_mask, 'total_seg_len']
    aunt_seg_count = m.loc[aunt_mask, 'seg_count']

    grandmother_mask = [is_grandmother(id1, id2, degree) for id1, id2, degree in zip(m.id1, m.id2, m.true_degree)]
    gm_seg_len = m.loc[grandmother_mask, 'total_seg_len']
    gm_seg_count = m.loc[grandmother_mask, 'seg_count']

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(17, 12))

    ax1.hist([aunt_seg_len, gm_seg_len], label=['Aunt/Niece', 'Grandmother/Granddaughter'])
    ax2.hist([aunt_seg_count, gm_seg_count], label=['Aunt/Niece', 'Grandmother/granddaughter'])
    ax1.legend()
    ax2.legend()
    ax1.set_title('Segments length distribution')
    ax2.set_title('Segments count distribution')
    plt.savefig(dist_plot_name)
    plt.close()


if __name__ == '__main__':
    logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG, format='%(levelname)s:%(asctime)s %(message)s')
    sns.set_style()

    interval_evaluate(snakemake.input['rel'],
             snakemake.input['fam'],
             snakemake.output['accuracy'],
             snakemake.output['pr'],
             snakemake.output['conf_matrix'],
             snakemake.output['updated_rel'],
             only_client=True,
             source='both',
             pedigree_plot_name=snakemake.output['pedigree_plot'],
             dist_plot_name=snakemake.output['dist_plot'])

    interval_evaluate(snakemake.input['rel'],
             snakemake.input['fam'],
             snakemake.output['king_accuracy'],
             snakemake.output['king_pr'],
             snakemake.output['king_conf_matrix'],
             snakemake.output['king_updated_rel'],
             only_client=True,
             source='king')

    interval_evaluate(snakemake.input['rel'],
             snakemake.input['fam'],
             snakemake.output['ersa_accuracy'],
             snakemake.output['ersa_pr'],
             snakemake.output['ersa_conf_matrix'],
             snakemake.output['ersa_updated_rel'],
             only_client=True,
             source='ersa')
