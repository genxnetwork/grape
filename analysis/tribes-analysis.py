import itertools
import pandas as pd
import seaborn as sns
import networkx as nx
import matplotlib.pyplot as plt
import math
from scripts.utils.relatives import relatives_to_graph, read_pedigree, get_kinship
import numpy


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
            if size > 0.25*max_count:
                normed_sizes.append(2000 * 0.25)
            else:
                normed_sizes.append(2000*size / max_count)

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
        elif degree == -1 and pred_degree != -1:
            false_positives[pred_degree] += 1

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
    if not plot_name:
        plt.show()
    else:
        plt.savefig(plot_name)
        plt.close()


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


def get_tribes_kinship(tribes_rel: str) -> nx.Graph:
    """
    Given a networkx object, get all the pairwise relationships

    Args:
        tribes_rel (str): path to tribes relations tsv file in form ID1	ID2	Degree	NoMeioses	NoAnc
    Returns:
        networkx.Graph: undirected graph with all pairwise relationships
    """
    # get all the descendants (0 length means self)

    tribes = pd.read_table(tribes_rel)
    relatives = nx.Graph()
    for id1, id2, degree in zip(tribes.ID1, tribes.ID2, tribes.Degree):
        print(degree, type(degree))
        if degree == 'NA':
            continue
        if degree == 'PO':
            degree = 1
        elif degree == 'FS':
            degree = 2
        elif not isinstance(degree, str) and math.isnan(degree):
            continue
        else:
            degree = int(degree)
        relatives.add_edge(id1, id2, common_ancestor='', degree=degree)

    return relatives


def map_tribes_degree(degree):
    mapped = []
    for d in degree:
        if d == 'PO':
            mapped.append(1)
        elif d == 'NA':
            mapped.append(numpy.nan)
        else:
            mapped.append(int(d))

    return mapped


def tribes_degree_to_interval(degree):
    interval_map = {
        1: (1, 1),
        2: (2, 3),
        3: (2, 4),
        4: (3, 5),
        5: (4, 7),
        6: (4, 8),
        7: (5, 10),
        8: (5, 11),
        9: (5, 13),
        10: (6, 14),
        11: (7, 14),
        12: (7, 14),
        13: (7, 14),
        14: (7, 14)
    }
    lower = [interval_map[d][0] for d in degree]
    upper = [interval_map[d][1] for d in degree]
    return lower, upper


def tribes_to_relatives(path, output_path):
    # Id1,Id2,IBD0.cM,IBD1.cM,IBD2.cM,EstDegree
    tribes = pd.read_csv(path)
    # id1 id2 king_degree king_relation shared_genome_proportion kinship kinship_degree ersa_degree ersa_lower_bound
    # ersa_upper_bound is_niece_aunt final_degree total_seg_len seg_count
    relatives = pd.DataFrame()
    relatives.loc[:, 'id1'] = tribes.Id1
    relatives.loc[:, 'id2'] = tribes.Id2
    relatives.loc[:, 'king_degree'] = map_tribes_degree(tribes.EstDegree)
    relatives.loc[:, 'king_relation'] = tribes.EstDegree
    relatives.loc[:, 'ersa_degree'] = map_tribes_degree(tribes.EstDegree)
    relatives.loc[:, 'ersa_lower_bound'], relatives.loc[:, 'ersa_upper_bound'] = tribes_degree_to_interval(relatives.ersa_degree)
    relatives.loc[:, 'final_degree'] = map_tribes_degree(tribes.EstDegree)
    relatives.to_csv(output_path, sep='\t', index=False)


def interval_evaluate(result, fam, plot_name, pr_plot_name, conf_matrix_plot_name, output_path, only_client=False, source='both'):
    iids, pedigree = read_pedigree(fn=fam)
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
    if not total:
        print('total is not total')
        return
    print('correct: ', correct)
    print('total: ', total)
    keys = sorted(list(confusion_matrix.keys()))
    for key in keys:
        print(f'{key}\t{confusion_matrix[key]}')
    print()

    compare(total, correct, plot_name)
    interval_precision_recall(kinship, inferred, clients, source, pr_plot_name)

    plot_confusion_matrix(confusion_matrix, conf_matrix_plot_name)

    kinship_frame = kinship_to_dataframe(kinship)
    predictions = pd.read_table(result, index_col=['id1', 'id2'])

    merged = predictions.merge(kinship_frame, left_index=True, right_index=True, how='outer')
    merged.to_csv(output_path, sep='\t')


if __name__ == '__main__':
    tribes_rel = 'analysis/data/tribes/input_GRM-allchr_FPI_IBD.csv'
    rel = 'analysis/data/tribes/relatives.tsv'
    fam = '/media_ssd/pipeline_data/big_ibis_affymetrix_ersa1/pedsim/simulated/reheaded_data.fam'
    accuracy_path = 'analysis/data/tribes/accuracy.png'
    pr_path = 'analysis/data/tribes/pr.png'
    conf_matrix_path = 'analysis/data/tribes/conf_matrix.png'
    updated_rel_path = 'analysis/data/tribes/updated_relatives.tsv'
    tribes_to_relatives(tribes_rel, rel)
    interval_evaluate(rel,
                      fam,
                      accuracy_path,
                      pr_path,
                      conf_matrix_path,
                      updated_rel_path,
                      only_client=True,
                      source='both')