import pandas as pd
import networkx as nx
import logging


if __name__ == '__main__':

    logging.basicConfig(filename='log.log', level=logging.DEBUG, format='%(levelname)s:%(asctime)s %(message)s')

    relatives_path = snakemake.input['relatives']
    # within families
    # across families
    output_path = snakemake.output['remove_list']
    
    relatives = pd.read_csv(relatives_path, sep = '\t')
    total_rel = len(list(set(list(relatives['id1'].unique()) + list(relatives['id2'].unique()))))
    relatives = relatives.loc[relatives.final_degree < 7]
    G = nx.Graph()
    G.add_weighted_edges_from(zip(relatives.id1, relatives.id2, relatives.final_degree))

    def rm_node(g):
      remove_list = []
      degree_dict = {node:val for (node, val) in g.degree()}
      while max(degree_dict.values()) > 0:
        degree_dict = {node:val for (node, val) in g.degree()}
        max_node = max(degree_dict, key=degree_dict.get)
        remove_list.append(max_node)
        g.remove_node(max_node)
      return remove_list

    remove_list = rm_node(G)
    print(f"{total_rel - len(remove_list)} out of {total_rel} relatives left after cleaning, {len(remove_list)} was removed!")
        
    with open(output_path, 'w') as rm:
      rm.write('\n'.join(remove_list))
