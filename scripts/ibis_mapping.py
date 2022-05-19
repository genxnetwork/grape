import pandas
import numpy
import os
import re

def bim_mapping(bim_file: str, cm_path: str, output_path: str):
    """bim file mapping.
        Args:
            bim_file: .bim file path
            cm_path: path to flder where all map file chromosomes are located
            output_path: output .bim file
    """
    def atoi(text):
        return int(text) if text.isdigit() else text

    def natural_keys(text):
        return [atoi(c) for c in re.split(r'(\d+)', text)]
    chromosomes_map = []
    cmss = []
    bim = pandas.read_table(bim_file, header=None)
    chromosomes = os.listdir(cm_path)
    chromosomes.sort(key = natural_keys)
    for chr in chromosomes:
        map = pandas.read_table(os.path.join(cm_path, chr))
        chromosomes_map.append(map)
    for chr in list(bim.iloc[:, 0].unique()):
        chrom = bim[bim.iloc[:, 0] == chr]
        pos = chrom.iloc[:, 3].values
        cms = numpy.interp(pos, chromosomes_map[chr-1].iloc[:, 1].values, chromosomes_map[chr-1].iloc[:, -1].values)

        eps = 1e-7
        new_cms = []
        for i in range(0, cms.shape[0]):
            if i == 0:
                new_cms.append(cms[i])
                continue

            if new_cms[-1] >= cms[i]:
                new_cms.append(new_cms[-1] + eps)
            else:
                new_cms.append(cms[i])

        print(f'cms of len {len(cms)} after:')

        cmss += new_cms
    bim.iloc[:, 2] = cmss
    bim.to_csv(output_path, sep="\t", header=False, index=False)

