import pandas as pd

if __name__ == '__main__':

  bim_mapped_path = snakemake.input['bim_mapped']
  bim_path = snakemake.output['bim']

  bim_mapped = pd.read_csv(bim_mapped_path, sep='\t', header=None)
  bim_mapped.iloc[:, 2] = 0
  bim_mapped.to_csv(bim_path, sep="\t", header=False, index=False)
