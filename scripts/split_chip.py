import os
import json


def split_chip():
    num_runs = int(snakemake.params['num_runs'])
    num_founders = int(snakemake.params['num_founders'])
    with open(snakemake.input['samples'], "r") as f:
        lines = f.readlines()
        samples_list = [line.rstrip() for line in lines]
        
    for i in range(0, num_runs*num_founders, num_founders):
        sample = samples_list[i:i + num_founders]
        with open("sample.txt", "w") as s:
            s.write("\n".join(sample))
        seg = int(i/num_founders)
        print(f"Working on segment{seg}.vcf.gz...")
        os.system(f"bcftools view -S sample.txt pedsim/phased/background.bcf.gz --force-samples -O z -o background/segment{seg}.vcf.gz")


if __name__ == '__main__':
    split_chip()
