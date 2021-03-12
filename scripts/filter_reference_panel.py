import numpy
import matplotlib.pyplot as plt


if __name__ == '__main__':
    file = '/media/ref/1000genome/allele_info/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.only_rs.biallelic.tab'
    acs = []
    frequent = []
    with open(file, 'r') as tab:
        bad_lines = 0
        good_lines = 0
        for i, line in enumerate(tab):
            try:
                items = line.split(' ')
                rsid = items[2]
                AC = int(items[-1].split(';')[0].split('=')[1])
                acs.append(AC/5008)
                good_lines += 1
                frequent.append(rsid)
            except ValueError as e:
                bad_lines += 1
                continue

    print(f'we have total of {good_lines} good SNPs and {bad_lines} bad snps')
    hist, edges = numpy.histogram(acs, bins=100, density=False)
    print(numpy.cumsum(hist))

    plt.figure(figsize=(17, 12))
    plt.grid()
    plt.hist(acs, bins=100, density=True)
    plt.savefig('AF_1000GENOMES.png')
