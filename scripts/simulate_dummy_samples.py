import logging
from collections import namedtuple
import allel
import numpy
import gzip


def write_vcf_file(file_name, headers, calldata, variant_names, sample_names, sample_genotypes):
    with gzip.open(file_name, 'wt') as vcf_file:
        # Writing metadata headers (e.g., file format, reference, etc.)
        for header in headers:
            vcf_file.write(header + "\n")
            
        
        # Writing column header
        columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
        columns.extend(sample_names)
        vcf_file.write("\t".join(columns) + "\n")
        
        # Writing variant data
        variant_count = len(variant_names)
        # for tqdm_idx in tqdm.tqdm(range(variant_count)):
        for tqdm_idx in range(variant_count):
            var_name, calldata_item, sample_genotype = variant_names[tqdm_idx], calldata[tqdm_idx], sample_genotypes[tqdm_idx, :]
            chrom, pos, id, ref, alt = var_name
            qual, filter, format = calldata_item
            info = 'AF=0.5'
            data_line = [chrom, pos, id, ref, alt, qual, filter, info, format]
            data_line.extend([f'{sg[0]}|{sg[1]}' for sg in sample_genotype])  # add the sample genotype data here
            vcf_file.write("\t".join(map(str, data_line)) + "\n")

'''
# Example usage:
# You can create your input data as per your requirement
headers = ["fileformat=VCFv4.2", "reference=file:///path_to_reference_genome"]
calldata = [(None, None, "GT"), (None, None, "GT")]  # Example calldata, replace with actual data
variant_names = [("chr1", "12345", ".", "A", "C"), ("chr1", "12346", ".", "G", "T")]
variant_info = ["NS=3;DP=14;AF=0.5;DB;H2", "NS=3;DP=11;AF=0.017"]
sample_names = ["sample1", "sample2"]
sample_genotypes = [["0/1", "0/0"], ["1/1", "0/1"]]

# Writing to file
write_vcf_file("output.vcf", headers, calldata, variant_names, variant_info, sample_names, sample_genotypes)
'''

if __name__ == '__main__':

    try:
        snakemake
    except NameError:
        Snakemake = namedtuple('Snakemake', ['input', 'output', 'params', 'log'])
        snakemake = Snakemake(
            input={'background': 'test_data/vcf/background_20.vcf.gz'},
            output=['test_data/vcf/augmentated_background_20.vcf.gz'],
            params={'sample_count': 100},
            log=['test_data/merge_king_ersa.log']
        )

    logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG, format='%(levelname)s:%(asctime)s %(message)s')

    # We read background.vcf.gz file
    
    vcf = allel.read_vcf(snakemake.input['background'])
    print(f'vcf has {vcf.keys()} fields')
    print(f'We read vcf with {len(vcf["samples"])} samples from {snakemake.input["background"]}')
    # the first dimension corresponds to the variants genotyped
    # the second dimension corresponds to the samples genotyped
    # the third dimension corresponds to the ploidy of the samples.
    gt = vcf['calldata/GT']
    samples = vcf['samples']
    # We sample pair of samples from the file
    indices = numpy.arange(gt.shape[0]) # number of variants
    new_samples = []
    new_gt = []
    # ugly hack because bcftools somehow cannot divide multiallelic sites to biallelic ones
    # and we need data to be biallelic only
    gt[gt >= 2] = 0
    for i in range(0, len(samples)):
        if len(new_samples) >= snakemake.params['sample_count']:
            break
        for j in range(i + 1, len(samples)):
            sample_gt_i = gt[:, i, :]
            sample_gg_j = gt[:, j, :]
            # We generate a mix of these two samples
            change_indices = numpy.random.choice(indices, size=gt.shape[1] // 2, replace=False)
            sample_gt_i[change_indices, :] = sample_gg_j[change_indices, :]
            new_samples.append(samples[i].replace('_', '-') + '_' + samples[j].replace('_', '-'))
            new_gt.append(numpy.expand_dims(sample_gt_i, axis=1))
            if len(new_samples) >= snakemake.params['sample_count']:
                break        
    
    print(f'generated {len(new_samples)} samples')
    new_samples = numpy.array(new_samples)
    new_gt = numpy.concatenate(new_gt, axis=1)
    
    # headers = ["fileformat=VCFv4.2", "reference=file:///path_to_reference_genome"]
    with gzip.open(snakemake.input['background'], 'rt') as vcf_file:
        headers = []
        for i, line in enumerate(vcf_file):
            if line.startswith('##'):
                headers.append(line.strip())
            if i > 2 and not line.startswith('##'):
                break
    headers.append('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">')
    # qual, filter, format
    calldata = [(qual, 'PASS', 'GT') for qual, fp in zip(vcf['variants/QUAL'], vcf['variants/FILTER_PASS'])]
    print(f'calldata {calldata[:10]}')
    variants = [(chrom, pos, _id, ref, alt[0]) for chrom, pos, _id, ref, alt in zip(vcf['variants/CHROM'], 
                                                                                 vcf['variants/POS'], 
                                                                                 vcf['variants/ID'], 
                                                                                 vcf['variants/REF'], 
                                                                                 vcf['variants/ALT'])]
    
    samples = numpy.concatenate([vcf['samples'], new_samples], axis=0)
    gt_to_write = numpy.concatenate([gt, new_gt], axis=1)
    write_vcf_file(snakemake.output[0], headers, calldata, variants, samples, gt_to_write)
    
    test_read_vcf = allel.read_vcf(snakemake.output[0])
        