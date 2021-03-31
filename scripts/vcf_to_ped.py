import allel
import zarr
import numpy
import os
import time


def write_ped_buffer(buffer: dict, out: str):
    start = time.time()
    with open(out, 'a') as ped:
        lines = []
        for sample, gt in buffer.items():
            #gt_str = numpy.array2string(gt, separator=' ', max_line_width=gt.shape[0], threshold=gt.shape[0])
            line = ' '.join(['0', sample, '0', '0', '0', '0']) + ' ' + ' '.join(gt) + '\n'
            print(sample, len(gt))
            lines.append(line)

        ped.writelines(lines)

    end = time.time()
    print(f'write_ped_buffer took {end - start} seconds')


def write_map_file(ids, chrom, positions, out: str):
    with open(out, 'w') as file:
        for _id, _chr, pos in zip(ids, chrom, positions):
            file.write('\t'.join([_chr, _id, '0', str(pos)]) + '\n')


if __name__ == '__main__':
    start = time.time()

    vcf = snakemake.input['vcf'][0]
    zarr_cache = snakemake.params['zarr']
    ped = snakemake.output['ped']
    map_file = snakemake.output['map_file']

    '''
    vcf = 'test_data/imputed_chr9.vcf.gz'
    zarr_cache = 'test_data/imputed_chr9.zarr'
    ped = 'test_data/imputed_chr9.ped'
    map_file = 'test_data/imputed_chr9.map'
    '''

    '''
    vcf = '/home/ag3r/vcf_to_ped/testsample.vcf.gz'
    zarr_cache = '/home/ag3r/vcf_to_ped/testsample.zarr'
    ped = '/home/ag3r/vcf_to_ped/testsample.ped'
    map_file = '/home/ag3r/vcf_to_ped/testsample.map'
    '''

    if os.path.exists(ped):
        os.remove(ped)

    allel.vcf_to_zarr(
        input=vcf,
        output=zarr_cache,
        overwrite=True,
        fields=['variants/CHROM', 'variants/ID', 'variants/POS', 'variants/REF', 'variants/ALT', 'calldata/GT', 'samples']
    )

    array = zarr.open_group(zarr_cache, mode='r')

    samples_in_chunk = 100
    ids = array['variants/ID'][:]
    chrom = array['variants/CHROM'][:]
    positions = array['variants/POS'][:]
    write_map_file(ids, chrom, positions, map_file)

    # only for the single chromosome!
    print(f'written total of {len(positions)} variants from chr{chrom[0]} to {map_file}')

    for sample_start in range(0, len(array['samples']), samples_in_chunk):
        sample_end = min(len(array['samples']), sample_start + samples_in_chunk)
        samples = array['samples'][sample_start: sample_end]
        ped_buffer = {}

        # [variants]
        ref = array['variants']['REF'][:].reshape(-1, 1)
        # [variants, alt_count]
        alt = array['variants']['ALT'][:]

        # [variants, alt_count + 1]
        alleles = numpy.hstack([ref, alt])
        alleles[alleles == ''] = ' '
        # [variants, samples, 2] # 2 because we use biallelic only
        gt = array['calldata']['GT'][:, sample_start: sample_end, :]
        for sample_idx, sample in enumerate(samples):

            # [variants, 2]
            line = numpy.zeros((gt.shape[0], 2), dtype=alleles.dtype)

            # [variants] contains 0 or 1, i.e. count of reference alleles in corresponding haplotype (fwd or bck)
            fwd_idx = gt[:, sample_idx, 0].astype(bool)
            bck_idx = gt[:, sample_idx, 1].astype(bool)


            #print(gt.shape, alleles.shape, gt[:, sample_idx, 0].shape)
            # [variants] ref_line contains ref variants, alt_line contains alf_variants
            ref_line, alt_line = alleles[~fwd_idx, 0], alleles[fwd_idx, 1]
            # we are setting first haplotype ref and alt variants in corresponding positions
            line[~fwd_idx, 0] = ref_line
            line[fwd_idx, 0] = alt_line

            ref_line, alt_line = alleles[~bck_idx, 0], alleles[bck_idx, 1]
            line[~bck_idx, 1] = ref_line
            line[bck_idx, 1] = alt_line

            line_array = line.ravel()
            #print(line_array[:5], line_array.shape)

            #for variant_idx in range(gt.shape[0]):
            #    fwd_idx, bck_idx = gt[variant_idx, sample_idx, 0], gt[variant_idx, sample_idx, 1]
            #    fwd = alleles[variant_idx, fwd_idx]
            #    ped_line.append(fwd)

            #    bck = alleles[variant_idx, bck_idx]
            #    ped_line.append(bck)

            ped_buffer[sample] = line_array

        write_ped_buffer(ped_buffer, ped)
        print(f'written {len(ped_buffer)} samples to {ped}')
    end = time.time()
    print(f'total execution time is {end - start:.4f} seconds')
    print(f'written total of {len(array["samples"])} and {len(array["variants"]["REF"])} variants to {ped}')
