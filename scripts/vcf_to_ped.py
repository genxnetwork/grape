import allel
import zarr
import numpy
import os


def write_ped_buffer(buffer: dict, out: str):
    with open(out, 'a') as ped:
        for sample, gt in buffer.items():
            line = ' '.join(['0', sample, '0', '0', '0', '0'] + gt)
            ped.write(line + '\n')


if __name__ == '__main__':

    vcf = snakemake.input['vcf'][0]
    zarr_cache = snakemake.params['zarr']
    ped = snakemake.output['ped']

    #vcf = 'test_data/test.vcf.gz'
    #zarr_cache = 'test_data/test.zarr'
    #ped = 'test_data/test.ped'

    if os.path.exists(ped):
        os.remove(ped)

    allel.vcf_to_zarr(
        input=vcf,
        output=zarr_cache,
        overwrite=True,
        fields=['variants/CHROM', 'variants/POS', 'variants/REF', 'calldata/GT', 'variants/ALT', 'samples']
    )

    array = zarr.open_group(zarr_cache, mode='r')

    samples_in_chunk = 100

    for sample_start in range(0, len(array['samples']), samples_in_chunk):
        sample_end = sample_start + samples_in_chunk
        samples = array['samples'][sample_start: sample_end]
        ped_buffer = {}

        # [variants, samples]
        ref = array['variants']['REF'][:].reshape(-1, 1)
        # [variants, samples]
        alt = array['variants']['ALT'][:]

        alleles = numpy.hstack([ref, alt])

        # [variants, samples, 2] # 2 because we use biallelic only
        gt = array['calldata']['GT'][:, sample_start: sample_end, :]
        for sample_idx, sample in enumerate(samples):
            ped_line = []
            for variant_idx in range(gt.shape[0]):
                fwd_idx, bck_idx = gt[variant_idx, sample_idx, 0], gt[variant_idx, sample_idx, 1]
                ped_line.append(alleles[variant_idx, fwd_idx])
                ped_line.append(alleles[variant_idx, bck_idx])

            ped_buffer[sample] = ped_line

        write_ped_buffer(ped_buffer, ped)
        print(f'written {len(ped_buffer)} samples to {ped}')

    print(f'written total of {len(array["samples"])} and {len(array["variants"]["REF"])} variants to {ped}')










