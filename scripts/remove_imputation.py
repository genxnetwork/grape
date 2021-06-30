import gzip
import logging


if __name__ == '__main__':
    genotyped_count = 0
    input_file = snakemake.input['vcf']
    output_file = snakemake.output['vcf']

    logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG, format='%(levelname)s:%(asctime)s %(message)s')

    with gzip.open(output_file, 'ab') as gzab:
        with gzip.open(input_file, 'r') as gzr:
            for cnt, line in enumerate(gzr):
                if (cnt + 1) % 100000 == 0:
                    logging.info(f'processed {cnt} variants')
                if b'IMPUTED' not in line:
                    genotyped_count += 1
                    gzab.write(line)
    logging.info(f'processed {cnt} variants, removed {cnt - genotyped_count} variants, {genotyped_count} variants remained')
