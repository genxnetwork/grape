import subprocess
import logging

if __name__ == '__main__':

    vcf = snakemake.input['vcf']
    bad_samples = snakemake.input['bad_samples']
    out = snakemake.params['out']
    batch = snakemake.params['batch']

    logging.basicConfig(filename=snakemake.log[0],
                        level=logging.DEBUG,
                        format='%(levelname)s:%(asctime)s %(message)s')

    p_freqx = subprocess.run(f'plink --vcf {vcf} --freqx --out plink/{out}',
                             shell=True,
                             capture_output=True)

    p_remove = subprocess.run(f'plink --vcf {vcf} --remove {bad_samples} '
                              f'--make-bed --keep-allele-order --out plink/{out}',
                              shell=True,
                              capture_output=True)

    stdout_freqx = p_freqx.stdout.decode()
    stderr_freqx = p_freqx.stderr.decode()
    stdout_remove = p_remove.stderr.decode()
    stderr_remove = p_remove.stdout.decode()

    logging.info(f'{stdout_remove}\n{stdout_freqx}'
                 f'\n{stderr_remove}\n{stderr_freqx}'.decode('utf-8'))

    if p_remove.returncode != 0 or p_freqx.returncode != 0:
        with open('pass_batches.list', 'r') as list:
            lines = list.readlines()
        with open('pass_batches.list', 'w') as list:
            for line in lines:
                if line.strip('\n') != f'{batch}':
                    list.write(line)
                else:
                    print(f'removed {line}!')
