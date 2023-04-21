import pandas
import logging
from utils.bcftools import bcftools_stats
from utils.bcftools import bcftools_query
from utils.bcftools import bcftools_view
from io import StringIO


def find_outliers(psc: pandas.DataFrame, outliers_file_path: str, keep_samples_file_path: str, alpha: float):

    q1 = psc.nNonMissing.quantile(0.25)
    q3 = psc.nNonMissing.quantile(0.75)
    iqr = q3-q1
    lower_bound_mask = psc.nNonMissing < (q1 - alpha*iqr)
    upper_bound_mask = psc.nNonMissing > (q3 + alpha*iqr)
    outliers = psc[lower_bound_mask | upper_bound_mask]
    if outliers.empty:
        outliers['exclusion_reason'] = []
    else:
        outliers.loc[lower_bound_mask, 'exclusion_reason'] = f'Sample was below 1st quantile - IQR*{alpha} ({q1 - alpha*iqr}) by Missing Samples'
        outliers.loc[upper_bound_mask, 'exclusion_reason'] = f'Sample was above 3rd quantile + IQR*{alpha} ({q3 + alpha*iqr}) by Missing Samples'
    keep_samples = pandas.concat([psc, outliers, outliers]).drop_duplicates(keep=False)

    outliers_list = list(outliers.sample_id)
    with open(outliers_file_path, 'w') as outliers_file:
        for outlier in outliers_list:
            outliers_file.write(outlier + '\n')

    keep_samples_list = list(keep_samples.sample_id)
    with open(keep_samples_file_path, 'w') as keep_samples_file:
        for sample in keep_samples_list:
            keep_samples_file.write(sample + '\n')

    return outliers


def get_stats(vcf_input, samples_path, psc_path=False, stats_path=''):

    bcftools_query(vcf_input, list_samples=True, save_output=samples_path)
    raw_table = bcftools_stats(vcf_input, samples_path, save_output=psc_path, save_stats=stats_path)
    names = [
        'PSC',
        'id',
        'sample_id',
        'nRefHom',
        'nNonRefHom',
        'nHets',
        'nTransitions',
        'nTransversions',
        'nIndels',
        'average depth',
        'nSingletons',
        'nHapRef',
        'nHapAlt',
        'nMissing'
    ]
    psc = pandas.read_table(StringIO(raw_table), header=None, names=names)

    psc.loc[:, 'nNonMissing'] = psc.nRefHom + psc.nNonRefHom + psc.nHets
    psc.loc[:, 'missing_share'] = psc.nMissing / (psc.nMissing + psc.nNonMissing)
    psc.loc[:, 'alt_hom_share'] = psc.nNonRefHom / psc.nNonMissing
    psc.loc[:, 'het_samples_share'] = psc.nHets / psc.nNonMissing

    return psc



if __name__ == '__main__':

    input_vcf = snakemake.input['vcf']
    stats_file = snakemake.output['stats']
    samples = snakemake.params['samples']
    psc_path = snakemake.params['psc']
    keep_samples = snakemake.params['keep_samples']
    bad_samples_path = snakemake.output['bad_samples']
    report_path = snakemake.output['report']
    outliers = snakemake.output['outliers']
    no_outliers_vcf = snakemake.output['no_outliers_vcf']
    missing_samples = float(snakemake.params['missing_samples'])
    alt_hom_samples = float(snakemake.params['alt_hom_samples'])
    het_samples = float(snakemake.params['het_samples'])
    alpha = float(snakemake.params['iqr_alpha'])

    logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG, format='%(levelname)s:%(asctime)s %(message)s')

    # get initial stats, that we do not save
    logging.info("Getting initial vcf stats...")
    psc = get_stats(vcf_input=input_vcf,
                    samples_path=samples,
                    psc_path=psc_path,
                    stats_path=stats_file)
    # filter outliers
    outliers = find_outliers(psc, outliers_file_path=outliers, keep_samples_file_path=keep_samples, alpha=alpha)
    logging.info(f"Found {len(outliers)} outliers!")
    if not outliers.empty:
        bcftools_view(input_vcf, no_outliers_vcf, keep_samples)
    else:
        with open(no_outliers_vcf, 'w') as dummy_no_outliers_vcf:
            dummy_no_outliers_vcf.write('')
        no_outliers_vcf = input_vcf

    # get final stats without outliers
    logging.info("Getting final vcf stats...")
    psc = get_stats(vcf_input=no_outliers_vcf,
                    samples_path=samples,
                    psc_path=psc_path,
                    stats_path=stats_file)

    bad_missing_samples_mask = (psc.missing_share >= missing_samples / 100) | (psc.nMissing + psc.nNonMissing == 0)

    bad_alt_hom_samples_mask = (psc.alt_hom_share < alt_hom_samples / 100) | (psc.nNonMissing == 0)

    bad_het_samples_mask = (psc.het_samples_share < het_samples / 100) | (psc.nNonMissing == 0)

    psc.loc[bad_het_samples_mask, 'exclusion_reason'] = f'Sample has <= {het_samples}% heterozygous SNPs'
    psc.loc[bad_alt_hom_samples_mask, 'exclusion_reason'] = f'Sample has <= {alt_hom_samples}% homozygous alternative SNPs'
    psc.loc[bad_missing_samples_mask, 'exclusion_reason'] = f'Sample has >= {missing_samples}% missing SNPs'
    bad_samples = psc.loc[(bad_alt_hom_samples_mask | bad_het_samples_mask | bad_missing_samples_mask),
                          ['sample_id', 'missing_share', 'alt_hom_share', 'het_samples_share', 'exclusion_reason']]
    psc = psc.append(outliers[['sample_id', 'missing_share', 'alt_hom_share', 'het_samples_share', 'exclusion_reason']], ignore_index=True)

    samples_only = bad_samples.loc[:, ['sample_id']].copy()
    # if input vcf file has iids in form fid_iid we split it, else we just assign fid equal to iid
    samples_only.loc[:, 'fid'] = [s.split('_')[0] if '_' in s else s for s in samples_only.sample_id]
    samples_only.loc[:, 'iid'] = [s.split('_')[1] if '_' in s else s for s in samples_only.sample_id]

    samples_only.loc[:, ['fid', 'iid']].to_csv(bad_samples_path, sep='\t', header=None, index=False)
    bad_samples.to_csv(report_path, sep='\t')

    log = f'We have total of {psc.shape[0]} samples\n' \
          f'{bad_missing_samples_mask.sum()} samples have >= {missing_samples}% missing share\n' \
          f'{bad_alt_hom_samples_mask.sum()} samples have <= {alt_hom_samples}% homozygous alternative variants\n' \
          f'{bad_het_samples_mask.sum()} samples have <= {het_samples}% heterozygous variants\n' \
          f'{len(outliers)} samples detected as outliers\n'

    logging.info(log)
    print(log)
