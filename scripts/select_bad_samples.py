import pandas
import logging


if __name__ == '__main__':

    psc_path = snakemake.input['psc']
    bad_samples_path = snakemake.output['bad_samples']
    report_path = snakemake.output['report']
    missing_samples = float(snakemake.params['missing_samples'])
    alt_hom_samples = float(snakemake.params['alt_hom_samples'])
    het_samples = float(snakemake.params['het_samples'])

    logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG, format='%(levelname)s:%(asctime)s %(message)s')

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
    psc = pandas.read_table(psc_path, header=None, names=names)

    psc.loc[:, 'nNonMissing'] = psc.nRefHom + psc.nNonRefHom + psc.nHets
    psc.loc[:, 'missing_share'] = psc.nMissing / (psc.nMissing + psc.nNonMissing)
    psc.loc[:, 'alt_hom_share'] = psc.nNonRefHom / psc.nNonMissing
    psc.loc[:, 'het_samples_share'] = psc.nHets / psc.nNonMissing

    bad_missing_samples_mask = (psc.missing_share >= missing_samples / 100) | (psc.nMissing + psc.nNonMissing == 0)

    bad_alt_hom_samples_mask = (psc.alt_hom_share <= alt_hom_samples/100) | (psc.nNonMissing == 0)

    bad_het_samples_mask = (psc.het_samples_share <= het_samples/100) | (psc.nNonMissing == 0)

    psc.loc[bad_het_samples_mask, 'exclusion_reason'] = f'Sample has <= {het_samples}% heterozygous SNPs'
    psc.loc[bad_alt_hom_samples_mask, 'exclusion_reason'] = f'Sample has <= {alt_hom_samples}% homozygous alternative SNPs'
    psc.loc[bad_missing_samples_mask, 'exclusion_reason'] = f'Sample has >= {missing_samples}% missing SNPs'
    bad_samples = psc.loc[(bad_alt_hom_samples_mask | bad_het_samples_mask | bad_missing_samples_mask),
                          ['sample_id', 'missing_share', 'alt_hom_share', 'het_samples_share', 'exclusion_reason']]

    samples_only = bad_samples.loc[:, ['sample_id']].copy()
    # if input vcf file has iids in form fid_iid we split it, else we just assign fid equal to iid
    samples_only.loc[:, 'fid'] = [s.split('_')[0] if '_' in s else s for s in samples_only.sample_id]
    samples_only.loc[:, 'iid'] = [s.split('_')[1] if '_' in s else s for s in samples_only.sample_id]

    samples_only.loc[:, ['fid', 'iid']].to_csv(bad_samples_path, sep='\t', header=None, index=False)
    bad_samples.to_csv(report_path, sep='\t')

    log = f'We have total of {psc.shape[0]} samples\n' \
          f'{bad_missing_samples_mask.sum()} samples have >= {missing_samples}% missing share\n' \
          f'{bad_alt_hom_samples_mask.sum()} samples have <= {alt_hom_samples}% homozygous alternative variants\n' \
          f'{bad_het_samples_mask.sum()} samples have <= {het_samples}% heterozygous variants'

    logging.info(log)
    print(log)
