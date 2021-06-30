import pandas
import logging


if __name__ == '__main__':

    psc_path = snakemake.input['psc']
    bad_samples_path = snakemake.output['bad_samples']
    report_path = snakemake.output['report']

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

    bad_missing_samples_mask = (psc.missing_share > 0.1) | (psc.nMissing + psc.nNonMissing == 0)

    bad_alt_hom_samples_mask = ((psc.nNonRefHom / psc.nNonMissing) <= 0.01) | (psc.nNonMissing == 0)

    bad_het_samples_mask = ((psc.nHets / psc.nNonMissing) <= 0.05) | (psc.nNonMissing == 0)

    psc.loc[:, 'alt_hom_share'] = psc.nNonRefHom / psc.nNonMissing
    psc.loc[:, 'het_samples_share'] = psc.nHets / psc.nNonMissing

    psc.loc[bad_het_samples_mask, 'exclusion_reason'] = 'Sample has < 5% heterozygous SNPs'
    psc.loc[bad_alt_hom_samples_mask, 'exclusion_reason'] = 'Sample has < 1% homozygous alternative SNPs'
    psc.loc[bad_missing_samples_mask, 'exclusion_reason'] = 'Sample has > 10% missing SNPs'
    bad_samples = psc.loc[(bad_alt_hom_samples_mask | bad_het_samples_mask | bad_missing_samples_mask),
                          ['sample_id', 'missing_share', 'alt_hom_share', 'het_samples_share', 'exclusion_reason']]

    samples_only = bad_samples.loc[:, ['sample_id']].copy()
    samples_only.loc[:, 'fid'] = [s.split('_')[0] for s in samples_only.sample_id]
    samples_only.loc[:, 'iid'] = [s.split('_')[1] for s in samples_only.sample_id]

    samples_only.loc[:, ['fid', 'iid']].to_csv(bad_samples_path, sep='\t', header=None, index=False)
    bad_samples.to_csv(report_path, sep='\t')

    logging.info(f'We have total of {psc.shape[0]} samples')
    logging.info(f'{bad_missing_samples_mask.sum()} samples have >10% missing share')
    logging.info(f'{bad_alt_hom_samples_mask.sum()} samples have < 1% homozygous alternative variants')
    logging.info(f'{bad_het_samples_mask.sum()} samples have < 5% heterozygous variants')

    print(f'We have total of {psc.shape[0]} samples')
    print(f'{bad_missing_samples_mask.sum()} samples have >10% missing share')
    print(f'{bad_alt_hom_samples_mask.sum()} samples have < 1% homozygous alternative variants')
    print(f'{bad_het_samples_mask.sum()} samples have < 5% heterozygous variants')
