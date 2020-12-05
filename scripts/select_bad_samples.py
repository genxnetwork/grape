import pandas


if __name__ == '__main__':

    psc_path = snakemake.input['psc']
    bad_samples_path = snakemake.output['bad_samples']

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

    bad_alt_hom_samples_mask = (psc.nNonRefHom / psc.nNonMissing) <= 0.01
    bad_het_samples_mask = (psc.nHets / psc.nNonMissing) <= 0.05

    bad_samples = psc.loc[(bad_alt_hom_samples_mask | bad_het_samples_mask), ['sample_id']]
    bad_samples.to_csv(bad_samples_path, sep='\t', header=None)

    print(f'We have total of {psc.shape[0]} samples')
    print(f'{bad_alt_hom_samples_mask.sum()} samples have < 1% homozygous alternative variants')
    print(f'{bad_het_samples_mask.sum()} samples have < 5% heterozygous variants')
