import numpy
import pandas
from utils.bcftools import bcftools_view


'''
Index(['Genetic ID', 'Master ID', 'Skeletal code', 'Skeletal element',
       'Year data from this individual was first published [for a present-day individuals we give the data of the data reported here; missing GreenScience 2010 (Vi33.15, Vi33.26), Olalde2018 (I2657), RasmussenNature2010 (Australian)]',
       'Publication',
       'Method for Determining Date; unless otherwise specified, calibrations use 95.4% intervals from OxCal v4.4.2 Bronk Ramsey (2009); r5; Atmospheric data from Reimer et al (2020)',
       'Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]',
       'Date standard deviation in BP [OxCal sigma for a direct radiocarbon date, and standard deviation of the uniform distribution between the two bounds for a contextual date]',
       'Full Date One of two formats. (Format 1) 95.4% CI calibrated radiocarbon age (Conventional Radiocarbon Age BP, Lab number) e.g. 2624-2350 calBCE (3990±40 BP, Ua-35016). (Format 2) Archaeological context range, e.g. 2500-1700 BCE',
       'Age at Death from physical anthropology', 'Group ID', 'Locality',
       'Political Entity', 'Lat.', 'Long.', 'Pulldown Strategy', 'Data source',
       'No. Libraries',
       '1240k coverage (taken from original pulldown where possible)',
       'SNPs hit on autosomal targets (Computed using easystats on 1240k snpset)',
       'SNPs hit on autosomal targets (Computed using easystats on HO snpset)',
       'Molecular Sex', 'Family ID and position within family',
       'Y haplogroup (manual curation in terminal mutation format)',
       'Y haplogroup (manual curation in ISOGG format)',
       'mtDNA coverage (merged data)', 'mtDNA haplogroup if >2x or published',
       'mtDNA match to consensus if >2x (merged data)',
       'Damage rate in first nucleotide on sequences overlapping 1240k targets (merged data)',
       'Sex ratio [Y/(Y+X) counts] (merged data)',
       'Library type (minus=no.damage.correction, half=damage.retained.at.last.position, plus=damage.fully.corrected, ds=double.stranded.library.preparation, ss=single.stranded.library.preparation)',
       'Libraries', 'ASSESSMENT',
       'ASSESSMENT WARNINGS (Xcontam interval is listed if lower bound is >0.005, "QUESTIONABLE" if lower bound is 0.01-0.02, "QUESTIONABLE_CRITICAL" or "FAIL" if lower bound is >0.02) (mtcontam confidence interval is listed if coverage >2 and upper bound is <0.'],
      dtype='object')
'''


if __name__ == '__main__':
    cols = [
        'genetic_id', 'master_id', 'skeletal_code', 'skeletal_element',
        'year_data_from_this_individual_was_first_published',
        'publication',
        'method_for_determining_date',
        'date_mean',
        'date_sd',
        'full_date',
        'age_at_death_from_physical_anthropology', 'group_id', 'locality', 
        'political_entity', 'lat', 'long', 'pulldown_strategy', 'data_source',
        'no_libraries',
        'coverage', 'snps_hit_autosomal_1240k', 'snps_hit_autosomal_ho',
        'sex', 'family_id', 'y_haplogroup_terminal', 'y_haplogroup_isogg', 'mt_coverage', 'mt_haplogroup', 'mt_dna_match', 'damage_rate', 'sex_ratio', 'library_type', 'libraries', 'assessment', 'assessment_warnings'
    ]
    anno = pandas.read_table('workdir/v54.1.p1_1240K_public.tsv', header=0, names=cols, index_col='genetic_id')
    present = anno[anno['full_date'] == 'present']
    # group by publication and print counts
    print(present.groupby('publication').size())
    
    # plot hist of snps hit on autosomal targets
    present.loc[:, 'snps_hit_autosomal_1240k'] = present['snps_hit_autosomal_1240k'].astype(float)
    (present['snps_hit_autosomal_1240k']/present['snps_hit_autosomal_1240k'].max()).hist(bins=20).get_figure().savefig('workdir/snps_hit_autosomal_1240k.png')
    present.loc[:, 'coverage'] = present['coverage'].apply(lambda x: float(x) if x != '..' else numpy.nan).astype(float)
    present.reset_index()['genetic_id'].to_csv('workdir/present.txt', index=False, header=None)
    print(f'{len(present)} individuals are modern')
    
    bcftools_view('workdir/aadr.vcf.gz',  'workdir/aadr.present.vcf.gz', 'workdir/present.txt')