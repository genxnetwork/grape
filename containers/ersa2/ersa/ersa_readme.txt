Estimation of Recent Shared Ancestry (ERSA)
ERSA version 2.0. January 31st, 2014. 
Copyright (c) 2011 Chad Huff and David Witherspoon. All rights reserved. Patent pending. 

If results from ERSA are used in a published work, please cite: 
Huff CD*, Witherspoon D*, Simonson T, Xing J, Watkins S, Zhang Y, Tuohy T, Neklason D, Burt R, Guthery S, Woodward S, Jorde LB, 2011.  Maximum-likelihood estimation of recent shared ancestry (ERSA). Genome Research 2011 Feb 8.
*co-first authors.

http://genome.cshlp.org/content/early/2011/02/08/gr.115972.110.abstract

When using region_masking or other new features in ERSA 2.0, please cite:

Li H, Glusman G, Hu H, Shankaracharya, Caballero J, Hubley R, Witherspoon D., Guthery SL, Mauldin DE, Jorde LB, Hood L, Roach JC, Huff CD, 2014. Relationship Estimation from Whole-Genome Sequence Data.  PLOS Genet 10(1): e1004144. doi:10.1371/journal.pgen.1004144.

http://www.plosgenetics.org/article/fetchObject.action?uri=info%3Adoi%2F10.1371%2Fjournal.pgen.1004144&representation=PDF

Please contact chad@hufflab.org or david.witherspoon@utah.edu with questions.

** EXAMPLES **

* EXAMPLE FILES *
family.match: Pairwise IBD segment data estimated by GERMLINE 1.4 from Affymetrix SNP 6.0 genotypes collected on 169 individuals from three large Utah pedigrees of northern European descent.

ceu.match: Pairwise IBD segment data estimated by GERMLINE 1.4 from Affymetrix SNP 6.0 genotypes collected on 60 parents in the HapMap CEU sample.

pairs.txt: a limited list of pairs of individuals with IBD segment data in the file family.match.

asc.match: Pairwise IBD segment data estimated by GERMLINE 1.4 from Affymetrix SNP 6.0 genotypes for a set of individuals.

chr5.map: Genetic map information for chromosome 5 (UCSC hg18.) 
results1.out, results2.out: files containing expected output for comparison with self-test results. 

* ESTIMATE RELATIONSHIPS *
In the ersa20 directory, type the following command: 

ersa --segment_files=data/family.match --control_files=data/ceu.match --pairs=data/pairs.txt --output_file=example.out

This command uses IBD segment data in ceu.match to estimate parameters for the distributions of the number and lengths of IBD segments for unrelated individuals in the HapMap CEU sample. It then reads IBD segment data from family.match and uses that to estimate the relationships between the pairs of individuals listed in pairs.txt. The results are written to the file example.out.

To estimate relationships for all pairs of individuals with sufficient data in family.match, remove '--pairs=pairs.txt' from the command.

* ESTIMATE RELATIONSHIPS GIVEN ASCERTAINMENT *
In the ersa20 directory, type the following command: 

ersa --segment_files=asc.match --control_files=ceu.match --ascertained_chromosome=chr5 --ascertained_position=113542695 --recombination_files=chr5.map --output_file=ascexample.out

This command estimates relationships between all pairs of individuals with data in asc.match, assuming they were ascertained for the presence of a diagnostic marker on chromosome 5 at physical position 113542695, and writes the results to the file ascexample.out. 

** OUTPUT FORMAT **
The results file  has a single header line followed by rows of tab-delimited numeric data. 

The columns are:
individual_1: Identifier of first individual in pair. 
individual_2: ID of second individual. IDs are sorted alphabetically within the pair.
est_number_of_shared_ancestors: Number of ancestors (0, 1, or 2) in the most likely model.
est_degree_of_relatedness: Degree of relationship between the pair in the most likely model (no_sig_rel indicates that no significant relationship was detected between the two individuals).
CI_2p_lower: Lower bound of (1-alpha) confidence interval for 2-ancestor models.
2p_upper: Upper bound of (1-alpha) confidence interval for 2-ancestor models.
1p_lower: Lower bound of (1-alpha) confidence interval for 1-ancestor models.
1p_upper: Upper bound of (1-alpha) confidence interval for 1-ancestor models.
0p_lower: Lower bound of (1-alpha) confidence interval for 0-ancestor (e.g. direct descendant) models.
0p_upper: Upper bound of (1-alpha) confidence interval for 0-ancestor (e.g. direct descendant) models.

** USAGE NOTES **
To see the full list of options, run ersa with the parameter '-help'
