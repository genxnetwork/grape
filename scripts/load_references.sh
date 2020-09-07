#wget -O- ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | \
            #  gzip -d > {{output.fasta}}
#samtools faidx {{output.fasta}}
            #wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr{{9..22}}_GRCh38.genotypes.20170504.vcf.gz{{,.tbi}}

(bcftools view --no-version -h ALL.chr${c}_GRCh38.genotypes.20170504.vcf.gz | \
   grep -v "^##contig=<ID=[GNh]" | sed 's/^##contig=<ID=MT/##contig=<ID=chrM/;s/^##contig=<ID=\([0-9XY]\)/##contig=<ID=chr\1/'; \
   bcftools view --no-version -H -c 2 ALL.chr${c}_GRCh38.genotypes.20170504.vcf.gz | sed 's/^/chr/') | \
bcftools norm --no-version -Ou -m -any | \
bcftools norm --no-version -Ob -o ALL.chr${c}_GRCh38.genotypes.20170504.bcf -d none -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna && \
bcftools index -f ALL.chr${c}_GRCh38.genotypes.20170504.bcf
