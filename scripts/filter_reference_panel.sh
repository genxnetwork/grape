#!/bin/bash
for chrom in {1..22}; do
    bcftools view --min-ac 100 "/media/1000genome/minimac/ALL.chr$chrom.phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes.vcf.gz" -O z -o "chr$chrom.ac100.vcf.gz"
    echo "chr$chrom filtered"
done