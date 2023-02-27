#!/bin/sh

for chr in `seq 1 22`; do
	echo $chr
	wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	bcftools view -i 'AF>0.001' -s "." ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz --force-samples |bcftools annotate -x ^INFO/AF -O z -o chr${chr}.vcf.gz
	rm ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
done

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz
bcftools view -i 'AF>0.001' -s "." ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz --force-samples |bcftools annotate -x ^INFO/AF -O v -o chrX.vcf.gz


files_list=""
for chr in `seq 1 22` X; do
	files_list="$files_list chr${chr}.vcf.gz"
done

bcftools concat $files_list | bcftools view -H -Ov -o 1000G.vcf
