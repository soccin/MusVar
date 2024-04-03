#!/bin/bash

#
# This script is used to create the `germline_resource` file.
#
#
#

bcftools () { 
    singularity exec -B /rtsess01 /home/socci/ROOT/opt/singularity/cachedir_nxf/depot.galaxyproject.org-singularity-bcftools-1.17--haef29d1_0.img bcftools "$@"; 
}
module load samtools/1.19.2


bcftools +fill-tags "/rtsess01/compute/juno/bic/ROOT/rscr/references/Mus_musculus/BIC_MSK/GRCm38/MouseGenomeProject/mgp.v5.merged.snps_all.dbSNP142.vcf.gz" -o mgp.v5.merged.snps_all.dbSNP142.vcf -- -t AN,AC,AF

cat mgp.v5.merged.snps_all.dbSNP142.vcf | sed 's/DP=.*;AN/AN/' | cut -f-8 >mgp.v5.merged.snps_all.dbSNP142.AF_InfoOnly.vcf
bgzip mgp.v5.merged.snps_all.dbSNP142.AF_InfoOnly.vcf
tabix -p vcf mgp.v5.merged.snps_all.dbSNP142.AF_InfoOnly.vcf.gz

bcftools +fill-tags "/rtsess01/compute/juno/bic/ROOT/rscr/references/Mus_musculus/BIC_MSK/GRCm38/MouseGenomeProject/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz" -o mgp.v5.merged.indels.dbSNP142.vcf -- -t AN,AC,AF &

cat mgp.v5.merged.indels.dbSNP142.vcf | sed 's/INDEL;DP=.*;AN/AN/' | cut -f-8 >mgp.v5.merged.indels.dbSNP142.AF_InfoOnly.vcf
bgzip mgp.v5.merged.indels.dbSNP142.AF_InfoOnly.vcf

tabix -p vcf mgp.v5.merged.indels.dbSNP142.vcf.gz

bcftools concat -a \
    mgp.v5.merged.snps_all.dbSNP142.AF_InfoOnly.vcf.gz \
    mgp.v5.merged.indels.dbSNP142.AF_InfoOnly.vcf.gz \
    | bgzip -c - >mgp.v5.merged.ALL.dbSNP142.AF_InfoOnly.vcf.gz

tabix mgp.v5.merged.ALL.dbSNP142.AF_InfoOnly.vcf.gz

rm mgp.v5.merged.snps_all.dbSNP142.vcf*
rm mgp.v5.merged.indels.dbSNP142.vcf*
rm mgp.v5.merged.snps_all.dbSNP142.AF_InfoOnly.vcf.gz
rm mgp.v5.merged.indels.dbSNP142.AF_InfoOnly.vcf.gz*