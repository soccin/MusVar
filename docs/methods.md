# Mouse Variant Pipeline

## Main pipeline - nf-core/sarek

Variant calling was done using the nf-core/sarek pipeline [Garcia2020,Ewels2020]

## Genome Resources

Most of the references files were taken from the included igenomes source with the following changes. The main genome fasta file was not the igenomes one but instead the standard GRCm38 reference file from Ensembl (Release 68, patch level 0). The contigs were reordered to be in the standard ordering (see the genome diction file in the github repository for details) and the index was rebuilt with BWA version 7. A germline reference file, used by mutect2, was created by reformatting and merging the Mouse Genome Project variant files from the igenomes Ensemble resource. Finally a panel of normals reference file was also created. These custom VCF files are available at [ZENODO]



## References:

Garcia M, Juhos S, Larsson M et al. Sarek: A portable workflow for whole-genome sequencing analysis of germline and somatic variants. F1000Research 2020, 9:63 (https://doi.org/10.12688/f1000research.16665.2)

Ewels, PA, Peltzer, A, Fillinger, S et al. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol 38, 276–278 (2020). (https://doi.org/10.1038/s41587-020-0439-x)