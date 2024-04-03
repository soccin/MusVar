# Mouse Variant Pipeline

## Main Pipeline - nf-core/sarek

Variant calling was performed using the nf-core/sarek pipeline [Garcia2020, Ewels2020]. The pipeline was executed in somatic mode, utilizing the following callers: Mutect2, Strelka, and FreeBayes. Version 3.4.0 was used with one modification to the code: the options to run Mutect2 were altered to disable the calling of multiple SNP sites (`--max-mnp-distance 0`). This adjustment facilitated the creation of the panel of normals using the GATK toolkit, which requires this mode. Additionally, the genomics resources were updated and modified as detailed below.

The output from Sarek was then post-processed to merge the calls from the various algorithms into one table. Two tables were generated: one at a higher stringency, which was the union of all events marked PASS by any caller, and the second at a lowered stringency to recover clustered events by allowing Mutect2 calls that were flagged as multiallelic or clustered and Strelka events flagged as low-evidence. The script used for the filtering, along with other accessory scripts, are available here: [https://github.com/soccin/MusVar](https://github.com/soccin/MusVar).

## Genomic Resources

Most of the reference files were sourced from the included igenomes with the following modifications: The main genome FASTA file was not the one from igenomes but instead the standard GRCm38 reference file from Ensembl (Release 68, patch level 0). The contigs were reordered to standard ordering (see the genome dictionary file in the GitHub repository for details), and the index was rebuilt with BWA version 7. A germline reference file, used by Mutect2, was created by reformatting and merging the Mouse Genome Project variant files from the igenomes Ensembl resource. Finally, a panel of normals reference file was also created. These custom VCF files are available at [https://zenodo.org/records/10914483](https://zenodo.org/records/10914483).

## References

Garcia M, Juhos S, Larsson M, et al. Sarek: A portable workflow for whole-genome sequencing analysis of germline and somatic variants. F1000Research 2020, 9:63 ([https://doi.org/10.12688/f1000research.16665.2](https://doi.org/10.12688/f1000research.16665.2))

Ewels, PA, Peltzer, A, Fillinger, S, et al. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol 38, 276–278 (2020). ([https://doi.org/10.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x))
