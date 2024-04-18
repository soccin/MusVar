#!/bin/bash

module load samtools

get_sm_tag () {
    samtools view -H $1 \
        | egrep "^@RG" \
        | head -1 \
        | tr '\t' '\n' \
        | egrep "^SM:" \
        | sed 's/^SM://'
}

SDIR="$( cd "$( dirname "$0" )" && pwd )"
VDIR=$SDIR/VarDict-1.8.3/bin

source $SDIR/genomes_BIC_MSK_GRCm38.sh

AF_THR="0.01" # minimum allele frequency

if [ "$#" != "4" ]; then
    echo -e "\n\n   usage: varDictPair.sh OUT_DIR REGIONS.bed NORMAL.bam TUMOR.bam"
    echo -e "\n"
    exit
fi

ODIR=$1
BED=$(realpath $2)
NORMAL=$(realpath $3)
TUMOR=$(realpath $4)

NID=$(get_sm_tag $NORMAL)
TID=$(get_sm_tag $TUMOR)

EXTENSION=$($SDIR/getCommonSuffix.py $(basename $TUMOR) $(basename $NORMAL))

TTAG=$(basename ${TUMOR/$EXTENSION/})
NTAG=$(basename ${NORMAL/$EXTENSION/})

TUID=$(date +"%Y%m%d_%H%M%S")_$(uuidgen | sed 's/-.*//')
TDIR=/scratch/$USER/vardict/$TUID
mkdir -p $TDIR

trap "rm -rf $TDIR" EXIT

samtools view -t $fasta -b $NORMAL >$TDIR/$(basename ${NORMAL/.cram/.bam}) &
samtools view -t $fasta -b $TUMOR >$TDIR/$(basename ${TUMOR/.cram/.bam})

wait

NORMAL=$TDIR/$(basename ${NORMAL/.cram/.bam})
TUMOR=$TDIR/$(basename ${TUMOR/.cram/.bam})

samtools index -@ 6 $NORMAL &
samtools index -@ 6 $TUMOR

wait

ODIR=$ODIR/variant_calling/vardict/${TTAG}_vs_${NTAG}
mkdir -p $ODIR

OVCF=$ODIR/${TTAG}_vs_${NTAG}.vardict.vcf

$VDIR/VarDict \
    -th 12 \
    -G $fasta \
    -f $AF_THR \
    -N $TID \
    -b "$TUMOR|$NORMAL" \
    -c 1 -S 2 -E 3 $BED \
    | $VDIR/testsomatic.R \
    | $VDIR/var2vcf_paired.pl \
         -N "$TID|$NID" \
         -f $AF_THR \
         -G $fasta \
         -b $SDIR/GRCm38.bed \
    > $OVCF

bgzip $OVCF
tabix -p vcf ${OVCF}.gz

