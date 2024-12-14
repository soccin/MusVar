#!/bin/bash

# BSUB: -o LSF/ -J C2B -n 18 -W 12:00 -R cmorsc1

if [ "$#" != "1" ]; then
    echo -e "\n\tusage: sarekCramToBam.sh FILE.cram\n"
    exit
fi

SDIR=$(dirname "$(readlink -f "$0")")

module load samtools

CRAM=$1

if [[ "$CRAM" =~ /preprocessing/ ]]; then
    ODIR=$(echo $CRAM | perl -pe 's|/preprocessing/.*||')/bam
else
    ODIR=$(dirname $CRAM)
fi

SM=$(
    samtools view -H $CRAM \
    | egrep "^@RG" \
    | head -1 \
    | tr '\t' '\n' \
    | fgrep SM: \
    | sed 's/SM://'
    )

ODIR=$ODIR/$SM
mkdir -p $ODIR

GENOME=$($SDIR/getGenomeBuildBAM.sh $CRAM)

case $GENOME in

    b37)
    GENOME_FILE=/juno/bic/depot/assemblies/H.sapiens/b37/b37.fasta
    ;;

    GRC_m38)
    GENOME_FILE=/rtsess01/compute/juno/bic/ROOT/rscr/references/Mus_musculus/BIC_MSK/GRCm38/Sequence/WholeGenomeFasta/GRCm38.fa
    ;;

    *)
    echo -e "\n\nUNKNOWN GENOME=[${GENOME}]\n\n"
    exit
    ;;

esac

samtools view -@ 16 -T $GENOME_FILE -b $CRAM -o $ODIR/${SM}.smap.bam
samtools index -@ 16 $ODIR/${SM}.smap.bam

