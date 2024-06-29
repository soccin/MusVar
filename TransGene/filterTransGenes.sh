#!/bin/bash

SDIR="$( cd "$( dirname "$0" )" && pwd )"

if [ "$#" != "1" ]; then
    echo -e "\n\tusage:filterTransGenes.sh FILE.sam\n\n"
    exit
fi

SAM=$1

SM=$(head -1000 $SAM | fgrep "@RG" | head -1 | tr '\t' '\n' | fgrep "SM" | sed 's/SM://' | sed 's/^s_//')

if [ "$SM" == "" ]; then
    echo -e "\n\tVALUE ERROR: Can not find SM Tag\n\n"
    exit -1
fi



ODIR=fastqTg/Sample_${SM}
mkdir -p $ODIR

python3 $SDIR/filterTransGenes.py $SAM \
    | picardV2 SamToFastq I=/dev/stdin \
        F=$ODIR/${SM}_R1_001.fastq.gz \
        F2=$ODIR/${SM}_R2_001.fastq.gz \
        FU=$ODIR/${SM}_UP_001.fastq.gz

