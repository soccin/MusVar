#!/bin/bash

SDIR="$( cd "$( dirname "$0" )" && pwd )"
GENOME=/rtsess01/compute/juno/bic/ROOT/rscr/references/Mus_musculus/BIC_MSK/GRCm38/Sequence/WholeGenomeFasta/GRCm38.fa

TARGET=$(cat out/pipeline_info/cmd.sh.log | fgrep TARGET: | awk '{print $2}')

CRAM=$1

module load samtools

SID=$(samtools view -H $CRAM | fgrep @RG | head -1 | tr '\t' '\n' | fgrep LB | sed 's/..://')

ODIR=post/metrics/$SID
mkdir -p $ODIR
TDIR=tmp/bam/$SID
mkdir -p $TDIR

samtools view -@ 8 -T $GENOME -b $CRAM -o $TDIR/${SID}.bam
samtools index -@ 8 $TDIR/${SID}.bam

picardV3 \
    CollectAlignmentSummaryMetrics \
    LEVEL=null \
    LEVEL=SAMPLE \
    I=$TDIR/${SID}.bam \
    O=$ODIR/${SID}.as.txt &

picardV3 \
    CollectInsertSizeMetrics \
    LEVEL=null \
    LEVEL=SAMPLE \
    I=$TDIR/${SID}.bam \
    O=$ODIR/${SID}.ins.txt \
    H=$ODIR/${SID}.ins.pdf &

$SDIR/collectHsMetrics.sh \
    $TARGET \
    $TDIR/${SID}.bam

wait
