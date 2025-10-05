#!/bin/bash

SDIR="$( cd "$( dirname "$0" )" && pwd )"
GENOME=/rtsess01/compute/juno/bic/ROOT/rscr/references/Mus_musculus/BIC_MSK/GRCm38/Sequence/WholeGenomeFasta/GRCm38.fa

TARGET=$(cat out/pipeline_info/cmd.sh.log | fgrep TARGET: | awk '{print $2}')

CRAM=$1

module load samtools

SID=$(samtools view -H $CRAM | fgrep @RG | head -1 | tr '\t' '\n' | fgrep LB | sed 's/..://')

ODIR=post/metrics/$SID
mkdir -p $ODIR

TUID=$(date +"%Y%m%d_%H%M%S")_$(uuidgen | sed 's/-.*//')
TDIR=/scratch/$USER/MusVar/$TUID
mkdir -p $TDIR

trap "rm -rf $TDIR" EXIT

CORES=${LSB_DJOB_NUMPROC:-8}
echo \$CORES=$CORES

samtools view -@ $CORES -T $GENOME -b $CRAM -o $TDIR/${SID}.bam
samtools index -@ $CORES $TDIR/${SID}.bam

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
