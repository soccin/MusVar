#!/bin/bash

set -ue

bsub () {
    if [ "${EXECUTOR:-}" == "bsub" ]; then
        exec bsub -o post/logs/LSF/ -J PIC.chm -n 3 -R "rusage[mem=12]" -W 359 \
        $@
    else
        exec $@
    fi
}

SDIR="$( cd "$( dirname "$0" )" && pwd )"
export RDIR=$SDIR/..

RBASE="/rtsess01/compute/juno/bic/ROOT/rscr/references"
fasta="${RBASE}/Mus_musculus/BIC_MSK/GRCm38/Sequence/WholeGenomeFasta/GRCm38.fa"
TARGETS=$1

TARGET_RESOURCES=$RDIR/assets/Targets/$TARGETS/target.resources.sh

if [ ! -e $TARGET_RESOURCES ]; then
    echo
    echo "ERROR: Targets resources not found: $TARGET_RESOURCES"
    echo
    ls -1d $RDIR/assets/Targets/* | perl -pe 's|.*/|      |'
    echo
    exit 1
fi

. $TARGET_RESOURCES

BAM=$2
SID=$(basename $BAM | sed 's/.recal.cram//')

mkdir -vp post/logs/LSF

exec >> post/logs/${SID}_hs_metrics.log
exec 2>&1

ODIR=post/metrics/collectHsMetrics
mkdir -vp $ODIR

bsub picardV3 CollectHsMetrics \
    COVERAGE_CAP=2500 \
    I=$BAM \
    O=$ODIR/${SID}_hs_metrics.txt \
    R=$fasta \
    BAIT_INTERVALS=$BAIT_INTERVALS \
    TARGET_INTERVALS=$TARGET_INTERVALS \
    PER_TARGET_COVERAGE=$ODIR/${SID}_hs_metrics_Targets.txt
