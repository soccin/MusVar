#!/bin/bash

bsub () {
    if [ "$EXECUTOR" == "local" ]; then
        exec $@
    else
        exec bsub -o post/logs/LSF/ -J PIC.chm -n 3 -R "rusage[mem=12]" -W 359 \
        $@
    fi
}

SDIR="$( cd "$( dirname "$0" )" && pwd )"

RBASE="/rtsess01/compute/juno/bic/ROOT/rscr/references"
fasta="${RBASE}/Mus_musculus/BIC_MSK/GRCm38/Sequence/WholeGenomeFasta/GRCm38.fa"
TARGET_DIR=$SDIR/../assets/Targets/
BAM=$1
SID=$(basename $1 | sed 's/.recal.cram//')

mkdir -vp post/logs/LSF

exec >> post/logs/${SID}_hs_metrics.log
exec 2>&1

ODIR=post/metrics/collectHsMetrics
mkdir -vp $ODIR

bsub picardV2 CollectHsMetrics \
    I=$BAM \
    O=$ODIR/${SID}_hs_metrics.txt \
    R=$fasta \
    BAIT_INTERVALS=$TARGET_DIR/M-IMPACT_v2_GRCm38_baits.ilist \
    TARGET_INTERVALS=$TARGET_DIR/M-IMPACT_v2_GRCm38_targets.ilist \
    PER_TARGET_COVERAGE=$ODIR/${SID}_hs_metrics_Targets.txt
