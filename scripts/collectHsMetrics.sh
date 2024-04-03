#!/bin/bash

RBASE="/rtsess01/compute/juno/bic/ROOT/rscr/references"
fasta="${RBASE}/Mus_musculus/BIC_MSK/GRCm38/Sequence/WholeGenomeFasta/GRCm38.fa"
TARGET_DIR=/rtsess01/compute/juno/bic/ROOT/work/socci/Users/SawyersC/RomeroR1/NextFlow/Sarak/resources/targets
BAM=$1
SID=$(basename $1 | sed 's/.recal.cram//')

exec >> output_${SID}_hs_metrics.log
exec 2>&1

picardV2 CollectHsMetrics \
    I=$BAM \
    O=output_${SID}_hs_metrics.txt \
    R=$fasta \
    BAIT_INTERVALS=$TARGET_DIR/M-IMPACT_v2_GRCm38/M-IMPACT_v2_GRCm38_baits.ilist \
    TARGET_INTERVALS=$TARGET_DIR/M-IMPACT_v2_GRCm38/M-IMPACT_v2_GRCm38_targets.ilist \
    PER_TARGET_COVERAGE=output_${SID}_hs_metrics_Targets.txt
