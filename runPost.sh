#!/bin/bash

if [ "$BSUB_ARGS" == "" ]; then BSUB_ARGS=""; fi

set -eu

SDIR="$( cd "$( dirname "$0" )" && pwd )"
export PATH=$SDIR/multicall/bin:$PATH

VEPVERSION=$(vep --help | fgrep ensembl-vep | awk '{print $3}')
SAREK_INPUT=$(cat out/pipeline_info/cmd.sh.log  | fgrep INPUT: | awk '{print $2}')

if [ "$VEPVERSION" != "102.0" ]; then
    echo -e "\n\n   vep not properly installed; see instructions\n"
    echo -e "VEPVERSION=[${VEPVERSION}]\n\n"
    exit 1
fi

TARGET=$(cat out/pipeline_info/cmd.sh.log | fgrep TARGET: | awk '{print $2}')
INTERVAL_BED_FILE=$(cat out/pipeline_info/cmd.sh.log | fgrep INTERVAL_BED_FILE: | awk '{print $2}')

#######################################################
#
# Get more stats from BAM
#

ls out/preprocessing/markduplicates/*/*cram \
    | xargs -n 1 bsub $BSUB_ARGS -o LSF.S_$$/ -J Stats_$$ -n 8 -R cmorsc1 \
        $SDIR/bin/getBAMStats.sh

#######################################################
#
# Call multicall/vardict and post
#

Rscript $SDIR/multicall/getSarekPairs.R \
    $SAREK_INPUT out/preprocessing/recalibrated/ \
    | xargs -n 2 \
        bsub $BSUB_ARGS -o LSF.V_$$/ -J VarD_$$ -n 16 -W 6:00 -R cmorsc1 \
            $SDIR/VarDict/varDictPaired.sh post \
            $INTERVAL_BED_FILE

bSync VarD_$$

echo -e "\n\n============================================================"
echo -e "\n\nDone with varDictPaired.sh\n\n"

Rscript $SDIR/multicall/getSarekPairs.R $SAREK_INPUT \
    | xargs -n 3 bsub $BSUB_ARGS -o LSF.PS/ -J PS_$$ -n 5 \
        $SDIR/multicall/postSarekPair.sh $TARGET out/variant_calling

bSync PS_$$

echo -e "\n\n============================================================"
echo -e "Done with postSarekPair.sh\n\n"

Rscript $SDIR/multicall/filter01.R

bSync Stats_$$

