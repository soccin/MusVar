#!/bin/bash

if [ "$BSUB_ARGS" == "" ]; then BSUB_ARGS=""; fi

set -eu

SDIR="$( cd "$( dirname "$0" )" && pwd )"
export PATH=$SDIR/multicall/bin:$PATH
RDIR=$SDIR

VEPVERSION=$(vep --help | fgrep ensembl-vep | awk '{print $3}')
SAREK_INPUT=$(cat out/pipeline_info/cmd.sh.log  | fgrep INPUT: | awk '{print $2}')

if [ "$VEPVERSION" != "102.0" ]; then
    echo -e "\n\n   vep not properly installed; see instructions\n"
    echo -e "VEPVERSION=[${VEPVERSION}]\n\n"
    exit 1
fi


TARGET=$(cat out/pipeline_info/cmd.sh.log | fgrep TARGET: | awk '{print $2}')
INTERVAL_BED_FILE=$(cat out/pipeline_info/cmd.sh.log | fgrep INTERVAL_BED_FILE: | awk '{print $2}')

. $RDIR/assets/Targets/$TARGET/target.resources.sh

if [ -v VARDICT_BED_FILE ]; then
    echo -e "\nUsing VarDict BED file for calling\n"
    INTERVAL_BED_FILE=$VARDICT_BED_FILE
fi

echo TARGET=$TARGET
echo INTERVAL_BED_FILE=$INTERVAL_BED_FILE

if [ "$TARGET" == "WGS_GRCm38" ]; then
  CORES=48
else
  CORES=24
fi

echo CORES=$CORES

#######################################################
#
# Get more stats from BAM
#

#ls out/preprocessing/markduplicates/*/*cram \
ls out/preprocessing/recalibrated/*/*cram \
    | xargs -n 1 bsub $BSUB_ARGS -o LSF.S_$$/ -J Stats_$$ -n $CORES -R cmorsc1 -W 24:00 \
        $SDIR/bin/getBAMStats.sh

#######################################################
#
# Call multicall/vardict and post
#

Rscript $SDIR/multicall/getSarekPairs.R \
    $SAREK_INPUT out/preprocessing/recalibrated/ \
    | xargs -n 2 \
        bsub $BSUB_ARGS -o LSF.V_$$/ -J VarD_$$ -n $CORES -W 48:00 -R cmorsc1 \
            $SDIR/VarDict/varDictPaired.sh post \
            $INTERVAL_BED_FILE

bSync VarD_$$

echo -e "\n\n============================================================"
echo -e "\n\nDone with varDictPaired.sh\n\n"

Rscript $SDIR/multicall/getSarekPairs.R $SAREK_INPUT \
    | xargs -n 3 bsub $BSUB_ARGS -o LSF.PS/ -J PS_$$ -n 5 -W 48:00 \
        $SDIR/multicall/postSarekPair.sh $TARGET out/variant_calling

bSync PS_$$

echo -e "\n\n============================================================"
echo -e "Done with postSarekPair.sh\n\n"

Rscript $SDIR/multicall/filter01.R
Rscript $SDIR/src/qcReport1.R

bSync Stats_$$

