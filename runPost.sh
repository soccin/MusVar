#!/bin/bash

set -eu

SDIR="$( cd "$( dirname "$0" )" && pwd )"
export PATH=$SDIR/multicall/bin:$PATH

VEPVERSION=$(vep --help | fgrep ensembl-vep | awk '{print $3}')
SAREK_INPUT=$(cat out/pipeline_info/cmd.sh.log  | fgrep Script: | awk '{print $3}')

if [ "$VEPVERSION" != "102.0" ]; then
    echo -e "\n\n   vep not properly installed; see instructions\n"
    echo -e "VEPVERSION=[${VEPVERSION}]\n\n"
    exit 1
fi

Rscript MusVar/multicall/getSarekPairs.R \
    $SAREK_INPUT out/preprocessing/recalibrated/ \
    | xargs -n 2 \
        bsub -o LSF.V_$$/ -J VarD_$$ -n 16 -W 6:00 -R cmorsc1 \
            MusVar/VarDict/varDictPaired.sh post MusVar/assets/Targets/M-IMPACT_v2_mm10_targets__1000pad.bed

bSync VarD_$$

echo -e "\n\n============================================================"
echo -e "\n\nDone with varDictPaired.sh\n\n"

Rscript MusVar/multicall/getSarekPairs.R $SAREK_INPUT \
    | xargs -n 3 bsub -o LSF.PS/ -J PS_$$ -n 5 \
        ./MusVar/multicall/postSarekPair.sh out/variant_calling

bSync PS_$$

echo -e "\n\n============================================================"
echo -e "Done with postSarekPair.sh\n\n"

Rscript MusVar/multicall/filter01.R
