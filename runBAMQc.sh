#!/bin/bash

if [ "$BSUB_ARGS" == "" ]; then BSUB_ARGS=""; fi

set -eu

SDIR="$( cd "$( dirname "$0" )" && pwd )"

#######################################################
#
# Get more stats from BAM
#

ls out/preprocessing/markduplicates/*/*cram \
    | xargs -n 1 bsub $BSUB_ARGS -o LSF.S_$$/ -J Stats_$$ -n 8 -R cmorsc1 \
        $SDIR/bin/getBAMStats.sh

bSync Stats_$$

