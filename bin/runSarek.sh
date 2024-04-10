#!/bin/bash

OPWD=$(pwd -P)
SDIR="$( cd "$( dirname "$0" )" && pwd )"

RDIR=$(realpath $SDIR/..)

export NXF_SINGULARITY_CACHEDIR=/rtsess01/compute/juno/bic/ROOT/opt/singularity/cachedir_socci
export TMPDIR=/scratch/socci
export PATH=$RDIR/bin:$PATH


if [ "$#" -lt "1" ]; then
    echo
    echo usage: runSarek.sh INPUT_SAREK.csv
    echo
    exit
fi

INPUT=$(realpath $1)

ODIR=$(pwd -P)/out

echo $INPUT > runSarek.log
echo $PID >> runSarek.log

#
# Need each instance to run in its own directory
#
TUID=$(date +"%Y%m%d_%H%M%S")_$(uuidgen | sed 's/-.*//')
WDIR=run/$TUID

mkdir -p $WDIR
cd $WDIR

LOG=runSarek.log

echo \$WDIR=$(realpath .) >$LOG
echo \$ODIR=$ODIR >>$LOG

nextflow run $RDIR/sarek/main.nf -ansi-log false \
    -profile singularity \
    -c $RDIR/assets/GenomeResources/genomes_BIC_MSK_GRCm38.config \
    --genome null --igenomes_ignore true \
    --tools freebayes,mutect2,strelka,manta \
    --intervals $RDIR/assets/Targets/M-IMPACT_v2_mm10_targets__1000pad.bed \
    --input $INPUT \
    --outdir $ODIR #\
    # > nextflow_sarek_${PID}.log \
    # 2> nextflow_sarek_${PID}.err

mkdir -p $ODIR/runlog

GTAG=$(git --git-dir=$RDIR/.git --work-tree=$RDIR describe --all --long --tags --dirty="-UNCOMMITED" --always)
GURL=$(git --git-dir=$RDIR/.git --work-tree=$RDIR config --get remote.origin.url)

# cat <<-END_VERSION > $ODIR/cmd.sh.log
# RDIR: $RDIR
# GURL: $GURL
# GTAG: $GTAG
# PWD: $OPWD
# WDIR: $WDIR

# Script: $0 $*

# END_VERSION
