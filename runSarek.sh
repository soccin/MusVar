#!/bin/bash

OPWD=$(pwd -P)
SDIR="$( cd "$( dirname "$0" )" && pwd )"

RDIR=$(realpath $SDIR)

export NXF_SINGULARITY_CACHEDIR=/rtsess01/compute/juno/bic/ROOT/opt/singularity/cachedir_socci
export TMPDIR=/scratch/socci
export PATH=$RDIR/bin:$PATH

usage() {
    echo
    echo "  " usage: runSarek.sh [-t TARGETS] INPUT_SAREK.csv
    echo
    exit
}

TARGETDIR=M-IMPACT_v2_GRCm38

while getopts "ht:" opt; do
    case $opt in
        t) TARGETDIR=$OPTARG ;;
        h|*) usage ;;
    esac
done

shift $((OPTIND-1))

if [ "$#" -lt "1" ]; then
    usage
fi

TARGETS=$RDIR/assets/Targets/$TARGETDIR/*__1000pad.bed

if [ ! -e $TARGETS ]; then
    echo
    echo "ERROR: Targets file not found: $TARGETDIR"
    echo 
    ls -1d $RDIR/assets/Targets/* | perl -pe 's|.*/|      |'
    echo
    exit 1
fi

INPUT=$(realpath $1)

ODIR=$(pwd -P)/out

#
# Need each instance to run in its own directory
#
TUID=$(date +"%Y%m%d_%H%M%S")_$(uuidgen | sed 's/-.*//')
WDIR=run/$TUID

mkdir -p $WDIR
cd $WDIR

LOG=sarekRun.log

echo \$WDIR=$(realpath .) >$LOG
echo \$ODIR=$ODIR >>$LOG

echo nextflow run $RDIR/sarek/main.nf -ansi-log false \
    -profile singularity \
    -c $RDIR/conf/genomes_BIC_MSK_GRCm38.config \
    -c $RDIR/conf/neo.config \
    --genome null --igenomes_ignore true \
    --tools freebayes,mutect2,strelka,manta \
    --intervals $TARGETS \
    --input $INPUT \
    --outdir $ODIR \
    > $LOG
    2> ${LOG/.log/.err}

CMD_LOG=$ODIR/pipeline_info/cmd.sh.log
mkdir -p $(dirname $CMD_LOG)

GTAG=$(git --git-dir=$RDIR/.git --work-tree=$RDIR describe --long --tags --dirty="-UNCOMMITED" --always)
GURL=$(git --git-dir=$RDIR/.git --work-tree=$RDIR config --get remote.origin.url)

cat <<-END_VERSION > $CMD_LOG
DATE: $(date)
RDIR: $RDIR
GURL: $GURL
GTAG: $GTAG
PWD: $OPWD
WDIR: $WDIR

Script: $0 $*

nextflow run $RDIR/sarek/main.nf -ansi-log false \
    -profile singularity \
    -c $RDIR/conf/genomes_BIC_MSK_GRCm38.config \
    -c $RDIR/conf/neo.config \
    --genome null --igenomes_ignore true \
    --tools freebayes,mutect2,strelka,manta \
    --intervals $TARGETS \
    --input $INPUT \
    --outdir $ODIR

END_VERSION
