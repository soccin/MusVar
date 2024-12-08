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


#
# Save full command line to print in log
#
CMDLINE="$@"

TARGET=M-IMPACT_v2_GRCm38

while getopts "ht:" opt; do
    case $opt in
        t) TARGET=$OPTARG ;;
        h|*) usage ;;
    esac
done

shift $((OPTIND-1))

if [ "$#" -lt "1" ]; then
    usage
fi

TARGET_RESOURCES=$RDIR/assets/Targets/$TARGET/target.resources.sh

if [ ! -e $TARGET_RESOURCES ]; then
    echo
    echo "ERROR: Targets resources not found: $TARGET_RESOURCES"
    echo 
    ls -1d $RDIR/assets/Targets/* | perl -pe 's|.*/|      |'
    echo
    exit 1
fi

. $TARGET_RESOURCES

echo
echo "  "Using target resources: ${TARGET}
echo

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

#
# If the script is running in a terminal, then set ANSI_LOG to true
#
case $(ps -o stat= -p $$) in
  *+*) ANSI_LOG="true" ;;
  *) ANSI_LOG="false" ;;
esac

nextflow run $RDIR/sarek/main.nf -ansi-log $ANSI_LOG \
    -resume \
    -profile singularity \
    -c $RDIR/conf/genomes_BIC_MSK_GRCm38.config \
    -c $RDIR/conf/${TARGET}.config \
    -c $RDIR/conf/neo.config \
    --genome null --igenomes_ignore true \
    --tools freebayes,mutect2,strelka,manta \
    --intervals $INTERVAL_BED_FILE \
    --input $INPUT \
    --outdir $ODIR \
    2> ${LOG/.log/.err} \
    | tee $LOG

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

Script: $0 $CMDLINE
TARGET: $TARGET
INTERVAL_BED_FILE: $INTERVAL_BED_FILE
INPUT: $INPUT
ODIR: $ODIR

nextflow run $RDIR/sarek/main.nf -ansi-log $ANSI_LOG \
    -resume \
    -profile singularity \
    -c $RDIR/conf/genomes_BIC_MSK_GRCm38.config \
    -c $RDIR/conf/${TARGET}.config \
    -c $RDIR/conf/neo.config \
    --genome null --igenomes_ignore true \
    --tools freebayes,mutect2,strelka,manta \
    --intervals $INTERVAL_BED_FILE \
    --input $INPUT \
    --outdir $ODIR

END_VERSION
