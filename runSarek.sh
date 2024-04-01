#!/bin/bash
SDIR="$( cd "$( dirname "$0" )" && pwd )"

. $SDIR/SETENV

PROJECT_ID=Proj_15301

INPUT=$(realpath $1)
PID=$(basename ${INPUT/_sarek_input*/} | sed 's/^pt_//')

echo $INPUT > ${PID}_runSarek.log
echo $PID >> ${PID}_runSarek.log

TDIR=$SDIR/scratch/$PID
mkdir -vp $TDIR >> ${PID}_runSarek.log

cd $TDIR
ODIR=$SDIR/out/${PROJECT_ID}/somatic/$PID

$SDIR/nextflow run $SDIR/sarek/main.nf -ansi-log false \
    -profile singularity \
    -c $SDIR/genomes_BIC_MSK_GRCm38.config \
    --genome null --igenomes_ignore true \
    --tools freebayes,mutect2,strelka,manta \
    --intervals $SDIR/M-IMPACT_v2_mm10_targets__1000pad.bed \
    --input $INPUT \
    --outdir $ODIR #\
    # > nextflow_sarek_${PID}.log \
    # 2> nextflow_sarek_${PID}.err

# cat <<-END_VERSION > $ODIR/cmd.sh.log
# Script: $(realpath $0)
# $SDIR/nextflow run $SDIR/sarek/main.nf -ansi-log false \
#     -profile singularity \
#     -c $SDIR/genomes_BIC_MSK_GRCm38.config \
#     --genome null --igenomes_ignore true \
#     --tools freebayes,mutect2,strelka,manta \
#     --intervals $SDIR/M-IMPACT_v2_mm10_targets__1000pad.bed \
#     --input $INPUT \
#     --outdir $ODIR
# END_VERSION
