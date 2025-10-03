#!/bin/bash

module load samtools

get_sm_tag () {
    samtools view -H $1 \
        | egrep "^@RG" \
        | head -1 \
        | tr '\t' '\n' \
        | egrep "^SM:" \
        | sed 's/^SM://'
}

SDIR="$( cd "$( dirname "$0" )" && pwd )"
VDIR=$SDIR/VarDict-1.8.3/bin

source $SDIR/genomes_BIC_MSK_GRCm38.sh

AF_THR="0.01" # minimum allele frequency

if [ "$#" != "4" ]; then
    echo -e "\n\n   usage: varDictPair.sh OUT_DIR REGIONS.bed NORMAL.bam TUMOR.bam"
    echo -e "\n"
    exit
fi

if [ "$LSB_DJOB_NUMPROC" == "" ]; then
  CORES=12
else
  CORES=$LSB_DJOB_NUMPROC
fi

# Number of parallel VarDict jobs to run
PARALLEL_JOBS=${PARALLEL_JOBS:-16}

ODIR=$1
BED=$(realpath $2)
NORMAL=$(realpath $3)
TUMOR=$(realpath $4)

NID=$(get_sm_tag $NORMAL)
TID=$(get_sm_tag $TUMOR)

EXTENSION=$($SDIR/getCommonSuffix.py $(basename $TUMOR) $(basename $NORMAL))

TTAG=$(basename ${TUMOR/$EXTENSION/})
NTAG=$(basename ${NORMAL/$EXTENSION/})

TUID=$(date +"%Y%m%d_%H%M%S")_$(uuidgen | sed 's/-.*//')
TDIR=/scratch/$USER/vardict/$TUID
mkdir -p $TDIR

trap "rm -rf $TDIR" EXIT

echo \$TDIR=$TDIR

if [[ "$TUMOR" == *.cram ]]; then
    samtools view -@ $((CORES / 2)) -t $fasta -b $NORMAL >$TDIR/$(basename ${NORMAL/.cram/.bam}) &
    samtools view -@ $((CORES / 2)) -t $fasta -b $TUMOR >$TDIR/$(basename ${TUMOR/.cram/.bam})

    wait

    NORMAL=$TDIR/$(basename ${NORMAL/.cram/.bam})
    TUMOR=$TDIR/$(basename ${TUMOR/.cram/.bam})

    samtools index -@ $((CORES / 2)) $NORMAL &
    samtools index -@ $((CORES / 2)) $TUMOR

    wait
fi

ODIR=$ODIR/variant_calling/vardict/$(basename ${BED/.bed/})/${TTAG}_vs_${NTAG}
mkdir -p $ODIR

OVCF=$ODIR/${TTAG}_vs_${NTAG}.vardict.vcf

# Split BED file by chromosome (column 1)
BEDDIR=$ODIR/bed_chunks
mkdir -p $BEDDIR

awk -v beddir="$BEDDIR" '{
    chr = $1
    count[chr]++
    file_num = int((count[chr] - 1) / 1000)
    filename = sprintf("%s/chr_%s_%03d.bed", beddir, chr, file_num)
    print > filename
}' $BED

# Get list of BED chunks (sorted by chromosome)
BED_CHUNKS=($BEDDIR/chr_*.bed)
TOTAL_CHUNKS=${#BED_CHUNKS[@]}

echo "Processing $TOTAL_CHUNKS chromosomes in parallel (max $PARALLEL_JOBS at a time)"

# Function to run VarDict on a single BED chunk
run_vardict_chunk() {
    local chunk_bed=$1
    local chunk_name=$(basename $chunk_bed .bed)
    local chunk_vcf=$ODIR/tmp/${chunk_name}.vcf.gz

    $SDIR/clipRegionsByCoverage.sh $TUMOR $chunk_bed

    timeout -k 30 3800 \
    $VDIR/VarDict \
        -th $((CORES / PARALLEL_JOBS)) \
        -G $fasta \
        -f $AF_THR \
        -N $TID \
        -b "$TUMOR|$NORMAL" \
        -c 1 -S 2 -E 3 ${chunk_bed/.bed/.clip.bed} \
        > ${chunk_vcf/.vcf.gz/.tbl}
    cat ${chunk_vcf/.vcf.gz/.tbl} \
        | $VDIR/testsomatic.R \
        | $VDIR/var2vcf_paired.pl \
             -N "$TID|$NID" \
             -f $AF_THR \
             -G $fasta \
             -b $SDIR/GRCm38.bed \
        | bgzip -c \
        > $chunk_vcf
    tabix -p vcf $chunk_vcf

}

export -f run_vardict_chunk
export VDIR TID TUMOR NORMAL CORES PARALLEL_JOBS fasta AF_THR NID SDIR TDIR ODIR

mkdir -p $ODIR/tmp

# Process chunks in parallel using GNU parallel
printf '%s\n' "${BED_CHUNKS[@]}" \
  | parallel -j $PARALLEL_JOBS --timeout 3600 --joblog $ODIR/parallel.log \
    'run_vardict_chunk {}'

# Check for failed parallel jobs
FAILED_JOBS=$(awk '$7 == -1 {print}' $ODIR/parallel.log)
if [ -n "$FAILED_JOBS" ]; then
    echo "ERROR: Some VarDict chunks failed to complete:"
    echo "$FAILED_JOBS"
    exit 1
fi

echo "All chunks completed. Merging VCFs..."

# Merge VCFs using bcftools concat
module load bcftools

# Sort chunk VCFs by chromosome order
CHUNK_VCFS=()
for chunk in "${BED_CHUNKS[@]}"; do
    chunk_name=$(basename $chunk .bed)
    chunk_vcf="$ODIR/tmp/${chunk_name}.vcf.gz"

    # Only include if file exists and has non-zero size
    if [ -s "$chunk_vcf" ]; then
        CHUNK_VCFS+=("$chunk_vcf")
    fi
done

bcftools concat -a "${CHUNK_VCFS[@]}" > $OVCF

bgzip $OVCF
tabix -p vcf ${OVCF}.gz

echo "VCF merging complete: ${OVCF}.gz"

