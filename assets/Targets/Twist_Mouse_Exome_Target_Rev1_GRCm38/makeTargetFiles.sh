#!/bin/bash

module load bedtools

#
# Original from
# https://www.twistbioscience.com/resources/data-files/twist-mouse-exome-panel-bed-file
#

#
# Remove random chromosomes, do not seem to make the
# GRCm38 assembly being used.
#

cat Twist_Mouse_Exome_Target_Rev1_7APR20.bed \
    | fgrep -v _ \
    | sed 's/^chr//' \
    | sort -k1,1V -k2,2n -k3,3n \
    >Twist_Mouse_Exome_Target_Rev1_7APR20_GRCm38_CLEAN.bed

#
# Add padding to the targets. This is the file that is passed
# to Sarek as the intervals file (--intervals)
#
bedtools slop \
    -i Twist_Mouse_Exome_Target_Rev1_7APR20_GRCm38_CLEAN.bed \
    -b 1000 \
    -g ~/lib/bedtools/genomes/mouse.mm10.genome \
    | bedtools merge -i - \
    > Twist_Mouse_Exome_Target_Rev1_7APR20_GRCm38_CLEAN__1000pad.bed

(
    cat GRCm38.dict; \
    cat Twist_Mouse_Exome_Target_Rev1_7APR20_GRCm38_CLEAN.bed \
        | awk '{print $1,$2-1,$3,"+","Target_"++s}' \
        | tr ' ' '\t'
) > Twist_Mouse_Exome_Target_Rev1_7APR20_GRCm38_targets.ilist

cp Twist_Mouse_Exome_Target_Rev1_7APR20_GRCm38_targets.ilist Twist_Mouse_Exome_Target_Rev1_7APR20_GRCm38_baits.ilist
