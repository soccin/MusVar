#!/bin/bash

if [ "#$" != "1" ]; then
    echo -e "\n   usage: deliver.sh /path/to/delivery/folder/r_00x\n\n"
fi

ODIR=$1

rsync -avP --exclude="markduplicates" --exclude="recalibrated" out $ODIR
rsync -avP post $ODIR


