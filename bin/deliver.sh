#!/bin/bash

SDIR=$(dirname "$(readlink -f "$0")")

if [ "$#" != "1" ]; then
    echo -e "\n   usage: deliver.sh /path/to/delivery/folder/r_00x\n"
    exit
fi

ODIR=$1

rsync -avP --exclude="markduplicates" --exclude="recalibrated" out/ $ODIR/sarek
rsync -avP post $ODIR

PROJNO=$(echo $ODIR | tr '/' '\n' | fgrep Proj_ | sed 's/Proj_//')

if [ "$PROJNO" != "" ]; then
    cat $SDIR/../assets/Email/delivery_email_template.txt \
        | perl -pe 's/{PROJNO}/'$PROJNO'/g' \
        | tee DELIVER_EMAIL_${PROJNO}.txt
fi
