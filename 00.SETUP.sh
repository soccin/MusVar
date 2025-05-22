#!/bin/bash

SDIR=$(dirname "$(readlink -f "$0")")

cd $SDIR/bin
curl -s https://get.nextflow.io | bash
cd $SDIR/multicall/resources/vep/
. CMDS.INSTALL
cd $SDIR
