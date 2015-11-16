#!/usr/bin/env bash

set -eu

OUTDIR=output
INNAME=${1-zprime-hbb.hep.gz}
OUTNAME=${INNAME%.hep.gz}.root
CARD=${2-cards/delphes_tracksmear.tcl}

function check-input() {
    if [[ ! -f $INNAME ]] ; then
	echo "no input file ${INNAME}... quitting" >&2
	return 1
    fi
}

check-input $INNAME
check-input $CARD

mkdir -p $OUTDIR
OUTPATH=${OUTDIR}/${OUTNAME}
if [[ -f $OUTPATH ]] ; then
    OLDDIR=${OUTDIR}/old
    mkdir -p $OLDDIR
    mv $OUTPATH ${OLDDIR}/$(date +%F-%R)-$OUTNAME
fi

gunzip $INNAME -c | ./DelphesSTDHEP $CARD $OUTPATH
