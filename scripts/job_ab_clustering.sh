#!/bin/sh

set -e
set -u
source $3


TMP=$1
TYPE=$2

OUTDIR=$OUTBAM/$TYPE

mkdir $TMP

$BIN/./ab_clustering -i $OUTDIR/raw_abnormal.bam -a $TMP/blocs.txt -b $TMP/blocs.bed
$BEDTOOLS sort -i $TMP/blocs.bed > $OUTDIR/raw_blocs.bed
mv $TMP/blocs.txt                  $OUTDIR/raw_blocs.txt

rm -fr $TMP
