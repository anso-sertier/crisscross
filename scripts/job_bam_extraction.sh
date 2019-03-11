#!/bin/sh

set -e
set -u
source $4

TMP=$1
BAM=$2
TYPE=$3

OUTDIR=$OUTBAM/$TYPE

mkdir $TMP

time $BIN/./bam_extract_ab_sc_check -i $BAM -r $REFERENCE -a $TMP/abnormal.bam -b $TMP/raw_sc.bam -c $TMP/raw_sc.fasta -d $TMP/raw_sc_summary.txt -e $TMP/raw_sc_partial.txt
$SAMTOOLS sort $TMP/abnormal.bam $TMP/abnormal_sort
mv $TMP/abnormal_sort.bam $TMP/raw_abnormal.bam
$SAMTOOLS index $TMP/raw_abnormal.bam
$SAMTOOLS index $TMP/raw_sc.bam

mv $TMP/raw_sc*  $TMP/raw_abnormal.bam*  $OUTDIR/.
rm -fr $TMP
