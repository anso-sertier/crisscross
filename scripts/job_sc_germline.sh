#!/bin/sh

set -e
set -u
source $2

TMP=$1

TUSC=$OUTBAM/$TUMOR/sc_summary_cons.txt
NOSC=$OUTBAM/$NORMAL/sc_summary_cons.txt
OUTTU=$SPEPAIR/${TUMOR}_sc
OUTNO=$SPEPAIR/${NORMAL}_sc
OUTGERM=$GERMPAIR/germline_sc

mkdir $TMP

grep -v -w "side" $TUSC | awk 'BEGIN{OFS="\t";FS="\t"}{print $1,$2-1,$2}' | uniq > $TMP/tumor_sc.bed
grep -v -w "side" $NOSC | awk 'BEGIN{OFS="\t";FS="\t"}{print $1,$2-1,$2}' | uniq > $TMP/normal_sc.bed
$BEDTOOLS intersect -a $TMP/tumor_sc.bed  -b $MAPPABLE1 -wa -wb > $TMP/overlap_tumor.bed
$BEDTOOLS intersect -a $TMP/normal_sc.bed -b $MAPPABLE1 -wa -wb > $TMP/overlap_normal.bed

time $BIN/./filter_sc_germline_mappable.py -i $TUSC -j $NOSC -k $TMP/overlap_tumor.bed -l $TMP/overlap_normal.bed -t $NBCPUTU -a $OUTTU.txt -b $OUTNO.txt -c $OUTGERM.txt

SUF="_maxsc$MAXSC"

awk '$6 >= '$MAXSC' {print}' $OUTTU.txt   > $OUTTU$SUF.txt
awk '$6 >= '$MAXSC' {print}' $OUTNO.txt   > $OUTNO$SUF.txt
awk '$6 >= '$MAXSC' {print}' $OUTGERM.txt > $OUTGERM$SUF.txt

rm -fr $TMP
