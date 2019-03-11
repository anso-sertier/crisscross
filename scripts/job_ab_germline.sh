#!/bin/sh

set -e
set -u
source $2

TMP=$1

BLOCTU=$OUTBAM/$TUMOR/raw_blocs.txt
BLOCNO=$OUTBAM/$NORMAL/raw_blocs.txt
BEDBLOCTU=$OUTBAM/$TUMOR/raw_blocs.bed
BEDBLOCNO=$OUTBAM/$NORMAL/raw_blocs.bed

TUFILE=${TUMOR}_blocs.txt
NOFILE=${NORMAL}_blocs.txt

mkdir $TMP

## Mappability filtering
$BEDTOOLS intersect -a $BEDBLOCTU -b $MAPPABLE1 -u -f 0.9 > $TMP/overlap_tumapp_90.bed
$BEDTOOLS intersect -a $BEDBLOCNO -b $MAPPABLE1 -u -f 0.9 > $TMP/overlap_nomapp_90.bed

if [[ ${MAPPABLE2:+1} ]]
then
$BEDTOOLS intersect -a $BEDBLOCTU -b $MAPPABLE1 $MAPPABLE2 -u -f 0.75 > $TMP/overlap_tumapp_75.bed
$BEDTOOLS intersect -a $BEDBLOCNO -b $MAPPABLE1 $MAPPABLE2 -u -f 0.75 > $TMP/overlap_nomapp_75.bed
else
$BEDTOOLS intersect -a $BEDBLOCTU -b $MAPPABLE1 -u -f 0.75 > $TMP/overlap_tumapp_75.bed
$BEDTOOLS intersect -a $BEDBLOCNO -b $MAPPABLE1 -u -f 0.75 > $TMP/overlap_nomapp_75.bed
fi

time $BIN/./filter_ab_mappable.py -i $BLOCTU -j $TMP/overlap_tumapp_90.bed -k $TMP/overlap_tumapp_75.bed -o $UNMAPPAIR/$TUFILE -t $TMP/tu_map.txt -b $TMP/tu_map.bed
time $BIN/./filter_ab_mappable.py -i $BLOCNO -j $TMP/overlap_nomapp_90.bed -k $TMP/overlap_nomapp_75.bed -o $UNMAPPAIR/$NOFILE -t $TMP/no_map.txt -b $TMP/no_map.bed

## Germline filtering

sort -k1,1 -k2,2n $TMP/tu_map.bed > $TMP/tu_map_sort.bed
sort -k1,1 -k2,2n $TMP/no_map.bed > $TMP/no_map_sort.bed

for ab in $ABTYPE
do
	echo $ab
	grep $ab $TMP/tu_map_sort.bed > $TMP/$ab.tu.bed
	grep $ab $TMP/no_map_sort.bed > $TMP/$ab.no.bed
	$BEDTOOLS intersect -a $TMP/$ab.tu.bed -b $TMP/$ab.no.bed -wa -wb >> $TMP/overlap_tu_no.bed
done

time $BIN/./filter_ab_germline.py -i $TMP/tu_map.txt -j $TMP/no_map.txt -k $TMP/overlap_tu_no.bed -o $SPEPAIR/$TUFILE -p $SPEPAIR/$NOFILE -q $GERMPAIR/germline_blocs.txt -r $GERMPAIR/overlaps_tu_no.txt


awk '$9+$10 > 1 || $1 == "type" {print}' $SPEPAIR/$NOFILE             > $SPEPAIR/${NORMAL}_blocs.selected.txt
awk '$9+$10 > 1 || $1 == "type" {print}' $SPEPAIR/$TUFILE             > $SPEPAIR/${TUMOR}_blocs.selected.txt
awk '$9+$10 > 2 || $1 == "type" {print}' $GERMPAIR/germline_blocs.txt > $GERMPAIR/germline_blocs.selected.txt

rm -fr $TMP
