#!/bin/sh

set -e
set -u
source $5

TMP=$1
NAME=$2
NBPROC=$3
LOGFILE=$4


if [[ $NAME == "germline" ]]
then
	SCREALI=$OUTBAM/$NORMAL/sc_reali.txt
	BLOCS=$GERMPAIR/${NAME}_blocs.selected.txt
	SCPOS=$GERMPAIR/${NAME}_sc_maxsc$MAXSC.txt
else
	SCREALI=$OUTBAM/${NAME}/sc_reali.txt
	BLOCS=$SPEPAIR/${NAME}_blocs.selected.txt
	SCPOS=$SPEPAIR/${NAME}_sc_maxsc$MAXSC.txt
fi

LINKS=$SVDIR/${NAME}_links_d${delta}.txt
FINAL=$SVDIR/${NAME}_sv_d${delta}.txt
FILTERED=$SVDIR/${NAME}_sv_filtered_d${delta}.txt

mkdir $TMP
cp $BLOCS $TMP/blocs.txt

grep -v -w "side" $SCPOS | awk 'BEGIN{OFS="\t";FS="\t"}{print $1,$2-1,$2,$3}' | $BEDTOOLS sort -i | uniq >  $TMP/scpos.bed
grep -P "FR" $TMP/blocs.txt | awk 'BEGIN{OFS="\t"}{print $3,$5-'$delta',$5+'$delta',$1$2"a\n"$6,$7-'$delta',$7+'$delta',$1$2"b"}' >  $TMP/blocs.bed
grep -P "RR" $TMP/blocs.txt | awk 'BEGIN{OFS="\t"}{print $3,$4-'$delta',$4+'$delta',$1$2"a\n"$6,$7-'$delta',$7+'$delta',$1$2"b"}' >> $TMP/blocs.bed
grep -P "FF" $TMP/blocs.txt | awk 'BEGIN{OFS="\t"}{print $3,$5-'$delta',$5+'$delta',$1$2"a\n"$6,$8-'$delta',$8+'$delta',$1$2"b"}' >> $TMP/blocs.bed
grep -P "RF" $TMP/blocs.txt | awk 'BEGIN{OFS="\t"}{print $3,$4-'$delta',$4+'$delta',$1$2"a\n"$6,$8-'$delta',$8+'$delta',$1$2"b"}' >> $TMP/blocs.bed
awk 'BEGIN{OFS="\t"}{if($2 < 0) {$2 = 0} ; print}' $TMP/blocs.bed | $BEDTOOLS sort -i  > $TMP/blocs_sort.bed
$BEDTOOLS intersect -a $TMP/blocs_sort.bed -b $TMP/scpos.bed -wa -wb -sorted > $TMP/overlap_d$delta.bed

time $BIN/./search_links.py -i $TMP/overlap_d$delta.bed -c $TMP/blocs.txt -s $SCPOS -a $SCREALI -r $REFERENCE -o $LINKS -l $LOGFILE.log -d $delta -t $NBPROC

grep "intraFF" $LINKS | awk 'BEGIN{OFS="\t"}{print $3,$8-'$delta',$8+'$delta',$1$2"a",$11,"\n"$5,$10-'$delta',$10+'$delta',$1$2"b",$11}' | $BEDTOOLS sort -i > $TMP/oriFF.bed
grep "intraRR" $LINKS | awk 'BEGIN{OFS="\t"}{print $3,$7-'$delta',$7+'$delta',$1$2"a",$11,"\n"$5, $9-'$delta', $9+'$delta',$1$2"b",$11}' | $BEDTOOLS sort -i > $TMP/oriRR.bed
$BEDTOOLS intersect -a $TMP/oriFF.bed -b $TMP/oriRR.bed -wa -wb -sorted > $TMP/overlap.bed

if [[ $NAME == "germline" ]]
then
$BIN/./combine_links_filtered.py -i $LINKS -j $TMP/overlap.bed -o $FINAL -f $FILTERED -g
else
$BIN/./combine_links_filtered.py -i $LINKS -j $TMP/overlap.bed -o $FINAL -f $FILTERED
fi


rm -fr $TMP


### remove intermediate file no more needed
rm -f $OUTBAM/$TUMOR/raw_sc_summary.txt  $OUTBAM/$TUMOR/raw_sc.fasta
rm -f $OUTBAM/$NORMAL/raw_sc_summary.txt $OUTBAM/$NORMAL/raw_sc.fasta



