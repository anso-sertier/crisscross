#!/bin/sh

set -e
set -u
source $4

TMP=$1
TYPE=$2
CPU=$3

OUTDIR=$OUTBAM/$TYPE

mkdir $TMP
cp $OUTDIR/raw_sc_summary.txt $OUTDIR/raw_sc.fasta $OUTDIR/raw_sc_partial.txt $TMP/.

time $BIN/./sc_consensus.py -i $TMP/raw_sc_summary.txt -j $TMP/raw_sc_partial.txt -f $TMP/raw_sc.fasta -o $TMP/sc_summary_cons.txt -c $CPU
$BIN/./sc_cons2fastq.py -i $TMP/sc_summary_cons.txt -o $TMP/consensus.fastq -n 1 -l 15
$BWA aln -l 15 -k 1 -t $CPU $REFBWA $TMP/consensus.fastq > $TMP/consensus.sai
$BWA samse -n 10 -f $TMP/consensus.sam $REFBWA $TMP/consensus.sai $TMP/consensus.fastq
time $BIN/./sc_analyse_realignment.py -i $TMP/consensus.sam -s $TMP/sc_summary_cons.txt -r $REFERENCE -o $TMP/sc_reali.txt -m 100

mv $TMP/sc_summary_cons.txt $TMP/consensus.sam $TMP/sc_reali.txt $OUTDIR/.

rm -r $TMP
