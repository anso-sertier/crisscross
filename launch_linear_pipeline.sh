#!/bin/sh

set -e
set -u

FILE=$1
source $FILE

##############################################################
#### STEP 1 : extract and process abnormal pairs and soft-clip
#### from normal and tumor bam
#### each bam (3 sub-steps) can be ran in parallel
##############################################################

for TYPE in $TUMOR $NORMAL
do

	if [[ $TYPE =~ $NORMAL ]]
		then
			BAM=${BAMNORMAL}
	elif [[ $TYPE =~ $TUMOR ]]
		then
			BAM=${BAMTUMOR}
	fi

	TMP=$TMP/bam.$ECHA.$TYPE/
	$BIN/./job_bam_extraction.sh $TMP $BAM $TYPE $FILE > $LOG/extract_bam.$ECHA.$TYPE.txt

	TMP=$TMP/realign.$ECHA.$TYPE/
	$BIN/./job_sc_realignment.sh $TMP $TYPE $CPU $FILE > $LOG/sc_realignment.$ECHA.$TYPE.txt

	TMP=$TMP/clustering.$ECHA.$TYPE/
	$BIN/./job_ab_clustering.sh $TMP $TYPE $FILE > $LOG/ab_clustering.$ECHA.$TYPE.txt

done # BAM

##############################################################
#### STEP 2 : in parallel
#### Split blocs and sc pos into specific and germline status
##############################################################

TMP=$TMP/ab_germline.$ECHA.$NORMAL.$TUMOR
$BIN/./job_ab_germline.sh $TMP $FILE > $LOG/ab_germline.$ECHA.$NORMAL.$TUMOR.txt

TMP=$TMP/sc_germline.$ECHA.$NORMAL.$TUMOR
$BIN/./job_sc_germline.sh $TMP $FILE > $LOG/sc_germline.$ECHA.$NORMAL.$TUMOR.txt

#################################################
#### STEP 3: in parallel
#### Search SV candidates for each splitted file
#################################################

TMP=$TMP/sv.d$delta.${ECHA}.${NORMAL}_${TUMOR}.somatic/
LOGFILE=$LOG/sv_detect.d$delta.${ECHA}.${NORMAL}_${TUMOR}.somatic.filt
$BIN/./job_sv_detect.sh $TMP $TUMOR $NBCPUTU $LOGFILE $FILE > $LOGFILE.txt

TMP=$TMP/sv.d$delta.${ECHA}.${NORMAL}_${TUMOR}.germline/
LOGFILE=$LOG/sv_detect.d$delta.${ECHA}.${NORMAL}_${TUMOR}.germline.filt
$BIN/./job_sv_detect.sh $TMP germline 1 $LOGFILE $FILE > $LOGFILE.txt

TMP=$TMP/sv.d$delta.${ECHA}.${NORMAL}_${TUMOR}.normal/
LOGFILE=$LOG/sv_detect.d$delta.${ECHA}.${NORMAL}_${TUMOR}.normal.filt
$BIN/./job_sv_detect.sh $TMP $NORMAL 1 $LOGFILE $FILE > $LOGFILE.txt

