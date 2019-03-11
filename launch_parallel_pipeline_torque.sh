#!/bin/sh

set -e
set -u

SUBDIR="breast"
HG=b37
ROOTDIR=/home/anneso/PROJECT/structural_variants/snooper_v2.4/$SUBDIR
mkdir -p $ROOTDIR

IDS=/home/anneso/PROJECT/ICGC/BREAST_analysis/IDs_bam.txt
BAMDIR=/data-ddn/warehouse/database/sampletracker

PATIENTS="P03023 S09155 P03257 S07625 S07464 P03274 S09431" #" #P03023" # 

BAMNORMAL=$(grep "${ECHA}" ${IDS} | cut -f 3)
BAMTUMOR=$(grep "${ECHA}" ${IDS} | cut -f 5)
TYPENO=$(grep "${ECHA}" ${IDS} | cut -f 2)
TYPETU=$(grep "${ECHA}" ${IDS} | cut -f 4)
	
FILE=$ROOTDIR/$ECHA/config.$TYPENO.$TYPETU.`date +%m-%d_%H%M%S`.sh
source scripts/./make_config.sh $ROOTDIR/$ECHA $SUBDIR $FILE $TYPETU $TYPENO $ECHA $HG
source $FILE

jobpairwise=""

##############################################################
#### STEP 2 : in parallel
#### Split blocs and sc pos into specific and germline status
##############################################################

for TYPE in $TYPETU $TYPENO
do

	if [[ $TYPE =~ $TYPENO ]]
		then
			BAM=${BAMDIR}/${ECHA}/${TYPENO}/${BAMNORMAL}
			CPU=$NBCPUNO
	elif [[ $TYPE =~ $TYPETU ]]
		then
			BAM=${BAMDIR}/${ECHA}/${TYPETU}/${BAMTUMOR}
			CPU=$NBCPUTU
	fi

TMP=/tmp/bam.$ECHA.$TYPE/
LOGFILE=$LOG/extract_bam.$ECHA.$TYPE
job=`qsub -N bam.$ECHA.$TYPE -l walltime=30:00:00,nodes=1:ppn=$CPU -j oe -o $LOGFILE.oe << eof
time $BIN/./job_bam_extraction.sh $TMP $BAM $TYPE $LOGFILE.collectl $FILE
eof`

jobid=`echo $job | cut -d "." -f 1`

TMP=/tmp/realign.$ECHA.$TYPE/
LOGFILE=$LOG/sc_realignment.$ECHA.$TYPE
job=`qsub -N reali.$ECHA.$TYPE -l walltime=30:00:00,nodes=1:ppn=$CPU -W depend=afterok:$jobid -j oe -o $LOGFILE.oe << eof
time $BIN/./job_sc_realignment.sh $TMP $TYPE $CPU $LOGFILE.collectl $FILE
eof`

jobpairwise=$jobpairwise:`echo $job | cut -d "." -f 1`

TMP=/tmp/clustering.$ECHA.$TYPE/
LOGFILE=$LOG/ab_clustering.$ECHA.$TYPE
job=`qsub -N abclus.$ECHA.$TYPE -l walltime=30:00:00,nodes=1:ppn=1 -W depend=afterok:$jobid -j oe -o $LOGFILE.oe << eof
time $BIN/./job_ab_clustering.sh $TMP $TYPE $LOGFILE.collectl $FILE
eof`

jobpairwise=$jobpairwise:`echo $job | cut -d "." -f 1`

done # BAM

$BIN/./pairwise_pipeline.sh $jobpairwise $FILE

##############################################################
#### STEP 2 : in parallel
#### Split blocs and sc pos into specific and germline status
##############################################################

LOGFILE=$LOG/step3/ab_germline.$ECHA.$NORMAL.$TUMOR
TMP=/tmp/ab_germline.$ECHA.$NORMAL.$TUMOR
job=`qsub -N ab.germ.$ECHA -l walltime=100:00:00,nodes=1:ppn=8 $Wdepend -j oe -o $LOGFILE.oe << eof
time $BIN/./job_ab_germline.sh $TMP $LOGFILE.collectl $FILE
eof`

jobgermline=$jobgermline:`echo $job | cut -d "." -f 1`

TMP=/tmp/sc_germline.$ECHA.$NORMAL.$TUMOR
LOGFILE=$LOG/step3/sc_germline.$ECHA.$NORMAL.$TUMOR
job=`qsub -N sc.germ.$ECHA -l walltime=10:00:00,nodes=1:ppn=5 $Wdepend -j oe -o $LOGFILE.oe << eof
time $BIN/./job_sc_germline.sh $TMP $LOGFILE.collectl $FILE
eof`

jobgermline=$jobgermline:`echo $job | cut -d "." -f 1`


#################################################
#### STEP 3: in parallel
#### Search SV candidates for each splitted file
#################################################

LOGFILE=$LOG/svdetect/sv_detect.d$delta.${ECHA}.${NORMAL}_${TUMOR}.somatic.filt
TMP=/tmp/sv.d$delta.${ECHA}.${NORMAL}_${TUMOR}.somatic/
qsub -N sv.som.$ECHA -l walltime=100:00:00,nodes=1:ppn=8 $Wdepend -j oe -o $LOGFILE.oe << eof
time $BIN/./job_sv_detect.sh $TMP $TUMOR $NBCPUTU $LOGFILE $FILE
eof


LOGFILE=$LOG/svdetect/sv_detect.d$delta.${ECHA}.${NORMAL}_${TUMOR}.germline.filt
TMP=/tmp/sv.d$delta.${ECHA}.${NORMAL}_${TUMOR}.germline/
qsub -N sv.germ.$ECHA -l walltime=100:00:00,nodes=1:ppn=1 $Wdepend -j oe -o $LOGFILE.oe << eof
time $BIN/./job_sv_detect.sh $TMP germline 1 $LOGFILE $FILE
eof


LOGFILE=$LOG/svdetect/sv_detect.d$delta.${ECHA}.${NORMAL}_${TUMOR}.normal.filt
TMP=/tmp/sv.d$delta.${ECHA}.${NORMAL}_${TUMOR}.normal/
qsub -N sv.normal.$ECHA -l walltime=100:00:00,nodes=1:ppn=1 $Wdepend -j oe -o $LOGFILE.oe << eof
time $BIN/./job_sv_detect.sh $TMP $NORMAL 1 $LOGFILE $FILE
eof

