#!/bin/sh

FILE=$1
##############################################################
#### VARIABLE CONFIGURATION
#### set up these variables
##############################################################

## SAMPLES IDENTIFICATION : used to manage output directory tree
TUMOR="tumor_id"
NORMAL="normal_id"
PATIENT="patient_id"

## BAM FILES
BAMNORMAL="/path/to/normal/bam/file"
BAMTUMOR="/path/to/tumor/bam/file"

## WORKING PATH
OUTDIR="/writable/path/where/output/are/stored/"
TMP="/writable/path/to/tmp/directory"
LOG="/writable/path/to/log/files"

## SOFTWARE PATH
BINSCRIPTS="/path/to/pipeline/scripts/dir"
SAMTOOLS="/path/to/samtools"
BEDTOOLS="/path/to/bedtools"
BWA="/path/to/bwa"

## REFERENCE FILES
REFERENCE="/path/to/human_g1k_v37.fasta" # fasta file of reference genome
MAPPABLE1="/path/to/unmappable_human_g1k_v37.m5.e4.s01k90_24chr.bed" # bedtools of unmappable regions used to filter sv
REFBWA="/path/to/bwa-0.7.5a/human_g1k_v37.fasta" # path to bwa index
##############################################################

WORKDIR=$OUTDIR/$PATIENT
mkdir -p $WORKDIR
mkdir -p $LOG
mkdir -p $TMPDIR

FILE=$WORKDIR/config.sh



#### WRITING CONFIG FILE

echo "#!/bin/sh" > $FILE
echo ""         >> $FILE

echo "###### VARIABLE SPECIFIC TO SAMPLE ######" >> $FILE
echo ""                                          >> $FILE
echo "WORKDIR=$WORKDIR"                          >> $FILE
echo "ECHA=$PATIENT"                             >> $FILE
echo "TUMOR=$TUMOR"                              >> $FILE
echo "NORMAL=$NORMAL"                            >> $FILE
echo "LOG=$LOG"                                  >> $FILE
echo "OUTBAM=$WORKDIR/bam_extraction"            >> $FILE
FILTER=$WORKDIR/filter_germline/${NORMAL}_vs_${TUMOR}
echo "SPEPAIR=$FILTER/specific"                  >> $FILE
echo "GERMPAIR=$FILTER/germline"                 >> $FILE
echo "UNMAPPAIR=$FILTER/unmappable"              >> $FILE
SVDIR=$WORKDIR/sv_candidates/${NORMAL}_vs_${TUMOR}
echo "SVDIR=$WORKDIR/sv_candidates/${NORMAL}_vs_${TUMOR}" >> $FILE
echo ""                                          >> $FILE

OUTBAM=$WORKDIR/bam_extraction

mkdir -p $OUTBAM/$TUMOR $OUTBAM/$NORMAL
mkdir -p $FILTER/specific $FILTER/germline $FILTER/unmappable
mkdir -p $SVDIR

echo "###### BIN AND PARAMETERS VARIABLES ######" >> $FILE
echo ""                                           >> $FILE
echo "BIN=$BINSCRIPTS"     >> $FILE
echo "SAMTOOLS=$SAMTOOLS"  >> $FILE
echo "BEDTOOLS=/$BEDTOOLS" >> $FILE
echo "BWA=$BWA"            >> $FILE
echo ""                    >> $FILE
echo "CPU=1"               >> $FILE
echo "MAXSC=5"             >> $FILE
echo "delta=100"           >> $FILE
echo "" >> $FILE
echo "ABTYPE=\"intraFR intraRF intraFF intraRR interFR interRF interFF interRR\"" >> $FILE

echo "###### REFERENCES FILES ######" >> $FILE

echo "REFERENCE=$REFERENCE" >> $FILE
echo "MAPPABLE1=$MAPPABLE"  >> $FILE
echo "REFBWA=$REFBWA"       >> $FILE


chmod u+x $FILE
