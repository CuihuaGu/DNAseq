# This script is written to find Chimeric reads using RNA-seq,(DNA-seq as well)
# 2015-04-09
# Cuihua Gu
##############################
## User input virables #######
##############################
LEFT_FASTQ=$1
#LEFT_FASTQ=/mnt/Storage2/home/liyx/HCC/Data/PrivaData/Raw/Nov21_2014/CNIBRandFudanjointproject/Desheng_Kong_20130507_CNIBR_and_Fudan_joint_project_Transcriptome_LncRNA/F13FTSECKF0092_HUMkdwO/CleanReads/ZS-RNA04_L1_1.fq.gz #$1
RIGHT_FASTQ=$2
#RIGHT_FASTQ=/mnt/Storage2/home/liyx/HCC/Data/PrivaData/Raw/Nov21_2014/CNIBRandFudanjointproject/Desheng_Kong_20130507_CNIBR_and_Fudan_joint_project_Transcriptome_LncRNA/F13FTSECKF0092_HUMkdwO/CleanReads/ZS-RNA04_L1_2.fq.gz #$2
OUTDIR=$3
#OUTDIR=/mnt/Storage2/home/guch/HCC/chimericPip/  #$3
VIRAL=$4
CPU=$5
VIRAL_MM=$6
PREFIX=$7

#################################
# creating output files/folders #
#################################
export PDF_DIR=$OUTDIR/pdfs && mkdir -p $PDF_DIR
READS_DIR=$OUTDIR/input_read_files && mkdir -p $READS_DIR 
GENOMIC_MAPPING_DIR=$OUTDIR/genome_mapping && mkdir -p $GENOMIC_MAPPING_DIR
#A=`basename ${LEFT_FASTQ}`
#B=`basename ${RIGHT_FASTQ}`
#PREFIX=`printf "%s\n%s\n"  $A $B | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/'`
GENOMIC_FEATURE_DIR=$OUTDIR/genomic_features && mkdir -p $GENOMIC_FEATURE_DIR
CHIMERIC_DIR=$OUTDIR/chimeric_reads && mkdir -p $CHIMERIC_DIR

echo "Sample: " $PREFIX
echo "output directory: " $3
echo "viral annotation: " $4
echo "CPU: " $5
echo "Mapping mismatch: " $6

##############################
## Common virables ###########
##############################
VIRAL_INDEX=/mnt/Storage/home/guch/lib/bowtie_index/${VIRAL}
HUMAN_INDEX=/mnt/Storage/home/guch/lib/hg19_bt_index/hg19
JUNC_INDEX=/mnt/Storage/home/guch/lib/bowtie_index/rRNA.multiT
CHIMERIC_INDEX=/mnt/Storage/home/guch/lib/bowtie_index/HBV73.HCV896.HDV105
VIRAL_LIST=/mnt/Storage/home/guch/lib/HBV73.HCV896.HDV105.list
STAR_INDEX=/mnt/Storage/home/tuser/fusl/pipipes/piPipes-master/common/hg19/STARIndex/
CURRENT_WD=`pwd`
echo `date` "====================================================START=================================================="

##############################
## 1. Pre mapping ###############
##############################
echo `date` " 1. Mapping input reads to rRNA and multiT with Bowtie"
# Due to a multi T region in HCV genome, first remove all-T sequence in input reads

bowtie2 -x ${JUNC_INDEX} -1 $LEFT_FASTQ -2 $RIGHT_FASTQ -k 1 --very-fast --no-mixed --no-discordant --un-conc ${READS_DIR}/${PREFIX}.x_rRNA_Ts.fq -p $CPU -S /dev/null 2>${READS_DIR}/${PREFIX}.rRNA_Ts.log

##############################
## 2. Custom genome pair-end mapping with bowtie :NEED TO BE CORRECTION!
##############################
# Custom genome contains:
	#2. 73 HBV genome
	#3. HCV genome
	#4. HDV genome
echo `date` "# OMITTED:2. Mapping x_rRNA pair end reads to custom genome with bowtie"
#bowtie $CHIMERIC_INDEX -1 ${READS_DIR}/${PREFIX}.x_rRNA_Ts.1.fq -2 ${READS_DIR}/${PREFIX}.x_rRNA_Ts.2.fq -a -v $VIRAL_MM -p $CPU > ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.viral 2>${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.viral.log

# Statistics
# Viral genome mapping distribution
#cut -f3 ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.viral|sort -n|uniq -c|sed 's/^ *//'|awk -F ' ' '{printf($2"\t"$1/2"\n")}'|sort -k1 > ${GENOMIC_MAPPING_DIR}/tmp
#join ${GENOMIC_MAPPING_DIR}/tmp $VIRAL_LIST |sed 's/ /\t/g'|sort -k2 -n > ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.viral_genome.distribution;rm ${GENOMIC_MAPPING_DIR}/tmp

#############################
## 3. hg19 genome mapping 
#############################
# Mapping pair end reads to genome, only keep the unmappable reads.
echo `date` "3. Mapping x_rRNA pair end reads to hg19"
xrRNA_LEFT_FQ=${READS_DIR}/${PREFIX}.x_rRNA_Ts.1.fq
xrRNA_RIGHT_FQ=${READS_DIR}/${PREFIX}.x_rRNA_Ts.2.fq

STAR \
		--runMode alignReads \
		--limitOutSAMoneReadBytes 1000000 \
		--genomeDir $STAR_INDEX \
		--readFilesIn ${xrRNA_LEFT_FQ} ${xrRNA_RIGHT_FQ} \
		--runThreadN $CPU \
		--outFilterScoreMin 0 \
		--outFilterScoreMinOverLread 0.72 \
		--outFilterMatchNmin 0 \
		--outFilterMatchNminOverLread 0.72 \
		--outFilterMultimapScoreRange 1 \
		--outFilterMultimapNmax -1 \
		--outFilterMismatchNmax 10 \
		--outFilterMismatchNoverLmax 0.05 \
		--alignIntronMax 0 \
		--alignIntronMin 21 \
		--outFilterIntronMotifs RemoveNoncanonicalUnannotated \
		--genomeLoad NoSharedMemory \
		--outFileNamePrefix $GENOMIC_MAPPING_DIR/${PREFIX}.x_rRNA_Ts.hg19. \
		--outSAMunmapped None \
		--outReadsUnmapped Fastx \
		--outSJfilterReads Unique \
		--seedSearchStartLmax 20 \
		--seedSearchStartLmaxOverLread 1.0 \
		--chimSegmentMin 0 2>&1 1> $GENOMIC_MAPPING_DIR/${PREFIX}.x_rRNA_Ts.${GENOME}.STAR.log

## processing sam files
samtools view -bS ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA_Ts.hg19.Aligned.out.sam > ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA_Ts.hg19.Aligned.out.bam
/mnt/Storage/home/tuser/fusl/pipipes/piPipes-master/bin/samtools sort -@ $CPU ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA_Ts.hg19.Aligned.out.bam ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA_Ts.hg19.sorted
samtools index ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA_Ts.hg19.sorted.bam
rm -rf  ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA_Ts.hg19.Aligned.out.bam
rm -rf ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA_Ts.hg19.Aligned.out.sam

##############################
## 4. Custom genome pair-end mapping with bowtie AGAIN
##############################
echo `date` "4. Mapping x_rRNA.x_hg19 pair end reads to custom genome with bowtie"
bowtie ${CHIMERIC_INDEX} -1 ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA_Ts.hg19.Unmapped.out.mate1 -2 ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA_Ts.hg19.Unmapped.out.mate2 -a -v $VIRAL_MM -p $CPU > ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA_Ts.x_hg19.viral 2>${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA_Ts.x_hg19.viral.log

# Statistics, Viral genome mapping distribution
cut -f4 ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA_Ts.x_hg19.viral|sort -n|uniq -c|sed 's/^ *//'|awk -F ' ' '{printf($2"\t"$1/2"\n")}'|sort -k1 > ${GENOMIC_MAPPING_DIR}/tmp
join ${GENOMIC_MAPPING_DIR}/tmp $VIRAL_LIST |sed 's/ /\t/g'|sort -k2 -n > ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA_Ts.x_hg19.viral_genome.distribution;rm ${GENOMIC_MAPPING_DIR}/tmp

##############################
## 5. Chimeric transcript finding
##############################
echo `date` "5. Find chimeric transcript"
cd ${CHIMERIC_DIR}
touch .viral.ID2 .viral.ID1 .hg.ID1 .hg.ID2 viral2.hg1 viral1.hg2

R1=${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA_Ts.hg19.Unmapped.out.mate1
R2=${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA_Ts.hg19.Unmapped.out.mate2
cat $R1 $R2 > .R12
bowtie ${VIRAL_INDEX} -q .R12 -a -v $VIRAL_MM -p $CPU > ${CHIMERIC_DIR}/${PREFIX}.x_rRNA_Ts.x_hg19.viral_chimeric 2>${CHIMERIC_DIR}/${PREFIX}.x_rRNA_Ts.x_hg19.viral_chimeric.log
bowtie ${HUMAN_INDEX} -q .R12 -a -v $VIRAL_MM -p $CPU > ${CHIMERIC_DIR}/${PREFIX}.x_rRNA_Ts.x_hg19.hg19_chimeric 2>${CHIMERIC_DIR}/${PREFIX}.x_rRNA_Ts.x_hg19.hg19_chimeric.log

awk -F '\t' '{split($1,a,"/");if(a[2]=="1"){printf(a[1]"/2\n") > ".viral.ID2"}else{printf(a[1]"/1\n") > ".viral.ID1"}}' ${CHIMERIC_DIR}/${PREFIX}.x_rRNA_Ts.x_hg19.viral_chimeric
awk -F '\t' '{split($1,a,"/");if(a[2]=="1"){printf(a[1]"/1\n") > ".hg.ID1"}else{printf(a[1]"/2\n") > ".hg.ID2"}}' ${CHIMERIC_DIR}/${PREFIX}.x_rRNA_Ts.x_hg19.hg19_chimeric
sort .hg.ID1 > hg.ID1;sort .hg.ID2 > hg.ID2;sort .viral.ID1 > viral.ID1;sort .viral.ID2 > viral.ID2
comm -12 hg.ID1 viral.ID1 > hg1.viral2
comm -12 hg.ID2 viral.ID2 > hg2.viral1
rm .R12 .hg.ID1 .hg.ID2 .viral.ID1 .viral.ID2 hg.ID1 hg.ID2 viral.ID1 viral.ID2

weedLines -invert hg1.viral2 ${CHIMERIC_DIR}/${PREFIX}.x_rRNA_Ts.x_hg19.hg19_chimeric hg1
weedLines -invert hg2.viral1 ${CHIMERIC_DIR}/${PREFIX}.x_rRNA_Ts.x_hg19.hg19_chimeric hg2
cat hg1.viral2 hg2.viral1 > m
awk -F '/' '{if($2=="1"){printf($1"/2\n") > "viral2.hg1"}else{printf($1"/1\n") > "viral1.hg2"}}' m;rm m
weedLines -invert viral1.hg2 ${CHIMERIC_DIR}/${PREFIX}.x_rRNA_Ts.x_hg19.viral_chimeric viral1
weedLines -invert viral2.hg1 ${CHIMERIC_DIR}/${PREFIX}.x_rRNA_Ts.x_hg19.viral_chimeric viral2

awk -F '/2' '{printf($1"\t"$2"\n")}' viral2|sort -k1 > viral2.tmp
awk -F '/2' '{printf($1"\t"$2"\n")}' hg2 |sort -k1> hg2.tmp
awk -F '/1' '{printf($1"\t"$2"\n")}' viral1|sort -k1 > viral1.tmp
awk -F '/1' '{printf($1"\t"$2"\n")}' hg1|sort -k1 > hg1.tmp

join viral1.tmp hg2.tmp|sed 's/ /\t/g'|cut -f1,3,4,5,11,12,13 > ${PREFIX}.x_rRNA_Ts.x_hg19.${VIRAL}_hg19_chimeric
join viral2.tmp hg1.tmp|sed 's/ /\t/g'|cut -f1,3,4,5,11,12,13 > ${PREFIX}.x_rRNA_Ts.x_hg19.hg19_${VIRAL}_chimeric

UNIQUE_CHIMERIC_PAIR=`wc -l hg*.viral*|tail -n1|awk '{print $1}'`
ALL_CHIMERIC_HG=`wc -l hg*.tmp|tail -n1|awk '{print $1}'`
ALL_CHIMERIC_VIRAL=`wc -l viral*.tmp|tail -n1|awk '{print $1}'`
echo "#unique chimeric pairs: " ${UNIQUE_CHIMERIC_PAIR} >> ${PREFIX}.x_rRNA_Ts.x_hg19.chimeric.log
echo "#all chimeric reads in hg19: " ${ALL_CHIMERIC_HG} >> ${PREFIX}.x_rRNA_Ts.x_hg19.chimeric.log
echo "#all chimeric reads in viral: " ${ALL_CHIMERIC_VIRAL} >> ${PREFIX}.x_rRNA_Ts.x_hg19.chimeric.log
echo "#Multiple reads ID in hg19:" >> ${PREFIX}.x_rRNA_Ts.x_hg19.chimeric.log
cut -f1 hg1 hg2|sort |uniq -d >> ${PREFIX}.x_rRNA_Ts.x_hg19.chimeric.log
echo "#Multiple reads ID in viral:" >> ${PREFIX}.x_rRNA_Ts.x_hg19.chimeric.log
cut -f1 viral1 viral2|sort |uniq -d >> ${PREFIX}.x_rRNA_Ts.x_hg19.chimeric.log

rm hg1 hg2 hg1.tmp hg2.tmp viral1 viral2 viral1.tmp viral2.tmp viral1.hg2 viral2.hg1 hg1.viral2 hg2.viral1
rm ${CHIMERIC_DIR}/${PREFIX}.x_rRNA_Ts.x_hg19.viral_chimeric*
rm ${CHIMERIC_DIR}/${PREFIX}.x_rRNA_Ts.x_hg19.hg19_chimeric*

cd ${CURRENT_WD}

##############################
## 6. Genomic features of chimeric reads
##############################
echo `date` "6. Genomic features of chimeric reads"
cd ${GENOMIC_FEATURE_DIR}

cat ${CHIMERIC_DIR}/${PREFIX}.x_rRNA_Ts.x_hg19.${VIRAL}_hg19_chimeric ${CHIMERIC_DIR}/${PREFIX}.x_rRNA_Ts.x_hg19.hg19_${VIRAL}_chimeric > ${CHIMERIC_DIR}/${PREFIX}.x_rRNA_Ts.x_hg19.chimeric
awk '{printf($6"\t"$7"\t"$7+90"\t"$1"\t0\t"$5"\t0\n")}' ${CHIMERIC_DIR}/${PREFIX}.x_rRNA_Ts.x_hg19.chimeric > ${GENOMIC_FEATURE_DIR}/${PREFIX}.bed2
sh /mnt/Storage2/home/guch/HCC/chimericPip/Code/genomic_feature.sh ${GENOMIC_FEATURE_DIR}/${PREFIX}.bed2 $PREFIX

cd ${CURRENT_WD}


echo `date` "====================================================FINISHED=================================================="


