#for i in ZS-DNA01
for i in ZS-DNA02 ZS-DNA03 ZS-DNA04 ZS-DNA15 ZS-DNA06 ZS-DNA07 ZS-DNA08 ZS-DNA13 ZS-DNA14 ZS-DNA16 ZS-DNA17 ZS-DNA18 ZS-DNA19 ZS-DNA20 ZS-DNA21 ZS-DNA22 ZS-DNA29 ZS-DNA28
do
VIRAL=AP011098.1
READ1=/mnt/Storage2/home/liyx/HCC/Data/PrivaData/Raw/Nov21_2014/CNIBRandFudanjointproject/Desheng_Kong_20130507_CNIBR_and_Fudan_joint_project_HBVintegrationsiteanalysis/F13FTSECKF2129_MUMzarX/${i}/clean_data/*_1.fq.gz
READ2=/mnt/Storage2/home/liyx/HCC/Data/PrivaData/Raw/Nov21_2014/CNIBRandFudanjointproject/Desheng_Kong_20130507_CNIBR_and_Fudan_joint_project_HBVintegrationsiteanalysis/F13FTSECKF2129_MUMzarX/${i}/clean_data/*_2.fq.gz
OUTDIR=/mnt/Storage2/home/tuser/guch/DNAseq/chimeric/${i}
sh chimericFinding.sh ${READ1} ${READ2} ${OUTDIR} ${VIRAL} 8 0 ${i}
done
