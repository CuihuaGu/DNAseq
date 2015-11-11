for i in ZS-DNA01
#for i in ZS-DNA05
do
VIRAL=AP011097.1
READ1=/mnt/Storage2/home/liyx/HCC/Data/PrivaData/Raw/Nov21_2014/CNIBRandFudanjointproject/Desheng_Kong_20130507_CNIBR_and_Fudan_joint_project_HBVintegrationsiteanalysis/F13FTSECKF2129_MUMzarX/${i}/clean_data/*_1.fq.gz
READ2=/mnt/Storage2/home/liyx/HCC/Data/PrivaData/Raw/Nov21_2014/CNIBRandFudanjointproject/Desheng_Kong_20130507_CNIBR_and_Fudan_joint_project_HBVintegrationsiteanalysis/F13FTSECKF2129_MUMzarX/${i}/clean_data/*_2.fq.gz
OUTDIR=/mnt/Storage2/home/tuser/guch/DNAseq/chimeric/${i}
sh chimericFinding2.sh ${READ1} ${READ2} ${OUTDIR} ${VIRAL} 8 0 ${i}
done
