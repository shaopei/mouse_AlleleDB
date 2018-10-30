#cd /workdir/sc2457/SeqDepth_AlleleHMM_performance_test
#mv /workdir/sc2457/GM_GroSeq_AlleleDB_20170606/test_in_total/alleledb_strandSpecific_9.sh .
#cd /workdir/sc2457/GM_GroSeq_AlleleDB_20170606/test_in_total/
#ln -s /workdir/sc2457/SeqDepth_AlleleHMM_performance_test/alleledb_strandSpecific_9.sh .

#cd /workdir/sc2457/SeqDepth_AlleleHMM_performance_test
#mv /workdir/sc2457/GM_GroSeq_AlleleDB_20170606/test_in_total/PIPELINE_StrandSpecific.mk .
#cd /workdir/sc2457/GM_GroSeq_AlleleDB_20170606/test_in_total/
#ln -s /workdir/sc2457/SeqDepth_AlleleHMM_performance_test/PIPELINE_StrandSpecific.mk .

for file in SRR1552485_*.fastq.gz
do 
f=`echo ${file} | rev | cut -d . -f 3 |rev`
echo $f
cd /workdir/sc2457/SeqDepth_AlleleHMM_performance_test
mkdir ${f}
cd ${f}
ln -s ../${f}.fastq.gz .
ln -s ../alleledb_strandSpecific_9.sh .

bash alleledb_strandSpecific_9.sh \
${f} \
NA12878_hg19_150109_w_transcriptome \
/workdir/sc2457/alleleseq.gersteinlab.org/NA12878_diploid/NA12878_diploid_2015_feb5_3versions/1kgp3-svs-pass_NA12878_hg19_150109_w_transcriptome \
/workdir/sc2457/SNP/1000genome_vol1.ftp.release.20130502/snp.calls.bed \
/workdir/sc2457/SeqDepth_AlleleHMM_performance_test/${f} \
${f} \
/workdir/sc2457/alleleDB/alleledb_pipeline \
/workdir/sc2457/SeqDepth_AlleleHMM_performance_test/PIPELINE_StrandSpecific.mk \
0.1 \
ase \
0 > ${f}_alleledb_strandSpecific_9.log  2>&1 &

done

