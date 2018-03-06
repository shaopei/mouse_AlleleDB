date > alleleDB.log
for f in *.fastq.gz ; 
do echo $f >> alleleDB.log
bash /workdir/sc2457/GM_GroSeq_AlleleDB_20170606/ENCODE_20170609/alleledb_strandSpecific_9_specifyFASTQ.sh \
${f:0:-9} \
NA12878_hg19_150109 \
/workdir/sc2457/alleleseq.gersteinlab.org/NA12878_diploid/NA12878_diploid_2015_feb5_3versions/1kgp3-svs-pass_NA12878_hg19_150109_w_transcriptome \
/workdir/sc2457/SNP/1000genome_vol1.ftp.release.20130502/snp.calls.bed \
$(pwd) \
${f:0:-9} \
/workdir/sc2457/alleleDB/alleledb_pipeline \
/workdir/sc2457/GM_GroSeq_AlleleDB_20170606/PIPELINE.mk \
0.1 \
asb \
0 >> alleleDB.log 2>&1
done
