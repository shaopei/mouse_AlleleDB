ln -s /workdir/sc2457/GM_GroSeq_AlleleDB_20170606/test_in_total/allelicbias-NA12878_hg19_150109_w_transcriptome-SRR1552485_total/1-alignment-SRR1552485_total/SRR1552485.scythed.sickled.fastq.gz SRR1552485_total.fastq.gz

zcat SRR1552485_total.fastq.gz > SRR1552485_total.fastq
fastq-sample -p 0.5 -o SRR1552485_sub2 SRR1552485_total.fastq   -s 100    #18
fastq-sample -p 0.5 -o SRR1552485_sub4 SRR1552485_sub2.fastq  -s 100  #9
fastq-sample -p 0.5 -o SRR1552485_sub8 SRR1552485_sub4.fastq   -s 100  #4.5
fastq-sample -p 0.5 -o SRR1552485_sub16 SRR1552485_sub8.fastq  -s 100  #2.25
wc -l *.fastq > 01_AlleleDB_Run.log &

for f in SRR1552485_sub*.fastq
do gzip $f &
done
wait
rm SRR1552485_total.fastq
