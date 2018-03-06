#What's in RunAlleleDN.sh
# date > alleleDB.log
# for f in *.fastq.gz ;
# do echo $f >> alleleDB.log
# bash /workdir/sc2457/GM_GroSeq_AlleleDB_20170606/ENCODE_20170609/alleledb_strandSpecific_9_specifyFASTQ.sh \
# ${f:0:-9} \
# NA12878_hg19_150109 \
# /workdir/sc2457/alleleseq.gersteinlab.org/NA12878_diploid/NA12878_diploid_2015_feb5_3versions/1kgp3-svs-pass_NA12878_hg19_150109_w_transcriptome \
# /workdir/sc2457/SNP/1000genome_vol1.ftp.release.20130502/snp.calls.bed \
# $(pwd) \
# ${f:0:-9} \
# /workdir/sc2457/alleleDB/alleledb_pipeline \
# /workdir/sc2457/GM_GroSeq_AlleleDB_20170606/PIPELINE.mk \
# 0.1 \
# asb \
# 0 >> alleleDB.log 2>&1
# done


####
# move the fastq.gz files to a folder and delete it MANUALLY afterward 
mkdir tmp
for f in *_interestingHets.txt ; 
do mv ${f:0:11}.fastq.gz tmp/.
done

# copy the counts.txt files to a folder
mkdir ENCODE_counts
for f in ENCODE_*Batch/*_interestingHets.txt ;
do cp ${f:0:-19}counts.txt ENCODE_counts/.

# move the counts.txt to 0_counts_files
cd ENCODE_counts
mv *counts.txt 0_counts_files/.

# make_1_CNVfiltered_counts_bed.sh
cd 0_counts_files
mkdir ../1_CNVFiltered_counts_bed
for f in *_counts.txt ;
do echo -e \#chrm$'\t'chromStart$'\t'snppos$'\t'${f:0:-11}.winning$'\t'${f:0:-11}.SymCls$'\t'${f:0:-11}.SymPval > ../1_CNVFiltered_counts_bed/${f:0:-11}.CNVFiltered_counts.bed;
cat $f | awk 'BEGIN{OFS="\t"} (NR >1 && $18 >0.5 && $18<1.5) {print "chr"$1,$2-1, $2,$14,$15,$16}' >> ../1_CNVFiltered_counts_bed/${f:0:-11}.CNVFiltered_counts.bed;
done
#make_1 end

# make_2_CNV.PeakFiltered_counts.bed.sh
# filter counts.txt with narrow or broad peak calls
python /workdir/sc2457/GM_GroSeq_AlleleDB_20170606/ENCODE_20170609/After_AlleleDB_pipeline/pair_up_ENCODE_fastq_bed.py
# make_2 end


#make_3_CNV.Peak.BlackList.Filtered_counts_bed.sh
#then filter against ENCODE black list
cd 2_CNV.PeakFiltered_counts_bed/
mkdir ../3_CNV.Peak.BlackList.Filtered_counts_bed
for f in *Filtered_counts.bed ;
do echo -e \#chrm$'\t'chromStart$'\t'snppos$'\t'${f:0:11}.winning$'\t'${f:0:11}.SymCls$'\t'${f:0:11}.SymPval > ../3_CNV.Peak.BlackList.Filtered_counts_bed/${f:0:-19}.BlackList.Filtered_counts.bed;
intersectBed -wa -v -a $f -b /workdir/sc2457/GM_GroSeq_AlleleDB_20170606/ENCODE_20170609/ENCFF001TDO.bed >> ../3_CNV.Peak.BlackList.Filtered_counts_bed/${f:0:-19}.BlackList.Filtered_counts.bed;
done
#make_3 end

# there is NO make_4 or make_5

# make_6_CNV.Peak.BlackList.Filtered_counts_bed_in_GMRegions_ConAndDis.sh
#filter using Concordant and Discordant region
echo GMRegions_path= $1 
echo GMRegions_Name= $2
cd 3_CNV.Peak.BlackList.Filtered_counts_bed
folder=../6_CNV.Peak.BlackList.Filtered_counts_bed_in_$2
mkdir ${folder}
grep Concordant $1 > ${folder}/Concordant.bed
grep Discordant $1 > ${folder}/Discordant.bed
for f in *Filtered_counts.bed ;
do echo -e chrm$'\t'chromStart$'\t'snppos$'\t'${f:0:11}.winning$'\t'${f:0:11}.SymCls$'\t'${f:0:11}.SymPval > ${folder}/${f:0:-20}.GM.ConcordantRegions.Filtered_counts.bed;
echo -e chrm$'\t'chromStart$'\t'snppos$'\t'${f:0:11}.winning$'\t'${f:0:11}.SymCls$'\t'${f:0:11}.SymPval > ${folder}/${f:0:-20}.GM.DiscordantRegions.Filtered_counts.bed;
intersectBed -wa -wb -a $f -b ${folder}/Concordant.bed >> ${folder}/${f:0:-20}.GM.ConcordantRegions.Filtered_counts.bed;
intersectBed -wa -wb -a $f -b ${folder}/Discordant.bed >> ${folder}/${f:0:-20}.GM.DiscordantRegions.Filtered_counts.bed;
done

cd ${folder}
python /workdir/sc2457/GM_GroSeq_AlleleDB_20170606/ENCODE_20170609/After_AlleleDB_pipeline/FisherExactTest_FromGM.DiscordantRegions.Filtered_counts.bed_AsymDomminat.py
#make_6 end

## how to run make_6
# eg:
bash make_6_CNV.Peak.BlackList.Filtered_counts_bed_in_GMRegions_ConAndDis.sh /workdir/sc2457/GM_GroSeq_AlleleDB_20170606/test_in_total/GM12878_GroSeq_d150_dRegion250_withStrandSpecific_MinCount5_MaxPvalue1_ConFiltered.bed GM12878_GroSeq_d150_dRegion250_withStrandSpecific_MinCount5_MaxPvalue1_ConFiltered
bash make_6_CNV.Peak.BlackList.Filtered_counts_bed_in_GMRegions_ConAndDis.sh /workdir/sc2457/Groseq.gm12878_Allele_seq_result/StranSpecificASE/GM12878_GroSeq_d150_dRegion250_withStrandSpecific_MinCount5_MaxPvalue1.txt Old_Groseq



##Only used to examine if new ENCODE_AllelDB looks similar to old ENCODE_AlleleSeq by removing some of the filters
bash make_6_from2_CNV.Peak.Filtered_counts_bed_in_GMRegions_ConAndDis.sh /workdir/sc2457/Groseq.gm12878_Allele_seq_result/StranSpecificASE/GM12878_GroSeq_d150_dRegion250_withStrandSpecific_MinCount5_MaxPvalue1.txt
#filter using Concordant and Discordant region
echo GMRegions_path= $1 
echo GMRegions_Name= $2
cd 2_CNV.PeakFiltered_counts_bed
folder=../6_from2_CNV.Peak.Filtered_counts_bed_in_$2
mkdir ${folder}
grep Concordant $1 > ${folder}/Concordant.bed
grep Discordant $1 > ${folder}/Discordant.bed
for f in *Filtered_counts.bed ;
do echo -e chrm$'\t'chromStart$'\t'snppos$'\t'${f:0:11}.winning$'\t'${f:0:11}.SymCls$'\t'${f:0:11}.SymPval > ${folder}/${f:0:11}.GM.ConcordantRegions.Filtered_counts.bed;
echo -e chrm$'\t'chromStart$'\t'snppos$'\t'${f:0:11}.winning$'\t'${f:0:11}.SymCls$'\t'${f:0:11}.SymPval > ${folder}/${f:0:11}.GM.DiscordantRegions.Filtered_counts.bed;
intersectBed -wa -wb -a $f -b ${folder}/Concordant.bed >> ${folder}/${f:0:11}.GM.ConcordantRegions.Filtered_counts.bed;
intersectBed -wa -wb -a $f -b ${folder}/Discordant.bed >> ${folder}/${f:0:11}.GM.DiscordantRegions.Filtered_counts.bed;
done




