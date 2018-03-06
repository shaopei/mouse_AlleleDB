cd 3_CNV.Peak.BlackList.Filtered_counts_bed
folder=../6_CNV.Peak.BlackList.Filtered_counts_bed_in_AlleleHMMClusterSwitch
mkdir ${folder}
grep No ../GM12878_AlleleHMM_cluster_switch.bed > ${folder}/Concordant.bed #NoSwitch.bed
grep Yes ../GM12878_AlleleHMM_cluster_switch.bed > ${folder}/Discordant.bed #YesSwitch.bed
for f in *Filtered_counts.bed ;
do echo -e chrm$'\t'chromStart$'\t'snppos$'\t'${f:0:11}.winning$'\t'${f:0:11}.SymCls$'\t'${f:0:11}.SymPval > ${folder}/${f:0:-20}.GM.ConcordantRegions.Filtered_counts.bed;
echo -e chrm$'\t'chromStart$'\t'snppos$'\t'${f:0:11}.winning$'\t'${f:0:11}.SymCls$'\t'${f:0:11}.SymPval > ${folder}/${f:0:-20}.GM.DiscordantRegions.Filtered_counts.bed;
intersectBed -wa -wb -a $f -b ${folder}/Concordant.bed >> ${folder}/${f:0:-20}.GM.ConcordantRegions.Filtered_counts.bed;
intersectBed -wa -wb -a $f -b ${folder}/Discordant.bed >> ${folder}/${f:0:-20}.GM.DiscordantRegions.Filtered_counts.bed;
done

cd ${folder}
#python ../FisherExactTest_FromGM.DiscordantRegions.Filtered_counts.bed_AsymDomminat.py
python /workdir/sc2457/GM_GroSeq_AlleleDB_20170606/test_in_total/After_AlleleDB_pipeline/FisherExactTest_FromGM.DiscordantRegions.Filtered_counts.bed_AsymDomminat.py

