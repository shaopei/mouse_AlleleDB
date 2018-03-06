cd 2_CNV.PeakFiltered_counts_bed/
mkdir ../3_CNV.Peak.BlackList.Filtered_counts_bed
for f in *Filtered_counts.bed ;
do echo -e \#chrm$'\t'chromStart$'\t'snppos$'\t'${f:0:11}.winning$'\t'${f:0:11}.SymCls$'\t'${f:0:11}.SymPval > ../3_CNV.Peak.BlackList.Filtered_counts_bed/${f:0:-19}.BlackList.Filtered_counts.bed;
intersectBed -wa -v -a $f -b /workdir/sc2457/GM_GroSeq_AlleleDB_20170606/ENCODE_20170609/ENCFF001TDO.bed >> ../3_CNV.Peak.BlackList.Filtered_counts_bed/${f:0:-19}.BlackList.Filtered_counts.bed;
done
