cd 0_counts_files
mkdir ../1_CNVFiltered_counts_bed
for f in *_counts.txt ;
do echo -e \#chrm$'\t'chromStart$'\t'snppos$'\t'${f:0:-11}.winning$'\t'${f:0:-11}.SymCls$'\t'${f:0:-11}.SymPval > ../1_CNVFiltered_counts_bed/${f:0:-11}.CNVFiltered_counts.bed;
cat $f | awk 'BEGIN{OFS="\t"} (NR >1 && $18 >0.5 && $18<1.5) {print "chr"$1,$2-1, $2,$14,$15,$16}' >> ../1_CNVFiltered_counts_bed/${f:0:-11}.CNVFiltered_counts.bed;
done
