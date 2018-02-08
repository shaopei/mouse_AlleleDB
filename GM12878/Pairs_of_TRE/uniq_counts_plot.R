#R --vanilla --slave --args $(pwd) "counts_hmm_regions_t*_interestingHets_5head_distance_toclosest-dReg_AccumulateCounts.txt" counts_hmm_regions_tX_interestingHets_5head_distance_toclosest-dReg_AccumulateCounts.pdf counts_hmm_regions_tX_interestingHets_5head_distance_toclosest-dReg_At5Kb.pdf < getFractionOfBlock_DistanceToNearestSites.R

#arguments here
args=(commandArgs(TRUE))
setwd(args[1])
f_pattern = glob2rx(args[2])


setwd("~/Box Sync/Danko_lab_work/mouse_AlleleDB/GM12878/Pairs_of_TRE")
col1 = 'dark green'#'springgreen4'
col2 ='dark orange'
f1p="same_par_accumulate_counts.txt"
f1=read.table(f1p,header=F)
plot(f1$V1, f1$V2/max(f1$V2), 
     main="to the nearest TSS on the opposite strand", 
     ylab="accumulate proportion",
     xlab="Distance to the nearest TSS on the opposite strand",
     col=col1)

f2p="oppo_par_accumulate_counts.txt"
f2=read.table(f2p,header=F)
points(f2$V1, f2$V2/max(f2$V2), col=col2)

legend("bottomright", legend = c("same parent", "diff parent"), col = c(col1, col2),
       pch=1
       )

#f3p="nearest_oppo_parentTss_accumulate_counts.txt"
f3p="nearest_oppo_parent_AlleleHMM_accumulate_counts.txt"
f3=read.table(f3p,header=F)
#f4p="nearest_same_parentTss_accumulate_counts.txt"
f4p="nearest_same_parent_AlleleHMM_accumulate_counts.txt"
f4=read.table(f4p,header=F)
f5p="nearest_otherStrand_AlleleHMM_accumulate_counts.txt"
f5=read.table(f5p,header=F)
plot(f5$V1, f5$V2, col="black")
points(f4$V1, f4$V2, col="blue")
points(f3$V1, f3$V2, col="red")
points(f2$V1, f2$V2, col=col2)
points(f1$V1, f1$V2, col=col1)

legend("bottomright", legend = c("to the nearest AlleleHMM with same parent",
                                 "to the nearest AlleleHMM with diff parent",
                                 "to the nearest AlleleHMM",
                                 "same parent", "diff parent"), col = c("blue","red","black",col1, col2),
       pch=1
)

legend("bottomright", legend = c("to the nearest TSS with same parent",
                                 "to the nearest TSS with diff parent",
                                 "to the nearest TSS",
                                 "same parent", "diff parent"), col = c("blue","red","black",col1, col2),
       pch=1
)


f6p="nearest_bothStrandAlleleHMM_accumulate_counts.txt"
f6=read.table(f6p,header=F)
#names(f6)[1] = "distance"
#names(f6)[2] = "count"
f7p="nearest_bothStrandAlleleHMM_RawCount.txt"
f7=read.table(f7p,header=F)

plot(f6$V1[f6$V1<=50000], f6$V2[f6$V1<=50000], col="black",
    xlab= "Distance to the nearest TSS on either strand",
    ylab= "frequency")

hist(f7$V1[f7$V1 <= 50000], col="blue", xlab='distance', breaks = 100)
hist(f7$V1[f7$V1 <= 10000], col="blue", xlab='distance', breaks = 100)
#View(f7)

f8p="nearest_bothStrandAlleleHMM_count.txt"
f8=read.table(f8p,header=F)
names(f8)[2] = "distance"
names(f8)[1] = "count"
plot(f8$distance[f8$distance <=50000], f8$count[f8$distance <=50000])#, type="l")
