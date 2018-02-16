setwd("~/Box Sync/Danko_lab_work/mouse_AlleleDB/GM12878/Pairs_of_TRE/")
f1p="counts_both_hmm_regions_t1e-05_interestingHets_IGV_clusterBy20000bp.bed"

setwd("~/Box Sync/Danko_lab_work/mouse_AlleleDB/GM12878/Pairs_of_TRE/B6_CAST/")
f1p="counts_both_hmm_regions_t1e-05_interestingHets_IGV_clusterBy10000bp.bed"
f1=read.table(f1p,header=F)


hist(unlist(lapply(1:max(f1$V7), function(i){sum(f1$V7==i)})), breaks=seq(1,200,1))

dummy=1
f1$clusterSize=unlist(lapply(f1$V7, function(i){sum(f1$V7==i)}))
f1$winP=unlist(lapply(f1$V4, function(i){strsplit(as.character(i),',',fixed=TRUE)[[1]][1]}))
f1$matCount=unlist(lapply(f1$V4, function(i){as.numeric(strsplit(as.character(i),',',fixed=TRUE)[[1]][2])}))
f1$patCount=unlist(lapply(f1$V4, function(i){as.numeric(strsplit(as.character(i),',',fixed=TRUE)[[1]][3])}))
f1$mat_pat_log2ratio=log2((f1$matCount+dummy)/(f1$patCount+dummy))

cluster2=f1[f1$clusterSize==2,]  
View(cluster2)

pair1=unlist(lapply(seq(1,dim(cluster2)[1],2), function(i){cluster2$mat_pat_log2ratio[i]}))
pair2=unlist(lapply(seq(2,dim(cluster2)[1],2), function(i){cluster2$mat_pat_log2ratio[i]}))
plot(pair1,pair2)
abline(h=0)
abline(v=0)

pair1=unlist(lapply(seq(1,dim(cluster2)[1],2), function(i){cluster2$winP[i]}))
pair2=unlist(lapply(seq(2,dim(cluster2)[1],2), function(i){cluster2$winP[i]}))
sum(pair1==pair2)
sum(pair1!=pair2)
length(pair1)

