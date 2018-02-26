setwd("~/Box Sync/Danko_lab_work/mouse_AlleleDB/GM12878/Pairs_of_TRE/")
f1p="counts_both_hmm_regions_t1e-05_interestingHets_IGV_clusterBy20000bp.bed"

setwd("~/Box Sync/Danko_lab_work/mouse_AlleleDB/GM12878/Pairs_of_TRE/B6_CAST/")
f1p="counts_both_hmm_regions_t1e-05_interestingHets_IGV_clusterBy10000bp.bed"
f1=read.table(f1p,header=F)


hist(unlist(lapply(1:max(f1$V7), function(i){sum(f1$V7==i)})), breaks=seq(1,200,1))

dummy=1
f1$clusterSize=unlist(lapply(f1$V7, function(i){sum(f1$V7==i)}))
f1$winP=unlist(lapply(f1$V4, function(i){strsplit(as.character(i),',',fixed=TRUE)[[1]][1]}))
# matCount is matReadsCount
f1$matReadsCount=unlist(lapply(f1$V4, function(i){as.numeric(strsplit(as.character(i),',',fixed=TRUE)[[1]][2])}))
f1$patReadsCount=unlist(lapply(f1$V4, function(i){as.numeric(strsplit(as.character(i),',',fixed=TRUE)[[1]][3])}))
f1$mat_pat_ReadsCount_log2ratio=log2((f1$matReadsCount+dummy)/(f1$patReadsCount+dummy))

# the start and end of each cluster
for (i in 1:max(f1$V7)){
  f1$minStart[f1$V7==i]<- min(f1$V2[f1$V7==i])
  f1$maxEnd[f1$V7==i]<- max(f1$V3[f1$V7==i])
}

# keep cluster with more than one HMM blocks
cluster_moreThan1=f1[f1$clusterSize>1,c(1,13,14,8)] 
cluster_moreThan1 = cluster_moreThan1[!duplicated(cluster_moreThan1),]
View(cluster_moreThan1)
write.table(cluster_moreThan1, paste(f1p,"_cluster_moreThan1.bed",sep="")
            , sep = "\t",row.names = F,quote = F,
            col.names = F)

# study correlations below
cluster2=f1[f1$clusterSize==2,]  
View(cluster2)

pair1_s=unlist(lapply(seq(1,dim(cluster2)[1],2), function(i){cluster2$V6[i]}))
pair2_s=unlist(lapply(seq(2,dim(cluster2)[1],2), function(i){cluster2$V6[i]}))
sum(pair1_s==pair2_s)
sum(pair1_s!=pair2_s)
pair1=unlist(lapply(seq(1,dim(cluster2)[1],2), function(i){cluster2$mat_pat_ReadsCount_log2ratio[i]}))
pair2=unlist(lapply(seq(2,dim(cluster2)[1],2), function(i){cluster2$mat_pat_ReadsCount_log2ratio[i]}))
plot(pair1[pair1_s!=pair2_s],pair2[pair1_s!=pair2_s]
     , xlab="mat_pat_log2ratio"
     , ylab="mat_pat_log2ratio"
     , col="orange")
points(pair1[pair1_s==pair2_s],pair2[pair1_s==pair2_s])
abline(h=0)
abline(v=0)

pair1=unlist(lapply(seq(1,dim(cluster2)[1],2), function(i){cluster2$winP[i]}))
pair2=unlist(lapply(seq(2,dim(cluster2)[1],2), function(i){cluster2$winP[i]}))
sum(pair1==pair2)
sum(pair1!=pair2)
length(pair1)

###
cluster3=f1[f1$clusterSize==3,]  
View(cluster3)
tri1s=unlist(lapply(seq(1,dim(cluster3)[1],3), function(i){cluster3$V6[i]}))
tri2s=unlist(lapply(seq(2,dim(cluster3)[1],3), function(i){cluster3$V6[i]}))
tri3s=unlist(lapply(seq(3,dim(cluster3)[1],3), function(i){cluster3$V6[i]}))

tri1=unlist(lapply(seq(1,dim(cluster3)[1],3), function(i){cluster3$mat_pat_ReadsCount_log2ratio[i]}))
tri2=unlist(lapply(seq(2,dim(cluster3)[1],3), function(i){cluster3$mat_pat_ReadsCount_log2ratio[i]}))
tri3=unlist(lapply(seq(3,dim(cluster3)[1],3), function(i){cluster3$mat_pat_ReadsCount_log2ratio[i]}))
plot(tri1[tri1s!=tri2s],tri2[tri1s!=tri2s]
     , xlab="mat_pat_log2ratio"
     , ylab="mat_pat_log2ratio"
     , col="orange")
points(tri1[tri1s!=tri3s],tri3[tri1s!=tri3s]
       , col="orange")
points(tri2[tri2s!=tri3s],tri3[tri2s!=tri3s]
       , col="orange")
points(tri1[tri1s==tri2s],tri2[tri1s==tri2s])
points(tri1[tri1s==tri3s],tri3[tri1s==tri3s])
points(tri2[tri2s==tri3s],tri3[tri2s==tri3s])
length(c(tri1[tri1s!=tri2s], tri2[tri2s!=tri3s], tri1[tri1s!=tri3s]))
length(c(tri1[tri1s==tri2s], tri2[tri2s==tri3s], tri1[tri1s==tri3s]))

sum(tri1[tri1s==tri2s]*tri2[tri1s==tri2s] <0)+
sum(tri1[tri1s==tri3s]*tri3[tri1s==tri3s] <0)+
sum(tri2[tri2s==tri3s]*tri3[tri2s==tri3s] <0)
