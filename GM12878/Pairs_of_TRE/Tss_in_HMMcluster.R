setwd("~/Box Sync/Danko_lab_work/mouse_AlleleDB/GM12878/Pairs_of_TRE/")
setwd("~/Box Sync/Danko_lab_work/mouse_AlleleDB/GM12878/Pairs_of_TRE/B6_CAST/")
f2p="counts_both_hmm_regions_t1e-05_interestingHets_IGV_clusterBy20000bp.bed"
f2p="counts_both_hmm_regions_t1e-05_interestingHets_IGV_clusterBy10000bp.bed"
f2=read.table(f2p,header=F)
f2$clusterSize=unlist(lapply(f2$V7, function(i){sum(f2$V7==i)}))
View(f2)
cluster_size=cbind.data.frame(f2$V7, f2$clusterSize)
cluster_size=cluster_size[!duplicated(cluster_size),]
colnames(cluster_size)= c( "ClusterID" ,"ClusterSize")
View(cluster_size)
#cluster_size=as.numeric(cluster_size[2])
hist(cluster_size$ClusterSize, breaks=seq(-0.99,200.99,1), col="blue")

###
# f1p is TSS that overlap with Allele specific HMM region.
f1p="counts_both_hmm_regions_t1e-05_interestingHets_IGV_clusterBy20000bp_tss_all.bed"
f1p="counts_both_hmm_regions_t1e-05_interestingHets_IGV_clusterBy10000bp_dReg.bed"
# include cluster regions with only one HMM block
f1=read.table(f1p,header=F)
View(f1)
f1$winP=unlist(lapply(f1$V4, function(i){strsplit(as.character(i),',',fixed=TRUE)[[1]][1]}))
#f1$matCount=unlist(lapply(f1$V4, function(i){as.numeric(strsplit(as.character(i),',',fixed=TRUE)[[1]][2])}))
#f1$patCount=unlist(lapply(f1$V4, function(i){as.numeric(strsplit(as.character(i),',',fixed=TRUE)[[1]][3])}))

# number of M or P TSS in cluster, M|P TSS is determined by the HMM block that ovewrlape with the TSS
for (i in 1:max(f1$V7)){
  f1$M_TSS_count_cluster[f1$V7==i]<- sum(f1$winP[f1$V7==i]=="M")
  f1$P_TSS_count_cluster[f1$V7==i]<- sum(f1$winP[f1$V7==i]=="P")
}

# the start and end of each cluster, define by the TSS at the two end
for (i in 1:max(f1$V7)){
  f1$minStart[f1$V7==i]<- min(f1$V2[f1$V7==i])
  f1$maxEnd[f1$V7==i]<- max(f1$V3[f1$V7==i])
}

# the order of the TSS sequence
for (i in 1:max(f1$V7)){
  f1$tssSeq[f1$V7==i]<- paste(f1$winP[f1$V7==i],collapse ="")
}

for (i in 1:max(f1$V7)){
  f1$tssSeqStrand[f1$V7==i]<- paste(f1$winP[f1$V7==i], f1$V6[f1$V7==i],collapse ="", sep="")
}


WinP_count_cluster=cbind.data.frame(f1$V1,f1$minStart, f1$maxEnd, f1$V7,f1$tssSeq,f1$tssSeqStrand , f1$M_TSS_count_cluster, f1$P_TSS_count_cluster)

colnames(WinP_count_cluster) = c("Chr","ClusterStart", "ClusterEnd",  "ClusterID","Tss_order","Tss_order_strand","M_TSS_count_cluster", "P_TSS_count_cluster")
WinP_count_cluster=WinP_count_cluster[!duplicated(WinP_count_cluster),]
WinP_count_cluster$totalTSSCount=WinP_count_cluster$M_TSS_count_cluster + WinP_count_cluster$P_TSS_count_cluster
View(WinP_count_cluster)
dim(WinP_count_cluster)
WinP_count_cluster$ClusterSize=unlist(lapply(WinP_count_cluster$ClusterID, 
                function(i){cluster_size$ClusterSize[i]}))
  
  
hist(WinP_count_cluster$totalTSSCount, breaks = seq(0,max(WinP_count_cluster$totalTSSCount),1), add=TRUE, col="green")
WinP_count_cluster$M_TssRatio=WinP_count_cluster$M_TSS_count_cluster/WinP_count_cluster$totalTSSCount
hist(WinP_count_cluster$M_TssRatio)
hist(WinP_count_cluster$M_TssRatio[WinP_count_cluster$M_TssRatio>0 & WinP_count_cluster$M_TssRatio<1])
length(WinP_count_cluster$M_TssRatio[WinP_count_cluster$M_TssRatio>0 & WinP_count_cluster$M_TssRatio<1])

# number of switch between M and P
WinP_count_cluster$switch=unlist(lapply(WinP_count_cluster$Tss_order, 
              function(i){
                t=as.character(i)
                s=0
                c=substr(t,1,1)
                if (nchar(t)>1){
                  for (j in 2:nchar(t)){
                    if (substr(t,j,j) != c){s=s+1}
                    c=substr(t,j,j)
                  }
                }
                return(s)
                }
              ))

#WinP_count_cluster=WinP_count_cluster[,c(1,2,3,4,10,5,6,12,7:ncol(WinP_count_cluster))]
hist(WinP_count_cluster$switch#[WinP_count_cluster$switch>0]
     #,breaks = 9
     , breaks = seq(-0.5,max(WinP_count_cluster$switch)+0.5,1)
     ,col="blue"
     ,xlab="number of switches"
     ,main=""
     ,xlim=c(0,10)
     )
sum(WinP_count_cluster$switch==0)
sum(WinP_count_cluster$switch==1)
hist(WinP_count_cluster$switch[WinP_count_cluster$switch>0]
     , breaks = seq(-0.5,max(WinP_count_cluster$switch)+0.5,1))

# keep cluster with more than one HMM blocks
switch0=WinP_count_cluster[WinP_count_cluster$switch==0 & WinP_count_cluster$ClusterSize>1, ]
switch1=WinP_count_cluster[WinP_count_cluster$switch>0 & WinP_count_cluster$ClusterSize>1, ]
switch0$switch="No"
switch1$switch="Yes"
s = rbind.data.frame(switch0[,c(1,2,3,12)], switch1[,c(1,2,3,12)])
write.table(s, "GM12878_AlleleHMM_cluster_switch.bed"
            , sep = "\t",row.names = F,quote = F,
            col.names = F)

hist(WinP_count_cluster$totalTSSCount[WinP_count_cluster$ClusterSize>1])
hist(WinP_count_cluster$totalTSSCount[WinP_count_cluster$ClusterSize>1 & WinP_count_cluster$switch==0], add=T, col=mycol1)
hist(WinP_count_cluster$totalTSSCount[WinP_count_cluster$ClusterSize>1 & WinP_count_cluster$switch>0], add=T, col=mycol2)

## some filter 
hist(WinP_count_cluster$M_TssRatio[WinP_count_cluster$ClusterSize>1])
length(WinP_count_cluster$M_TssRatio[WinP_count_cluster$ClusterSize>1])
length(WinP_count_cluster$M_TssRatio[WinP_count_cluster$ClusterSize>1 & WinP_count_cluster$M_TssRatio>0 & WinP_count_cluster$M_TssRatio<1])

WinP_count_cluster = WinP_count_cluster[WinP_count_cluster$switch>0,]
WinP_count_cluster$ClusterSize_kb=WinP_count_cluster$ClusterEnd- WinP_count_cluster$ClusterStart
hist(WinP_count_cluster$ClusterSize_kb)
#WinP_count_cluster=WinP_count_cluster[WinP_count_cluster$M_TssRatio>0 & WinP_count_cluster$M_TssRatio<1,]

#WinP_count_cluster=WinP_count_cluster[WinP_count_cluster$totalTSSCount>5,]
#WinP_count_cluster=WinP_count_cluster[WinP_count_cluster$M_TSS_count_cluster>0,]

