setwd("~/Box Sync/Danko_lab_work/mouse_AlleleDB/GM12878/Pairs_of_TRE/")
setwd("~/Box Sync/Danko_lab_work/mouse_AlleleDB/GM12878/Pairs_of_TRE/B6_CAST/")

f2p="20000_cluster_region_count.txt"
f2p="counts_both_hmm_regions_t1e-05_interestingHets_IGV_clusterBy10000bp_region_count.txt"
f2=read.table(f2p,header=F)
hist(f2$V1,col=mycol1
     , breaks = seq(0,200,1)
     , freq = F
     #, xlab="HMM blocks in the cluster"
     #, add=T
     , xlim=c(0,50)
     ,xlab="Number of HMM blocks or TSSs in each cluster"
     ,main=""
)
lamda = with(f2, 1/mean(V1))
x=seq(0,200,1)
lines(x, lamda*exp(-lamda*x), col="dark orange",lwd=1.5)
#lines(density(f2$V1), col="light blue") 
legend("topright", 
       legend = paste("Exponential pdf lamda = ", format(lamda, digits=2, nsmall=2),sep=""),
       lty=c(0,1,0,1), lwd=1.5, col="orange", bty = "n")

f3p="counts_both_hmm_regions_t1e-05_interestingHets_IGV_clusterBy20000bp_tssPerClusterCount.txt" 
f3p="counts_both_hmm_regions_t1e-05_interestingHets_IGV_clusterBy10000bp_dRegPerClusterCount.txt"
f3=read.table(f3p,header=F)
#hist(f3$V1, freq = F, col=mycol2)
hist(f3$V1, freq = F
     , breaks = seq(0,200,1)
     , add=T
     , xlab="Numer of TSSs in the cluster"
     , col=mycol2)
#lines(density(f3$V1), col="green") 
lamda2 = with(f3, 1/mean(V1))
lines(x, lamda2*exp(-lamda2*x), col="red",lwd=1.5)
legend("topright", 
       legend = paste("Exponential pdf lamda = ", format(lamda2, digits=2, nsmall=2),sep=""),
       lty=1, lwd=1.5, col="red", bty = "n")

legend("topright", 
       legend = c("AlleleHMM",paste("AlleleHMM Exponential pdf lamda = ", format(lamda, digits=2, nsmall=2),sep="")
                  ,"TSS",paste("TSS Exponential pdf lamda = ", format(lamda2, digits=2, nsmall=2),sep="")),
       pch=c(15,NA,15,NA),
       cex=0.8, 
       lty=c(0,1,0,1),
       #bty="n",
       lwd=1.5, 
       col=c(mycol1,"dark orange",mycol2, "red"), bty = "n")


###color transparant
t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #       percent = % transparency
  #          name = an optional name for the color
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100,
               names = name)
  ## Save the color
  invisible(t.col)
}

mycol1 <- t_col("blue", perc = 50, name = "lt.blue")
mycol2 <- t_col("green", perc = 50, name = "lt.green")
## END