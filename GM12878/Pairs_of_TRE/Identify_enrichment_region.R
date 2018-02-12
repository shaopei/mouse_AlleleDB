setwd("~/Box Sync/Danko_lab_work/mouse_AlleleDB/GM12878/Pairs_of_TRE")
f1p="counts_both_hmm_regions_t1e-05_interestingHets_IGV.bed_distanceToNearestRegion.txt"
#f1p="counts_both_hmm_regions_t1e-05_interestingHets_IGV.bed_pairwise_distance.txt"
f1=read.table(f1p,header=F)

hist(f1$V1, col="blue",
     breaks = seq(0,5000000,25000), freq = FALSE,
     xlab="distance(bp)",
     main="distanceToNearestRegion") 

Uplimit=100000
hist(f1$V1[f1$V1<=Uplimit], col="blue", xlim=c(0,Uplimit),
     breaks = seq(0,Uplimit,2500), freq = FALSE,
     xlab="distance(bp)",
     main="distanceToNearestRegion")

# use MLE 
lamda = with(f1, 1/mean(V1))
lines(x, lamda*exp(-lamda*x), col="orange",lwd=1.5)
legend("topright", 
       legend = paste("Exponential pdf lamda = ", format(lamda, digits=2, nsmall=2),sep=""),
       lty=1, lwd=1.5, col="orange", bty = "n")





### use fitdistr
dis=f1$V1
fit1 <- fitdistr(dis, "geometric") 
fit2 <- fitdistr(dis, "exponential") 
hist(dis, freq = FALSE, breaks = 5000, xlim = c(0, 200000))
curve(dexp(x, rate = fit1$estimate), col = "red", add = TRUE)
curve(dexp(x, rate = fit2$estimate), col = "blue", add = TRUE)
curve(dexp(x, rate = lamda), col = "orange", add = TRUE)

## use lm 
Uplimit=20000
breaks=5000
dummy=1
a=with (f1, hist(V1, col="blue", 
                 xlab='distance (base)',
                 breaks = breaks))

count=a$counts+dummy
dis=a$mids
exponential.model = lm(log(count) ~ dis)
summary(exponential.model)
x=seq(0,max(dis),max(dis)/5000)
y=exp(coef(exponential.model)[1]+x*coef(exponential.model)[2])
lines(x, y, col="red")