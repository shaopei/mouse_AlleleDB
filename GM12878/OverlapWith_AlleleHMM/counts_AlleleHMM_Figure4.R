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



setwd("~/Box Sync/Danko_lab_work/mouse_AlleleDB/GM12878/OverlapWith_AlleleHMM")
A=read.table("interestingHets_AlleleDB_in_AlleleHMM_MP_switch_counts.txt")
D=read.table("interestingHets_AlleleDB_out_AlleleHMM_MP_switch_counts.txt")
H=read.table("counts_noX_MinCount1_inAlleleHMM_t1e-05_interestingHets_outAlleleDB_switch_counts.txt")
xmax = max(c(A$V1, D$V1, H$V1))
par(mfrow=c(3,1))
par(cex.lab=2.2, cex.axis=2.2)
u=10
H$V1[H$V1>=u]=u
A$V1[A$V1>=u]=u
D$V1[D$V1>=u]=u
hist(H$V1-1, breaks = seq(-0.5,xmax,1), freq=F,col="red", xlim=c(0,10), ylim=c(0,1), main=NA, xlab="number of switches")
hist(A$V1-1, breaks = seq(-0.5,xmax,1), freq=F,col="purple", xlim=c(0,10), ylim=c(0,1), main=NA, xlab="number of switches") #,density=20,angle=180
hist(D$V1-1, breaks = seq(-0.5,xmax,1), freq=F,col="blue", xlim=c(0,10), ylim=c(0,1), main=NA, xlab="number of switches") #,density=20,angle=45






##
# Simple Pie Chart
lbls <- c("Concordant", "Discordant", "Symmetric")
H <- c( 7957, 399, 2250)
A <- c( 5274, 512, 1993)
D <- c( 15926, 7681, 34356)
#par(mfrow=c(1,3))
pie(H, labels = lbls, main="H")
pie(A, labels = lbls, main="A")
pie(D, labels = lbls, main="D")
legend(labels=lbls)

legend("topright", 
       title = " ",
       legend =lbls,
       #cex=1, 
       #lty=c(0,0),
       #lwd=1.5, 
       #fill=c(mycol1,mycol2, "red"), 
       bty = "n")

## read counts of SNPs in H (AlleleHMM, not AlleleDB), A (AlleleHMM and AlleleDB), and D (AlleleDB, not AlleleHMM)
H_Con=read.table("H_Concordant_counts.txt")
H_Dis=read.table("H_Discordant_counts.txt")
H_Sym=read.table("H_Symmetric_counts.txt")
A_Con=read.table("A_Concordant_counts.txt")
A_Dis=read.table("A_Discordant_counts.txt")
A_Sym=read.table("A_Symmetric_counts.txt")
D_Con=read.table("D_Concordant_counts.txt")
D_Dis=read.table("D_Discordant_counts.txt")
D_Sym=read.table("D_Symmetric_counts.txt")
H=c(H_Con$V1, H_Dis$V1, H_Sym$V1)
A=c(A_Con$V1, A_Dis$V1, A_Sym$V1)
D=c(D_Con$V1, D_Dis$V1, D_Sym$V1)

u=100
H[H>=u]=u
A[A>=u]=u
D[D>=u]=u

par(mfrow=c(3,1))
boxplot(D)
boxplot(A)
boxplot(H)

b=5
hist(H, freq=F, col="red", breaks =  seq(1,u+b,b), xlab="Read counts per SNP", ylab="fraction of SNPs", ylim=c(0,0.2))
hist(A, freq=F, col="purple", breaks =  seq(1,u+b,b), xlab="Read counts per SNP", ylab="fraction of SNPs", ylim=c(0,0.2))
hist(D, freq=F, col="blue", breaks =  seq(1,u+b,b), xlab="Read counts per SNP", ylab="fraction of SNPs", ylim=c(0,0.2))
#,density=20,angle=135

hist(H_Con$V1,freq=F, col="green",density=20,angle=135)
hist(H_Dis$V1,freq=F, add=T, col="orange",density=20,angle=135)
hist(H_Sym$V1,freq=F, add=T, col="blue",density=20,angle=135)

hist(A_Sym$V1,freq=F, col="blue",density=20,angle=135, breaks =  seq(0,15000,100), xlim=c(0,500))
hist(A_Con$V1,freq=F, add=T, col="green",density=20,angle=135, breaks =  seq(0,15000,100))
hist(A_Dis$V1,freq=F, add=T, col="orange",density=20,angle=135, breaks =  seq(0,15000,100))

hist(D_Sym$V1,freq=F, col="blue",density=20,angle=135, breaks =  seq(0,51000,100), xlim=c(0,500))
hist(D_Con$V1,freq=F, add=T, col="green",density=20,angle=135, breaks =  seq(0,15000,100))
hist(D_Dis$V1,freq=F, add=T, col="orange",density=20,angle=135, breaks =  seq(0,15000,100))

## read counts of SNPs in H (AlleleHMM, not AlleleDB), A (AlleleHMM and AlleleDB), and D (AlleleDB, not AlleleHMM)
## end


hist(c(rep(0,EBlock_with_no_gene), Em_withG$V1, Ep_withG$V1)
     , breaks = seq(-0.5,200,1)
     #,density=20,angle=135
     , freq = F ,ylim=c(0,0.6)
     # ,ylim=c(0,2500)
     , xlim=c(0,15)
     ,xlab="Number of genes in each HMM block"
     ,ylab="AlleleHMM block count"
     ,col= e_col
     ,main=""
     
     ,las=1
     ,bty="l"
)
hist(c(rep(0,GM_Block_with_no_gene), GM_m_wGene$V1, GM_m_wGene$V1)
     , breaks = seq(-0.5,200,1)
     , xlim=c(0,20)
     , freq = F
     ,density=25,angle=45
     ,col= g_col
     ,add=T
)