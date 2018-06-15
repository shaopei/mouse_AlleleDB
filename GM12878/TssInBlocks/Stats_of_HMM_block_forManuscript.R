#functions
hist_of_HMM_block_from_bed6 <- function(m_all, p_all, add=F, col="blue"){
  mp_all = rbind.data.frame(m_all, p_all)
  mp_all$blockSize = mp_all$V3 - mp_all$V2
  hist(log10(mp_all$blockSize),
       #freq=F,
       col = col,
       xlab="log10(HMM block size in bp)",
       ylab="fraction",
       main="The size distribution of HMM blocks",
       add=add)
    }


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
mycol3 <- t_col("red", perc = 50, name = "lt.red")




#GM12878
setwd("~/Box Sync/Danko_lab_work/mouse_AlleleDB/GM12878/TssInBlocks/")
m_all=read.table("counts_minus_hmm_regions_t1e-05_interestingHets_IGV.bed")
p_all=read.table("counts_plus_hmm_regions_t1e-05_interestingHets_IGV.bed")
hist_of_HMM_block_from_bed6(m_all, p_all, F, t_col("blue", perc = 50))
p_all$blockSize = p_all$V3 - p_all$V2
m_all$blockSize = m_all$V3 - m_all$V2
blockSize=c(m_all$blockSize,p_all$blockSize)

m_wTss=read.table("counts_minus_hmm_regions_t1e-05_interestingHets_IGV_tssallPerBlock.txt")
p_wTss=read.table("counts_plus_hmm_regions_t1e-05_interestingHets_IGV_tssallPerBlock.txt")
m_wTss=read.table("counts_minus_hmm_regions_t1e-05_interestingHets_IGV_dRegPeakPerBlock.txt")
p_wTss=read.table("counts_plus_hmm_regions_t1e-05_interestingHets_IGV_dRegPeakPerBlock.txt")

mp_wTss = rbind.data.frame(m_wTss, p_wTss)
mp_wTss$blockSize = mp_wTss$V4 - mp_wTss$V3
hist(log10(mp_wTss$blockSize),
     #breaks = seq(0,7,0.5),
     #freq=F, 
     col=t_col("green",  perc = 50),
     xlab="log10(HMM block size in bp)",
     ylab="fraction",
     main="The size distribution of HMM blocks"
     ,add=T
     )
legend("topleft", 
       title = "F1 SRR4041366_dedup_2",
       legend = c("all blocks", "blocks with TSS"),
       cex=1, 
       lty=c(0,0),
       lwd=1.5, 
       fill=c(mycol1,mycol2), bty = "n")

legend("topright", 
       title = "GM12878 read amount",
       legend = c("188M (all)", "47 M", "20M"),
       cex=1, 
       lty=c(0,0),
       lwd=1.5, 
       fill=c(mycol1,mycol2, "red"), bty = "n")

Em_withTSS=read.table("SRR4041366_dedup_2/counts_minus_hmm_regions_t1e-05_interestingHets_IGV_tssallPerBlock.txt")
Ep_withTSS=read.table("SRR4041366_dedup_2/counts_plus_hmm_regions_t1e-05_interestingHets_IGV_tssallPerBlock.txt")
EBlock_with_no_TSS= 3482- dim(Em_withTSS)[1]-dim(Ep_withTSS)[1]
Emp_wTss = rbind.data.frame(Em_withTSS, Ep_withTSS)
Emp_wTss$blockSize = Emp_wTss$V4 - Emp_wTss$V3
hist(log10(Emp_wTss$blockSize),
     #freq=F, 
     col=t_col("green",  perc = 50),
     xlab="log10(HMM block size in bp)",
     ylab="fraction",
     main="The size distribution of HMM blocks in Engreitz_F1"
     ,add=T
)

legend("topleft", 
       title = "Engreitz_F1",
       legend = c("all blocks", "blocks with TSS"),
       #pch=c(15,15),
       cex=1, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       #density=c(20,20), #angle=45,
       #angle=c(135,45),
       fill=c(mycol1,mycol2), bty = "n")



#Engreitz_F1
m_all_SRR4041366 = read.table("~/Box Sync/Danko_lab_work/mouse_AlleleDB/Engreitz_F1_20180523/allelicbias-PersonalGenome_P.CAST_M.129S1-SRR4041366_dedup_2/SRR4041366_dedup_2_counts_minus_hmm_regions_t1e-05_interestingHets_IGV.bed")
p_all_SRR4041366 = read.table("~/Box Sync/Danko_lab_work/mouse_AlleleDB/Engreitz_F1_20180523/allelicbias-PersonalGenome_P.CAST_M.129S1-SRR4041366_dedup_2/SRR4041366_dedup_2_counts_plus_hmm_regions_t1e-05_interestingHets_IGV.bed")
hist_of_HMM_block_from_bed6(m_all_SRR4041366, p_all_SRR4041366, add=T, col=mycol2)
hist_of_HMM_block_from_bed6(m_all_SRR4041366, p_all_SRR4041366, col=mycol2)
hist_of_HMM_block_from_bed6(m_all_SRR4041366, p_all_SRR4041366, col=mycol1)

m_all_SRR4041367 = read.table("~/Box Sync/Danko_lab_work/mouse_AlleleDB/Engreitz_F1_20180523/allelicbias-PersonalGenome_P.CAST_M.129S1-SRR4041367_dedup_2/SRR4041367_dedup_2_counts_minus_hmm_regions_t1e-05_interestingHets_IGV.bed")
p_all_SRR4041367 = read.table("~/Box Sync/Danko_lab_work/mouse_AlleleDB/Engreitz_F1_20180523/allelicbias-PersonalGenome_P.CAST_M.129S1-SRR4041367_dedup_2/SRR4041367_dedup_2_counts_plus_hmm_regions_t1e-05_interestingHets_IGV.bed")
hist_of_HMM_block_from_bed6(m_all_SRR4041367, p_all_SRR4041367, add=T, col=mycol1)

m_all_SRR4041367_2R = read.table("~/Box Sync/Danko_lab_work/mouse_AlleleDB/Engreitz_F1_20180523/allelicbias-PersonalGenome_P.CAST_M.129S1-SRR4041367_dedup_2_and_RC1/SRR4041367_dedup_2_and_RC1_counts_minus_hmm_regions_t1e-05_interestingHets_IGV.bed")
p_all_SRR4041367_2R = read.table("~/Box Sync/Danko_lab_work/mouse_AlleleDB/Engreitz_F1_20180523/allelicbias-PersonalGenome_P.CAST_M.129S1-SRR4041367_dedup_2_and_RC1/SRR4041367_dedup_2_and_RC1_counts_plus_hmm_regions_t1e-05_interestingHets_IGV.bed")
hist_of_HMM_block_from_bed6(m_all_SRR4041367, p_all_SRR4041367_2R, add=T, col=mycol3)

# LEP_ZYG
m_all_LZ = read.table("~/Box Sync/Danko_lab_work/mouse_AlleleDB/Engreitz_F1_20180523/allelicbias-PersonalGenome_P.CAST_M.B6-LEP_ZYG_ATGCA_forAlleleDB/LEP_ZYG_counts_minus_hmm_regions_t1e-05_interestingHets_IGV.bed")
p_all_LZ = read.table("~/Box Sync/Danko_lab_work/mouse_AlleleDB/Engreitz_F1_20180523/allelicbias-PersonalGenome_P.CAST_M.B6-LEP_ZYG_ATGCA_forAlleleDB/LEP_ZYG_counts_plus_hmm_regions_t1e-05_interestingHets_IGV.bed")
hist_of_HMM_block_from_bed6(m_all_LZ, p_all_LZ, T, "orange")


## gene per block
#GM12878
setwd("~/Box Sync/Danko_lab_work/mouse_AlleleDB/GM12878/GeneInBlocks/")
m_withG=read.table("counts_minus_hmm_regions_t1e-05_interestingHets_IGV_genePerBlock.txt")
p_withG=read.table("counts_plus_hmm_regions_t1e-05_interestingHets_IGV_genePerBlock.txt")
Block_with_no_gene= 4171 -  3123   

m_withG_atleastHalf = read.table("counts_minus_hmm_regions_t1e-05_interestingHets_IGV_geneMoreThanHalfPerBlock.txt")
p_withG_atleastHalf = read.table("counts_plus_hmm_regions_t1e-05_interestingHets_IGV_geneMoreThanHalfPerBlock.txt")
Block_with_no_gene_atleastHalf = 4171 -885

m_withG_whole = read.table("counts_minus_hmm_regions_t1e-05_interestingHets_IGV_wholeGenePerBlock.txt")
p_withG_whole = read.table("counts_plus_hmm_regions_t1e-05_interestingHets_IGV_wholeGenePerBlock.txt")
Block_with_no_gene_whole =4171- dim(m_withG_whole)[1]-dim(p_withG_whole)[1]


hist(c(rep(0,Block_with_no_gene), m_withG$V1, p_withG$V1)
     ,col= g_col
     #,density=20,angle=45
     , breaks = seq(-0.0001,200,1)
     , freq = F
     , xlim=c(0,20)
     ,xlab="Number of genes in each HMM block"
     ,main=" "
     , add=T
)

hist(c(rep(0,Block_with_no_gene_atleastHalf), m_withG_atleastHalf$V1, p_withG_atleastHalf$V1)
     ,col=t_col("green", perc = 50)
     #,density=20,angle=180
     , breaks = seq(-0.0001,200,1)
     , freq = F
     , xlim=c(0,20)
     ,xlab="Number of genes in each HMM block"
     ,main="GM12878"
     , add=T
)
hist(c(rep(0,Block_with_no_gene_whole), m_withG_whole$V1, p_withG_whole$V1)
     ,col="black" #g_col
     #,density=20,angle=180
     , breaks = seq(-0.0001,200,1)
     , freq = F
     , xlim=c(0,20)
     ,xlab="Number of genes in each HMM block"
     ,main="GM12878"
     #, add=T
)

legend("topright", 
       legend = c("F1 hybrid 129/castaneus", "GM12878"),
       #pch=c(15,15),
       cex=1, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(20,20), #angle=45,
       angle=c(135,45),
       fill=c("red","blue"), bty = "n")
##SRR4041366_dedup_2
Em_withG=read.table("SRR4041366_dedup_2/counts_minus_hmm_regions_t1e-05_interestingHets_IGV_genePerBlock.txt")
Ep_withG=read.table("SRR4041366_dedup_2/counts_plus_hmm_regions_t1e-05_interestingHets_IGV_genePerBlock.txt")
EBlock_with_no_gene= 3482-2522
Em_withG_atleastHalf = read.table("SRR4041366_dedup_2/counts_minus_hmm_regions_t1e-05_interestingHets_IGV_geneMoreThanHalfPerBlock.txt")
Ep_withG_atleastHalf = read.table("SRR4041366_dedup_2/counts_plus_hmm_regions_t1e-05_interestingHets_IGV_geneMoreThanHalfPerBlock.txt")
EBlock_with_no_gene_atleastHalf =3482- dim(Em_withG_atleastHalf)[1]-dim(Ep_withG_atleastHalf)[1]
Em_withG_whole = read.table("SRR4041366_dedup_2/counts_minus_hmm_regions_t1e-05_interestingHets_IGV_wholeGenePerBlock.txt")
Ep_withG_whole = read.table("SRR4041366_dedup_2/counts_plus_hmm_regions_t1e-05_interestingHets_IGV_wholeGenePerBlock.txt")
EBlock_with_no_gene_whole =3482- dim(Em_withG_whole)[1]-dim(Ep_withG_whole)[1]



hist(c(rep(0,EBlock_with_no_gene), Em_withG$V1, Ep_withG$V1)
     , breaks = seq(-0.0001,200,1)
     #,density=20,angle=135
     , freq = F
     , xlim=c(0,20)
     #,ylim=c(0,0.6)
     ,xlab="Number of genes in each HMM block"
     ,main=""
     ,col= e_col
    # ,add =T
    
)
hist(c(rep(0,EBlock_with_no_gene_atleastHalf), Em_withG_atleastHalf$V1, Ep_withG_atleastHalf$V1)
     ,col=t_col("green", perc = 50)
     #,density=20,angle=180
     , breaks = seq(-0.0001,200,1)
     #, freq = F
     , xlim=c(0,20)
     ,xlab="Number of genes in each HMM block"
     ,main="F1 hybrid 129/castaneus"
     , add=T
)
hist(c(rep(0,EBlock_with_no_gene_whole), Em_withG_whole$V1, Ep_withG_whole$V1)
     ,col=t_col("black", perc = 50)
     #,density=20,angle=180
     , breaks = seq(-0.0001,200,1)
     #, freq = F
     , xlim=c(0,20)
     ,xlab="Number of genes in each HMM block"
     ,main="F1 hybrid 129/castaneus"
     #, add=T
)

legend("topright", 
       legend = c("100%", "at least 50%", "any"),
       #pch=c(15,15),
       cex=1, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       #density=c(20,20), #angle=45,
       #angle=c(135,45),
       fill=c("black","green","red"), bty = "n")

legend("topright", 
       legend = c("F1 hybrid 129/castaneus", "GM12878"),
       #pch=c(15,15),
       cex=1, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       #density=c(20,20), #angle=45,
       #angle=c(135,45),
       fill=c(g_col,e_col)
       #fill=c("red","blue")
       , bty = "n")

### try to fit a distribution, not finish yet ###
a <- hist(c(m_withG$V1, p_withG$V1),col=mycol1
          , breaks = seq(0,200,1)
          , freq = F
          , xlim=c(0,20)
          ,xlab="Number of genes in each HMM block"
          ,main=""
)
library(fitdistrplus)
f <-fitdist (c(m_withG$V1, p_withG$V1), "geom")#,method = c("mle", "mme", "qme", "mge"))
hist(rpois(1000, f$estimate),add=T, freq = F)

lamda2 = 1/mean(c(m_withG$V1, p_withG$V1))
x=seq(0,200,1)
lines(x, lamda2*exp(-lamda2*x), col="red",lwd=1.5)
### try to fit a distribution, not finish yet ###

######### figure for paper ########
# figure patameter
# export  7.92 x 5.92 inches
e_col= "red" #t_col("red", perc = 20)
g_col="blue" #t_col("blue", perc = 20)
par(mar=c(5.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=1.8, cex.axis=1.8)


# AlleleHMM Block size distribution
setwd("~/Box Sync/Danko_lab_work/mouse_AlleleDB/GM12878/TssInBlocks/")
hist_of_HMM_block_from_bed6 <- function(m_all, p_all, add=F, col="blue"){
  mp_all = rbind.data.frame(m_all, p_all)
  mp_all$blockSize = mp_all$V3 - mp_all$V2
  hist(log10(mp_all$blockSize),
       #freq=F,
       col = col,
       xlim=c(0,8),
       xlab="log10(AlleleHMM block size in bp)",
       #ylab="number of AlleleHMM blocks",
       main="",
       add=add,
       #cex.lab=1.8, cex.axis=1.8,
       las=1,
       bty="l")
}
hist_of_HMM_block_from_bed6(m_all_SRR4041366, p_all_SRR4041366, col= e_col)

hist_of_HMM_block_from_bed6_layer2 <- function(m_all, p_all, add=F, col="blue"){
  mp_all = rbind.data.frame(m_all, p_all)
  mp_all$blockSize = mp_all$V3 - mp_all$V2
  hist(log10(mp_all$blockSize),
       #freq=F,
       col = col,
       density=30,angle=45,
       xlab="log10(HMM block size in bp)",
       ylab="fraction",
       main="The size distribution of HMM blocks",
       add=add)
}

#Engreitz_F1
m_all_SRR4041366 = read.table("SRR4041366_dedup_2/counts_minus_hmm_regions_t1e-05_interestingHets_IGV.bed")
p_all_SRR4041366 = read.table("SRR4041366_dedup_2/counts_plus_hmm_regions_t1e-05_interestingHets_IGV.bed")
hist_of_HMM_block_from_bed6(m_all_SRR4041366, p_all_SRR4041366, col= e_col)
#hist_of_HMM_block_from_bed6(m_all_SRR4041366, p_all_SRR4041366, col= e_col)

#GM12878
GM_m_all=read.table("SRR1552485_total/counts_minus_hmm_regions_t1e-05_interestingHets_IGV.bed")
GM_p_all=read.table("SRR1552485_total/counts_plus_hmm_regions_t1e-05_interestingHets_IGV.bed")
hist_of_HMM_block_from_bed6_layer2(GM_m_all, GM_p_all, T, g_col)

legend("topright", 
       legend = c("F1", "GM"),
       #pch=c(15,15),
       cex=1.8, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(10000,30),
       angle=c(180,45),
       #angle=45,
       fill=c(e_col,g_col)
       , bty = "n"
       )

# gene number in each block
##SRR4041366_dedup_2
Em_withG=read.table("SRR4041366_dedup_2/counts_minus_hmm_regions_t1e-05_interestingHets_IGV_genePerBlock.txt")
Ep_withG=read.table("SRR4041366_dedup_2/counts_plus_hmm_regions_t1e-05_interestingHets_IGV_genePerBlock")
EBlock_with_no_gene= dim(m_all_SRR4041366)[1]+dim(p_all_SRR4041366)[1] - dim(Em_withG)[1]-dim(Ep_withG)[1]

#GM12878
GM_m_wGene=read.table("SRR1552485_total/counts_minus_hmm_regions_t1e-05_interestingHets_IGV_genePerBlock.txt")
GM_p_wGene=read.table("SRR1552485_total/counts_plus_hmm_regions_t1e-05_interestingHets_IGV_genePerBlock.txt")
GM_Block_with_no_gene= dim(GM_m_all)[1]+dim(GM_p_all)[1] - dim(GM_m_wGene)[1]-dim(GM_p_wGene)[1]

#par(mar=c(5.1, 6.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
#par(mgp=c(3,1,0))
#par(mgp=c(axis.title.position axis.label.position axis.line.position))
hist(c(rep(0,EBlock_with_no_gene), Em_withG$V1, Ep_withG$V1)
     , breaks = seq(-0.0001,200,1)
     #,density=20,angle=135
     #, freq = F
     , xlim=c(0,15)
     ,ylim=c(0,2500)
     ,xlab="Number of genes in each HMM block"
     ,ylab="AlleleHMM block count"
     ,col= e_col
     ,main=""
     
     ,las=1
     ,bty="l"
)
hist(c(rep(0,GM_Block_with_no_gene), GM_m_wGene$V1, GM_m_wGene$V1)
     , breaks = seq(-0.0001,200,1)
     , xlim=c(0,20)
     ,density=30,angle=45
     ,col= g_col
     ,add=T

)
legend("topright", 
       legend = c("F1", "GM"),
       #pch=c(15,15),
       cex=1.8, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=c(10000,30),
       angle=c(180,45),
       #angle=45,
       fill=c(e_col,g_col)
       , bty = "n"
)
