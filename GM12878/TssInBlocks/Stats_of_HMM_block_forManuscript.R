#functions

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

hist_of_HMM_block_from_bed6 <- function(m_all, p_all, add=F, col=mycol1){
  mp_all = rbind.data.frame(m_all, p_all)
  mp_all$blockSize = mp_all$V3 - mp_all$V2
  hist(log10(mp_all$blockSize),
       #breaks = seq(0,1e+05,1000),
       freq=F,col = col,
       xlab="log10 of HMM block size (bp)",
       main="The size distribution of HMM blocks",add=add)
}



#GM12878
setwd("~/Box Sync/Danko_lab_work/mouse_AlleleDB/GM12878/TssInBlocks/")
m_all=read.table("counts_minus_hmm_regions_t1e-05_interestingHets_IGV.bed")
p_all=read.table("counts_plus_hmm_regions_t1e-05_interestingHets_IGV.bed")

hist_of_HMM_block_from_bed6(m_all, p_all, T, t_col("black", perc = 50))
#mp_all = rbind.data.frame(m_all, p_all)
#mp_all$blockSize = mp_all$V3 - mp_all$V2
#hist(log10(mp_all$blockSize),
#     freq=F,col = mycol1,
#     xlab="log10 of HMM block size (bp)",
#     main="The size distribution of HMM blocks")

#Engreitz_F1
m_all_SRR4041366 = read.table("~/Box Sync/Danko_lab_work/mouse_AlleleDB/Engreitz_F1_20180523/allelicbias-PersonalGenome_P.CAST_M.129S1-SRR4041366_dedup_2/SRR4041366_dedup_2_counts_minus_hmm_regions_t1e-05_interestingHets_IGV.bed")
p_all_SRR4041366 = read.table("~/Box Sync/Danko_lab_work/mouse_AlleleDB/Engreitz_F1_20180523/allelicbias-PersonalGenome_P.CAST_M.129S1-SRR4041366_dedup_2/SRR4041366_dedup_2_counts_plus_hmm_regions_t1e-05_interestingHets_IGV.bed")
hist_of_HMM_block_from_bed6(m_all_SRR4041366, p_all_SRR4041366, add=T, col=mycol2)
hist_of_HMM_block_from_bed6(m_all_SRR4041366, p_all_SRR4041366)

m_all_SRR4041367 = read.table("~/Box Sync/Danko_lab_work/mouse_AlleleDB/Engreitz_F1_20180523/allelicbias-PersonalGenome_P.CAST_M.129S1-SRR4041367_dedup_2/SRR4041367_dedup_2_counts_minus_hmm_regions_t1e-05_interestingHets_IGV.bed")
p_all_SRR4041367 = read.table("~/Box Sync/Danko_lab_work/mouse_AlleleDB/Engreitz_F1_20180523/allelicbias-PersonalGenome_P.CAST_M.129S1-SRR4041367_dedup_2/SRR4041367_dedup_2_counts_plus_hmm_regions_t1e-05_interestingHets_IGV.bed")
hist_of_HMM_block_from_bed6(m_all_SRR4041367, p_all_SRR4041367, add=T, col=mycol1)

m_all_SRR4041367_2R = read.table("~/Box Sync/Danko_lab_work/mouse_AlleleDB/Engreitz_F1_20180523/allelicbias-PersonalGenome_P.CAST_M.129S1-SRR4041367_dedup_2_and_RC1/SRR4041367_dedup_2_and_RC1_counts_minus_hmm_regions_t1e-05_interestingHets_IGV.bed")
p_all_SRR4041367_2R = read.table("~/Box Sync/Danko_lab_work/mouse_AlleleDB/Engreitz_F1_20180523/allelicbias-PersonalGenome_P.CAST_M.129S1-SRR4041367_dedup_2_and_RC1/SRR4041367_dedup_2_and_RC1_counts_plus_hmm_regions_t1e-05_interestingHets_IGV.bed")
hist_of_HMM_block_from_bed6(m_all_SRR4041367, p_all_SRR4041367_2R, add=T, col=mycol3)

# LEP_ZYG
m_all_LZ = read.table("~/Box Sync/Danko_lab_work/mouse_AlleleDB/Engreitz_F1_20180523/allelicbias-PersonalGenome_P.CAST_M.B6-LEP_ZYG_ATGCA_forAlleleDB/LEP_ZYG_counts_minus_hmm_regions_t1e-05_interestingHets_IGV.bed")
p_all_LZ = read.table("~/Box Sync/Danko_lab_work/mouse_AlleleDB/Engreitz_F1_20180523/allelicbias-PersonalGenome_P.CAST_M.B6-LEP_ZYG_ATGCA_forAlleleDB/LEP_ZYG_counts_plus_hmm_regions_t1e-05_interestingHets_IGV.bed")
hist_of_HMM_block_from_bed6(m_all_LZ, p_all_LZ, F, "orange")
