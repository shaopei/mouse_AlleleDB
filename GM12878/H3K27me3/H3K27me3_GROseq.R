setwd("~/Box Sync/Danko_lab_work/mouse_AlleleDB/GM12878/H3K27me3")

m_input_f="H3K27me3_minus_GROseq_hmm_regions_t1e-05.merged_cov_binomtest.bed"
p_input_f="H3K27me3_plus_GROseq_hmm_regions_t1e-05.merged_cov_binomtest.bed"
m_input_f="H3K27me3_minus_GROseqSig_hmm_regions_t1e-05.merged_cov_binomtest.bed"
p_input_f="H3K27me3_plus_GROseqSig_hmm_regions_t1e-05.merged_cov_binomtest.bed"



mStrand_count = read.table(m_input_f,header=T,sep = "\t",comment.char = "\\")
pStrand_count = read.table(p_input_f,header=T,sep = "\t",comment.char = "\\")

dummy = 0.1
mStrand_count$H_mat_pat_ratio = with(mStrand_count, log2((mat_allele_count+dummy)/(pat_allele_count+dummy)))
mStrandGROseq <- data.frame(do.call('rbind', strsplit(as.character(mStrand_count$GROseq_hmm_state),',',fixed=TRUE)))
names(mStrandGROseq) <- c("hmm_state","pat_count","mat_count")
mStrandGROseq$pat_count <- as.numeric(as.character(mStrandGROseq$pat_count))
mStrandGROseq$mat_count <- as.numeric(as.character(mStrandGROseq$mat_count))
mStrand_count$G_mat_pat_ratio = with(mStrandGROseq, log2((mat_count+dummy)/(pat_count+dummy)))
mStrand_count$H_mat_pat_sum = with(mStrand_count, mat_allele_count+pat_allele_count)
mStrand_count$G_mat_pat_sum = with(mStrandGROseq, mat_count+pat_count)

pStrand_count$H_mat_pat_ratio = with(pStrand_count, log2((mat_allele_count+dummy)/(pat_allele_count+dummy)))
pStrandGROseq <- data.frame(do.call('rbind', strsplit(as.character(pStrand_count$GROseq_hmm_state),',',fixed=TRUE)))
names(pStrandGROseq) <- c("hmm_state","pat_count","mat_count")
pStrandGROseq$pat_count <- as.numeric(as.character(pStrandGROseq$pat_count))
pStrandGROseq$mat_count <- as.numeric(as.character(pStrandGROseq$mat_count))
pStrand_count$G_mat_pat_ratio = with(pStrandGROseq, log2((mat_count+dummy)/(pat_count+dummy)))
pStrand_count$H_mat_pat_sum = with(pStrand_count, mat_allele_count+pat_allele_count)
pStrand_count$G_mat_pat_sum = with(pStrandGROseq, mat_count+pat_count)


#newdata1 <-pStrand_count[order(pStrand_count$Binom_p_value,decreasing = T),]
#newdata1$col=ifelse(newdata1$Binom_p_value < 0.05, "red", "pink")
with (pStrand_count, plot(G_mat_pat_ratio[H_mat_pat_sum >=1 & G_mat_pat_sum >=5], H_mat_pat_ratio[H_mat_pat_sum >=1 & G_mat_pat_sum >=5],
                          xlim=c(-10,10), ylim=c(-10,10),
                          xlab="GRO-seq log2(mat/pat)",
                          ylab="H3K27me3 ChIP-seq log2(mat/pat)",
                     col="pink",
                     #col=colfunc(dim(newdata1)[1]),
                     #pch=20
                     ))
colfunc <- colorRampPalette(c("red" , "pink" ))
#plot(rep(1,10),col=colfunc(10),pch=19,cex=3)
#with (pStrand_count, points(G_mat_pat_ratio[H_mat_pat_sum >=1 & G_mat_pat_sum >=5], H_mat_pat_ratio[H_mat_pat_sum >=1 & G_mat_pat_sum >=5],
#                          col=colfunc(3)[2],
#                          pch=20
#))

colfunc <- colorRampPalette(c("blue" , "light blue" ))
with (mStrand_count, points(G_mat_pat_ratio[H_mat_pat_sum >=1 & G_mat_pat_sum >=5], H_mat_pat_ratio[H_mat_pat_sum >=1 & G_mat_pat_sum >=5], 
                            #pch=20, 
                            col="light blue"))
#with (mStrand_count, points(G_mat_pat_ratio[H_mat_pat_sum >=1 & G_mat_pat_sum >=5], H_mat_pat_ratio[H_mat_pat_sum >=1 & G_mat_pat_sum >=5], 
#                            pch=20, col=colfunc(3)[2]))

fdr10_p=0.011292  #GRO-seq Sig region
fdr10_p=0.005  # GRO-seq all HMM regions
with (pStrand_count, points(G_mat_pat_ratio[H_mat_pat_sum >=1 & G_mat_pat_sum >=5 & Binom_p_value <= fdr10_p], H_mat_pat_ratio[H_mat_pat_sum >=1 & G_mat_pat_sum >=5& Binom_p_value <= fdr10_p],
                     col="red",
                     #pch=20
                     ))
with (mStrand_count, points(G_mat_pat_ratio[H_mat_pat_sum >=1 & G_mat_pat_sum >=5 & Binom_p_value <= fdr10_p], H_mat_pat_ratio[H_mat_pat_sum >=1 & G_mat_pat_sum >=5& Binom_p_value <= fdr10_p],
                            col="blue",
                            #pch=20
                            ))


#colfunc <- colorRampPalette(c("blue" , "light blue" ))
#newdata <-mStrand_count[order(mStrand_count$Binom_p_value),]
#with (newdata, points(G_mat_pat_ratio[H_mat_pat_sum >=1 & G_mat_pat_sum >=5], H_mat_pat_ratio[H_mat_pat_sum >=1 & G_mat_pat_sum >=5], 
#                            pch=20, col=colfunc(dim(newdata)[1])))

abline(v=0)
abline(h=0)

newdatap=with(pStrand_count, pStrand_count[H_mat_pat_sum >=1 & G_mat_pat_sum >=5 & Binom_p_value <= fdr10_p,])
newdatam=with(mStrand_count, mStrand_count[H_mat_pat_sum >=1 & G_mat_pat_sum >=5 & Binom_p_value <= fdr10_p,])
newdata=rbind(newdatam, newdatap)
with (newdata, points(G_mat_pat_ratio, H_mat_pat_ratio,
                            col="black",
                            pch=20))
with(newdata, cor(G_mat_pat_ratio, H_mat_pat_ratio, method = "spearman"))
with(newdata, cor(G_mat_pat_ratio, H_mat_pat_ratio, method = "pearson"))

newdata.lm = lm(G_mat_pat_ratio ~ H_mat_pat_ratio, data=newdata)
summary(newdata.lm)



cor(pStrand_count$G_mat_pat_ratio, pStrand_count$H_mat_pat_ratio, method = "spearman")
pStrand_count.lm = lm(G_mat_pat_ratio ~ H_mat_pat_ratio, data=pStrand_count)
summary(pStrand_count.lm)
mStrand_count.lm = lm(G_mat_pat_ratio ~ H_mat_pat_ratio, data=mStrand_count)
summary(mStrand_count.lm)

hist(c(pStrand_count$H_mat_pat_ratio, mStrand_count$H_mat_pat_ratio))
hist(c(pStrand_count$G_mat_pat_ratio, mStrand_count$G_mat_pat_ratio))
