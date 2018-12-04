## wrap into a function to test sensitivity 

#Simulation of one HMM block
#BlockLength=20
#ExpressionLevel=3
#MatBinoP=0.1
sim_HMM_block <- function(BlockLength, ExpressionLevel, MatBinoP){
  total_reads=rpois(BlockLength, ExpressionLevel)
  mat_reads=rbinom(BlockLength, total_reads, MatBinoP)
  SimState = rep(ifelse(MatBinoP<0.49, "P", ifelse(MatBinoP>0.51,"M", "S")),BlockLength)
  MatP= rep(round(MatBinoP, digits = 3),BlockLength)
  block=data.frame(total_reads, mat_reads,MatBinoP=MatP,SimState)
  return(block)
}
#b1=sim_HMM_block(10,5, 0.2)
#b1


numberOfBlock <- 3
#legnth of each block
l <- c(10,100,10)
# expression level of each block
e <- c(10,10,10)
mat_p <- c(0.5,0.9,0.5)

Sensitivity<- function(l, e, mat_p){
  # Simulation of 3 HMM blocks (Sym-Mat-Sym)
  i=1
  blockList=data.frame(sim_HMM_block(l[i],e[i], mat_p[i]), blockID=i)
  for (i in 2:numberOfBlock){
    blockList= rbind.data.frame(blockList, cbind(sim_HMM_block(l[i],e[i], mat_p[i]),blockID=i))
  }
  blockList=blockList[blockList$total_reads>0,c("blockID", "SimState", "MatBinoP", "total_reads",  "mat_reads")]
  #b_tmp=tempfile(tmpdir ="/workdir/sc2457/HMM_simulation/toremove",fileext = ".txt")
  b_tmp=tempfile(fileext = ".txt")
  write.table(blockList, file = b_tmp, quote = F, sep = "\t")
  #View(blockList)
  
  # naive model:treating SNPs independently.
  blockList$AlleleDB_p_value <- 1
  i=1
  for (i in 1:dim(blockList)[1]){
    a = binom.test(blockList$mat_reads[i], blockList$total_reads[i], 0.5)
    blockList$AlleleDB_p_value[i] <- a$p.value    
  }

 ## HMM prediction
  system(paste("python HMM_prediction.py","-i",b_tmp), wait=TRUE)
 ## HMM binomial test
 
 # combine the reads in the same HMM blocks and perform bionomial test
 hmm_result=read.table(paste(substr(b_tmp,1,nchar(b_tmp)-4),"_AlleleHMM.txt",sep = ""),header=T)
 hmm_result$tag=0
 hmm_result$accum_total_reads=0
 hmm_result$accum_mat_reads=0
 hmm_result$p_value = 1
 i=1
 accum_total_reads=hmm_result$total_reads[i]
 accum_mat_reads=hmm_result$mat_reads[i]
 hmm_result$tag[i]=1
 for (i in 2:dim(hmm_result)[1]){
   if (hmm_result$hmm_state[i] == hmm_result$hmm_state[i-1]){
     hmm_result$tag[i]=1
     accum_total_reads <- accum_total_reads + hmm_result$total_reads[i]
     accum_mat_reads <- accum_mat_reads + hmm_result$mat_reads[i]
   }
   else{
     hmm_result$accum_total_reads[hmm_result$tag==1] <- accum_total_reads
     hmm_result$accum_mat_reads[hmm_result$tag==1] <- accum_mat_reads
     a = binom.test(accum_mat_reads, accum_total_reads, 0.5)
     hmm_result$p_value[hmm_result$tag==1] <- a$p.value
     
     hmm_result$tag=0
     hmm_result$tag[i]=1
     accum_total_reads <- hmm_result$total_reads[i]
     accum_mat_reads   <- hmm_result$mat_reads[i]
   }
 }
 hmm_result$accum_total_reads[hmm_result$tag==1] <- accum_total_reads
 hmm_result$accum_mat_reads[hmm_result$tag==1] <- accum_mat_reads
 a = binom.test(accum_mat_reads, accum_total_reads, 0.5)
 hmm_result$p_value[hmm_result$tag==1] <- a$p.value
 #print(dim(blockList)[1] == dim(hmm_result)[1])
 
 blockList$HMM_p_value <- hmm_result$p_value
 blockList$HMM_accum_total_reads <- hmm_result$accum_total_reads
 blockList$HMM_accum_mat_reads<- hmm_result$accum_mat_reads
 blockList$HMM_state<- hmm_result$hmm_state
 ##test Sensitivity
 if (sum(blockList$SimState!="S") !=0) {
   SNP_Sensitivity = sum(blockList$AlleleDB_p_value <0.05 & blockList$SimState!="S") /sum(blockList$SimState!="S") 
   AlleleHMM_Sensitivity = sum(blockList$HMM_p_value <0.05 & blockList$SimState!="S") /sum(blockList$SimState!="S") 
 }  else{
    SNP_Sensitivity = sum(blockList$AlleleDB_p_value <0.05 & blockList$SimState!="S") /(sum(blockList$SimState!="S")+1) 
   AlleleHMM_Sensitivity = sum(blockList$HMM_p_value <0.05 & blockList$SimState!="S") /(sum(blockList$SimState!="S")+1) 
 }
 
 # precision = TP/TP+FP
  if (sum(blockList$AlleleDB_p_value <0.05) !=0) {
   SNP_precision = sum(blockList$AlleleDB_p_value <0.05 & blockList$SimState!="S") /sum(blockList$AlleleDB_p_value <0.05)
  }  else{   
      SNP_precision = sum(blockList$AlleleDB_p_value <0.05 & blockList$SimState!="S") /(sum(blockList$AlleleDB_p_value <0.05)+1)
  }
 if (sum(blockList$HMM_p_value <0.05) !=0){
   AlleleHMM_precision = sum(blockList$HMM_p_value <0.05 & blockList$SimState!="S") /sum(blockList$HMM_p_value <0.05)
 } else {
      AlleleHMM_precision = sum(blockList$HMM_p_value <0.05 & blockList$SimState!="S") /(sum(blockList$HMM_p_value <0.05)+1)
  }
 
 
 
 
 # specificity
    SNP_Specificity = sum(blockList$AlleleDB_p_value >=0.05 & blockList$SimState=="S") /sum(blockList$SimState=="S") 
   AlleleHMM_Specificity = sum(blockList$HMM_state=="S" & blockList$SimState=="S") /sum(blockList$SimState=="S") 
 #print(AlleleHMM_Specificity)
   
  return (list(SNP_Sensitivity=SNP_Sensitivity, 
               AlleleHMM_Sensitivity=AlleleHMM_Sensitivity,
               SNP_Specificity = SNP_Specificity,
               AlleleHMM_Specificity = AlleleHMM_Specificity,
               SNP_precision = SNP_precision,
               AlleleHMM_precision = AlleleHMM_precision 
               #, data=blockList 
               ))
}

Sensitivity_iter<- function(iteration,l, e, mat_p){
  SNP_sen_list=c()
  HMM_sen_list=c()
  SNP_spec_list=c()
  HMM_spec_list=c()
  SNP_prec_list=c()
  HMM_prec_list=c()
  for (i in 1:iteration) {
    #print (i)
    s=Sensitivity(l, e, mat_p)
    SNP_sen_list=c(SNP_sen_list, s$SNP_Sensitivity)
    HMM_sen_list=c(HMM_sen_list, s$AlleleHMM_Sensitivity)
    SNP_spec_list=c(SNP_spec_list, s$SNP_Specificity)
    HMM_spec_list=c(HMM_spec_list, s$AlleleHMM_Specificity)
    SNP_prec_list=c(SNP_prec_list, s$SNP_precision)
    HMM_prec_list=c(HMM_prec_list, s$AlleleHMM_precision)
  }
  return (c(mean(SNP_sen_list),
            mean(HMM_sen_list),
            sd(SNP_sen_list)/sqrt(length(SNP_sen_list)),
            sd(HMM_sen_list)/sqrt(length(HMM_sen_list)),
            mean(SNP_spec_list),
            mean(HMM_spec_list),
            sd(SNP_spec_list)/sqrt(length(SNP_spec_list)),
            sd(HMM_spec_list)/sqrt(length(HMM_spec_list)),
            mean(SNP_prec_list),
            mean(HMM_prec_list),
            sd(SNP_prec_list)/sqrt(length(SNP_prec_list)),
            sd(HMM_prec_list)/sqrt(length(HMM_prec_list))))
  #return (c(mean(SNP_sen_list), mean(HMM_sen_list)))
}

Get_table_from_simulation <- function(test){

result=unlist(test)
outout_length=12+1

m1=data.frame(result[seq(1,length(result),outout_length)],
              result[seq(2,length(result),outout_length)], 
              result[seq(4,length(result),outout_length)],
              result[seq(6,length(result),outout_length)],
              result[seq(8,length(result),outout_length)],
              result[seq(10,length(result),outout_length)],
              result[seq(12,length(result),outout_length)])
names(m1)=c("trait","senMean", "senSE", "specMean", "specSE", "precMean","precSE")
m1$method="SNP independantly"
#View(m)
m2=data.frame(result[seq(1,length(result),outout_length)],
              result[seq(3,length(result),outout_length)], 
              result[seq(5,length(result),outout_length)],
              result[seq(7,length(result),outout_length)],
              result[seq(9,length(result),outout_length)],
              result[seq(11,length(result),outout_length)],
              result[seq(13,length(result),outout_length)])
names(m2)=c("trait","senMean", "senSE", "specMean", "specSE", "precMean","precSE")
m2$method="Allele HMM"
m=rbind.data.frame(m1,m2)
}

# matP
## test multi-threads with MatP

numberOfBlock <- 3
#legnth of each block
l <- c(10,100,10)
# expression level of each block
e <- c(10,10,10)
mat_p <- c(0.5,0.9,0.5)

library(parallel)
m_test <- function(mat_p_2, iteration){
  #mat_p_list = c(mat_p_list, mat_p_2) 
  mat_p=c(0.5,mat_p_2,0.5)
  #print(mat_p)
  aveS = Sensitivity_iter(iteration,l, e, mat_p)
  return(c(mat_p_2,aveS))
}

test=mclapply(seq(0.1,0.9,0.05),FUN=function(idx){m_test(idx,1000)}, mc.cores = 50)

m=Get_table_from_simulation(test)
write.table(m, file = "MatP_SMS_OD0_e10_l100_sensitivity_and_precision_iter1K.txt", quote = F, sep = "\t",row.names = F)

#m=read.table("MatP_SMS_OD0_e10_l100_sensitivity_and_precision_iter1K.txt", header=T,sep = "\t")
library(ggplot2)
pdf("MatP_SMS_OD0_e10_l100_sensitivity_iter1k.pdf")
ggplot(m, aes(x=trait, y=senMean, colour=method)) + 
    geom_errorbar(aes(ymin=senMean-senSE, ymax=senMean+senSE), width=.01) +
    geom_line() +
    geom_point() +    
   scale_color_manual(values=c("red", "blue")) +
    theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"))+
  xlab("Mat reads fraction")
    #theme_bw()
dev.off()
pdf("MatP_SMS_OD0_e10_l100_specificity_iter1k.pdf")
ggplot(m, aes(x=trait, y=specMean, colour=method)) + 
    geom_errorbar(aes(ymin=specMean-specSE, ymax=specMean+specSE), width=.01) +
    geom_line() +
    geom_point()+    
   scale_color_manual(values=c("red", "blue")) +
    theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"))+
  xlab("Mat reads fraction")
    #theme_bw()
dev.off()


pdf("MatP_SMS_OD0_e10_l100_precision_iter1K.pdf")
ggplot(m, aes(x=trait, y=precMean, colour=method)) + 
    geom_errorbar(aes(ymin=precMean-precSE, ymax=precMean+precSE), width=.01) +
    geom_line() +
    geom_point() +
    scale_color_manual(values=c("red", "blue")) +
    theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"))+
  xlab("Mat reads fraction")
    #theme_bw()
dev.off()


pdf("MatP_SMS_OD0_e10_l100_precision_recall_iter1K.pdf")
ggplot(m, aes(x=senMean, y=precMean, colour=method)) + 
#    geom_errorbar(aes(ymin=precMean-precSE, ymax=precMean+precSE), width=.01) +
    geom_line() +
    geom_point() +
    scale_color_manual(values=c("red", "blue")) +
    theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"))+
    xlab("recall")+
    ylab("precision")+
    ggtitle("Mat reads fraction")
    #theme_bw()
dev.off()

# m=read.table("binomial/SMS/MatP_SMS_OD0_e10_l100_sensitivity_and_precision_iter1K.txt", header=T,sep = "\t")
# mplot<- ggplot(m, aes(x=trait, y=senMean, colour=method)) + 
#   geom_errorbar(aes(ymin=senMean-senSE, ymax=senMean+senSE), width=.05) +
#   geom_line() +
#   geom_point()+
#   xlab("Mat Binomial P")+
#   ylab("Sensitivity")+  
#   ylim(c(0,1)) 

# mplot+ geom_line(size = 1) +  geom_point(size = 4) +
#   theme(text = element_text(size=30)) +
#   theme(legend.position="none")+
#   scale_color_manual(values=c("red", "blue"))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))


## test multi-threads with expression level
numberOfBlock <- 3
#legnth of each block
l <- c(10,100,10)
# expression level of each block
e <- c(10,10,10)
mat_p <- c(0.5,0.9,0.5)

library(parallel)
e_test <- function(e_2, iteration){
  e=c(10, e_2, 10)
  aveS = Sensitivity_iter(iteration,l, e, mat_p)
  return(c(e_2,aveS))
}
e_test_out=mclapply(seq(1,50,1),FUN=function(idx){e_test(idx,1000)}, mc.cores = 50)

m=Get_table_from_simulation(e_test_out)
write.table(m, file = "expression_SMS_OD0_l100_sensitivity_and_precision_iter1K.txt", quote = F, sep = "\t",row.names = F)

pdf("exp_SMS_OD0_l100_precision_iter1K.pdf")
ggplot(m, aes(x=trait, y=precMean, colour=method)) + 
    geom_errorbar(aes(ymin=precMean-precSE, ymax=precMean+precSE), width=.01) +
    geom_line() +
    geom_point() +
    scale_color_manual(values=c("red", "blue")) +
    theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"))+
  xlab("read counts")
    #theme_bw()
dev.off()


pdf("exp_SMS_OD0_l100_precision_recall_iter1K.pdf")
ggplot(m, aes(x=senMean, y=precMean, colour=method)) + 
#    geom_errorbar(aes(ymin=precMean-precSE, ymax=precMean+precSE), width=.01) +
    geom_line() +
    geom_point() +
    scale_color_manual(values=c("red", "blue")) +
    theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"))+
    xlab("recall")+
    ylab("precision")+
    ggtitle("read count")
    #theme_bw()
dev.off()



## test multi-threads with block length
numberOfBlock <- 3
#legnth of each block
l <- c(10,100,10)
# expression level of each block
e <- c(10,10,10)
mat_p <- c(0.5,0.9,0.5)

library(parallel)
l_test <- function(l_2, iteration){
  l=c(10, l_2, 10)
  aveS = Sensitivity_iter(iteration,l, e, mat_p)
  return(c(l_2,aveS))
}

l_test_out=mclapply(seq(1,50,1),FUN=function(idx){l_test(idx,1000)}, mc.cores = 50)


m=Get_table_from_simulation(l_test_out)
write.table(m, file = "length_SMS_OD0_e10_sensitivity_and_precision_iter1K.txt", quote = F, sep = "\t",row.names = F)

library(ggplot2)
pdf("length_SMS_OD0_e10_sensitivity_iter1K.pdf")
ggplot(m, aes(x=trait, y=senMean, colour=method)) + 
    geom_errorbar(aes(ymin=senMean-senSE, ymax=senMean+senSE), width=.01) +
    geom_line() +
    geom_point()+
    scale_color_manual(values=c("red", "blue")) +
    theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"))+
  xlab("length") #+
    #theme_bw()
dev.off()
pdf("length_SMS_OD0_e10_specificity_iter1K.pdf")
ggplot(m, aes(x=trait, y=specMean, colour=method)) + 
    geom_errorbar(aes(ymin=specMean-specSE, ymax=specMean+specSE), width=.01) +
    geom_line() +
    geom_point() +
    scale_color_manual(values=c("red", "blue")) +
    theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"))+
  xlab("length") 
    #theme_bw()
dev.off()

pdf("length_SMS_OD0_e10_precision_iter1K.pdf")
ggplot(m, aes(x=trait, y=precMean, colour=method)) + 
    geom_errorbar(aes(ymin=precMean-precSE, ymax=precMean+precSE), width=.01) +
    geom_line() +
    geom_point() +
    scale_color_manual(values=c("red", "blue")) +
    theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"))+
  xlab("length")
    #theme_bw()
dev.off()

pdf("length_SMS_OD0_e10_precision_recall_iter1K.pdf")
ggplot(m, aes(x=senMean, y=precMean, colour=method)) + 
#    geom_errorbar(aes(ymin=precMean-precSE, ymax=precMean+precSE), width=.01) +
    geom_line() +
    geom_point() +
    scale_color_manual(values=c("red", "blue")) +
    theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"))+
    xlab("recall")+
    ylab("precision")+
    ggtitle("Center Block length")
    #theme_bw()
dev.off()