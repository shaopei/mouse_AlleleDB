## wrap into a function to test sensitivity 
#Simulation of one HMM block using Beta-Binomial, instead of Binomial

#BlockLength=20
#ExpressionLevel=3
#MatBinoP=0.1
#rbetabinom(n, size, m, s)
#rbinom(n, size, prob)
library(VGAM)
sim_HMM_block_betaBino <- function(BlockLength, ExpressionLevel, MatBinoP, OD){
  total_reads=rpois(BlockLength, ExpressionLevel)
  mat_reads=rbetabinom(BlockLength, total_reads, MatBinoP, OD)
  SimState = rep(ifelse(MatBinoP<0.49, "P", ifelse(MatBinoP>0.51,"M", "S")),BlockLength)
  MatP= rep(round(MatBinoP, digits = 3),BlockLength)
  Rho= rep(round(OD, digits = 3),BlockLength)
  block=data.frame(total_reads, mat_reads,SimState,MatBinoP=MatP,Rho)
  return(block)
}


#b1=sim_HMM_block_betaBino(10,5, 0.9, 0.01)
#b1


#numberOfBlock <- 3
#legnth of each block
#l <- c(10,100,10)
# expression level of each block
#e <- c(10,10,10)
#mat_p <- c(0.5,0.9,0.5)
#od <- c(0.25, 0.25, 0.25)

Sensitivity_betaBino<- function(l, e, mat_p, od, t){
  # Simulation of 3 HMM blocks (Sym-Mat-Sym)
  i=1
  blockList=data.frame(sim_HMM_block_betaBino(l[i],e[i], mat_p[i], od[i]), blockID=i)
  for (i in 2:numberOfBlock){
    blockList= rbind.data.frame(blockList, cbind(sim_HMM_block_betaBino(l[i],e[i], mat_p[i], od[i]),blockID=i))
  }
  #blockList=blockList[blockList$total_reads>0,c("blockID", "SimState", "MatBinoP", "total_reads",  "mat_reads")]
  blockList=blockList[blockList$total_reads>0,c("blockID", "SimState", "MatBinoP","Rho", "total_reads",  "mat_reads")]
  #b_tmp=tempfile(tmpdir ="/workdir/sc2457/HMM_simulation/toremove",fileext = ".txt")
  b_tmp=tempfile(fileext = ".txt")
  #b_tmp="sim_tmp.txt"
  write.table(blockList, file = b_tmp, quote = F, sep = "\t", row.names = F)
  #View(blockList)
  
  # naive model:treating SNPs independently.
  blockList$AlleleDB_p_value <- 1
  i=1
  for (i in 1:dim(blockList)[1]){
    a = binom.test(blockList$mat_reads[i], blockList$total_reads[i], 0.5)
    blockList$AlleleDB_p_value[i] <- a$p.value    
  }

 blockList$AlleleDB_state <- "S"
 blockList$AlleleDB_state[blockList$AlleleDB_p_value <0.05 & blockList$mat_reads/blockList$total_reads > 0.5] <- "M"
 blockList$AlleleDB_state[blockList$AlleleDB_p_value <0.05 & blockList$mat_reads/blockList$total_reads < 0.5] <- "P" 
 
 
  ## HMM prediction
  system(paste("python HMM_prediction.py","-i",b_tmp, "-t", t), wait=TRUE)
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

### calculate TP, FP, TN, N  # includes all SNPs at three blocks
  if (mat_p[2] >0.5){
     P=sum(blockList$SimState =="M") 
     SNP_TP = sum(blockList$AlleleDB_p_value < 0.05 & blockList$AlleleDB_state=="M" & blockList$SimState=="M")
     HMM_TP = sum(blockList$HMM_p_value <0.05 & blockList$HMM_state=="M" & blockList$SimState=="M") 
 } else{
   P=sum(blockList$SimState =="P") 
   SNP_TP = sum(blockList$AlleleDB_p_value < 0.05 & blockList$AlleleDB_state=="P" & blockList$SimState=="P")
   HMM_TP = sum(blockList$HMM_p_value <0.05 & blockList$HMM_state=="P" & blockList$SimState=="P")
   }

 SNP_FP = sum(blockList$AlleleDB_p_value <0.05 & blockList$SimState =="S")   #incldue padding blocks
 HMM_FP = sum(blockList$HMM_p_value <0.05 & blockList$HMM_state !="S" & blockList$SimState=="S") #incldue padding blocks
 
 N=sum(blockList$SimState=="S")   #incldue padding blocks
 
SNP_TN=sum(blockList$AlleleDB_p_value >=0.05 & blockList$SimState=="S")   #incldue padding blocks
HMM_TN=sum(blockList$HMM_state=="S" & blockList$SimState=="S")   #incldue padding blocks

 
##test Sensitivity
 if (P !=0) {
   SNP_Sensitivity = SNP_TP/P
   AlleleHMM_Sensitivity =HMM_TP /P
 }  else{
    SNP_Sensitivity = 0
   AlleleHMM_Sensitivity = 0
 }
 
 # specificity
    SNP_Specificity = SNP_TN /N
   AlleleHMM_Specificity = HMM_TN /N
 
  return (list(SNP_Sensitivity=SNP_Sensitivity, 
               AlleleHMM_Sensitivity=AlleleHMM_Sensitivity,
               SNP_Specificity = SNP_Specificity,
               AlleleHMM_Specificity = AlleleHMM_Specificity,
               P_Counts=c(P,SNP_TP, HMM_TP, SNP_FP, HMM_FP),
               N_Counts=c(N, SNP_TN, HMM_TN)
               #, data=blockList 
               ))
}



Precision_recall_specificity_betaBino_iter<- function(iteration,l, e, mat_p, od, t=1e-5){
  SNP_sen_list=c()
  HMM_sen_list=c()
  SNP_spec_list=c()
  HMM_spec_list=c()
  SNP_prec_list=c()
  HMM_prec_list=c()
  
  for (i in 1:iteration) {
    #print (i)
    m = Sensitivity_betaBino(l, e, mat_p, od, t)  #SMS blocks
    s = Sensitivity_betaBino(l, e, c(0.5,0.5,0.5), od, t) #SSS blocks
    P = m$P_Counts[1] #SMS    
    SNP_TP = m$P_Counts[2] #SMS
    HMM_TP = m$P_Counts[3] #SMS  
    
    SNP_FP = s$P_Counts[4]  #SSS
    HMM_FP = s$P_Counts[5]  #SSS
    N = s$N_Counts[1] #SSS
    SNP_TN = s$N_Counts[2] #SSS
    HMM_TN = s$N_Counts[3] #SSS

 ## Sensitivity TP/P
 if (P !=0) {
   SNP_Sensitivity = SNP_TP/P
   AlleleHMM_Sensitivity =HMM_TP /P
 }  else{
    SNP_Sensitivity = 0
   AlleleHMM_Sensitivity = 0
 }
 
 # precision = TP/TP+FP
    if (SNP_TP == 0){
      SNP_precision = 0
    } else{
      SNP_precision = SNP_TP /(SNP_TP+SNP_FP)
    }
    if (HMM_TP == 0){
      AlleleHMM_precision = 0
    } else{
      AlleleHMM_precision = HMM_TP /(HMM_TP+HMM_FP)
    }
    
 # specificity
  SNP_Specificity = SNP_TN /N
  AlleleHMM_Specificity = HMM_TN /N
 #print(AlleleHMM_Specificity)


    
    SNP_sen_list=c(SNP_sen_list, SNP_Sensitivity)
    HMM_sen_list=c(HMM_sen_list, AlleleHMM_Sensitivity)
    SNP_spec_list=c(SNP_spec_list, SNP_Specificity)
    HMM_spec_list=c(HMM_spec_list, AlleleHMM_Specificity)
    SNP_prec_list=c(SNP_prec_list, SNP_precision)
    HMM_prec_list=c(HMM_prec_list, AlleleHMM_precision)
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


## test multi-threads with MatP

numberOfBlock <- 3
#legnth of each block
l <- c(10,100,10)
# expression level of each block
e <- c(10,10,10)
mat_p <- c(0.5,0.9,0.5)
od <- c(0.25, 0.25, 0.25)

library(parallel)
m_test <- function(mat_p_2, iteration){
  #mat_p_list = c(mat_p_list, mat_p_2) 
  mat_p=c(0.5,mat_p_2,0.5)
  #print(mat_p)
  aveS = Precision_recall_specificity_betaBino_iter(iteration,l, e, mat_p, od)
  return(c(mat_p_2,aveS))
}

test=mclapply(seq(0,1,0.05),FUN=function(idx){m_test(idx,1000)}, mc.cores = 50)


m=Get_table_from_simulation(test)
write.table(m, file = "MatP_SMS_OD0.25_e10_l100_sensitivity_and_precision_iter1K.txt", quote = F, sep = "\t",row.names = F)

#m=read.table("MatP_SMS_OD0_e10_l100_sensitivity_and_precision_iter1K.txt", header=T,sep = "\t")
library(ggplot2)
pdf("MatP_SMS_OD0.25_e10_l100_sensitivity_iter1k.pdf")
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
pdf("MatP_SSS_OD0.25_e10_l100_specificity_iter1k.pdf")
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


pdf("MatP_SMS_OD0.25_e10_l100_precision_iter1K.pdf")
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


pdf("MatP_SMS_OD0.25_e10_l100_precision_recall_iter1K.pdf")
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



### expression level###
numberOfBlock <- 3
#legnth of each block
l <- c(10,100,10)
# expression level of each block
e <- c(10,10,10)
mat_p <- c(0.5,0.9,0.5)
od <- c(0.25, 0.25, 0.25)

library(parallel)
e_test <- function(e_2, iteration){
  e=c(10, e_2, 10)
  aveS = Precision_recall_specificity_betaBino_iter(iteration,l, e, mat_p, od, t=1e-5) 
  return(c(e_2,aveS))
}
e_test_out=mclapply(seq(1,50,1),FUN=function(idx){e_test(idx,1000)}, mc.cores = 50)

m=Get_table_from_simulation(e_test_out)
write.table(m, file = "exp_OD0.25_l100_sensitivity_and_precision_iter1K.txt", quote = F, sep = "\t",row.names = F)

#m=read.table("MatP_SMS_OD0_e10_l100_sensitivity_and_precision_iter1K.txt", header=T,sep = "\t")
library(ggplot2)
pdf("exp_OD0.25_l100_sensitivity_iter1k.pdf")
ggplot(m, aes(x=trait, y=senMean, colour=method)) + 
    geom_errorbar(aes(ymin=senMean-senSE, ymax=senMean+senSE), width=.01) +
    geom_line() +
    geom_point() +    
   scale_color_manual(values=c("red", "blue")) +
    theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"))+
  xlab("read count")
    #theme_bw()
dev.off()
pdf("exp_OD0.25_e10_l100_specificity_iter1k.pdf")
ggplot(m, aes(x=trait, y=specMean, colour=method)) + 
    geom_errorbar(aes(ymin=specMean-specSE, ymax=specMean+specSE), width=.01) +
    geom_line() +
    geom_point()+    
   scale_color_manual(values=c("red", "blue")) +
    theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"))+
  xlab("read count")
    #theme_bw()
dev.off()


pdf("exp_OD0.25_l100_precision_iter1K.pdf")
ggplot(m, aes(x=trait, y=precMean, colour=method)) + 
    geom_errorbar(aes(ymin=precMean-precSE, ymax=precMean+precSE), width=.01) +
    geom_line() +
    geom_point() +
    scale_color_manual(values=c("red", "blue")) +
    theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"))+
  xlab("read count")
    #theme_bw()
dev.off()


pdf("exp_OD0.25_l100_precision_recall_iter1K.pdf")
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


### length ####
numberOfBlock <- 3
#legnth of each block
l <- c(10,100,10)
# expression level of each block
e <- c(10,10,10)
mat_p <- c(0.5,0.9,0.5)
od <- c(0.25, 0.25, 0.25)

l_test <- function(l_2, iteration){
  l=c(10, l_2, 10)
  aveS = Precision_recall_specificity_betaBino_iter(iteration,l, e, mat_p, od)
  return(c(l_2,aveS))
}

l_test_out=mclapply(seq(1,50,1),FUN=function(idx){l_test(idx,1000)}, mc.cores = 50)


m=Get_table_from_simulation(l_test_out)
write.table(m, file = "length_SMS_OD0.25_e10_sensitivity_and_precision_iter1K.txt", quote = F, sep = "\t",row.names = F)

library(ggplot2)
pdf("length_SMS_OD0.25_e10_sensitivity_iter1K.pdf")
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
pdf("length_SMS_OD0.25_e10_specificity_iter1K.pdf")
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

pdf("length_SMS_OD0.25_e10_precision_iter1K.pdf")
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

pdf("length_SMS_OD0.25_e10_precision_recall_iter1K.pdf")
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


## test multi-threads with MatP

numberOfBlock <- 3
#legnth of each block
l <- c(10,100,10)
# expression level of each block
e <- c(10,10,10)
mat_p <- c(0.5,0.9,0.5)
od <- c(0.1, 0.1, 0.1)

library(parallel)
od_test <- function(od_2, iteration){
  od = c(od_2, od_2, od_2)
  aveS = Precision_recall_specificity_betaBino_iter(iteration,l, e, mat_p, od)
  return(c(od_2,aveS))
}

test=mclapply(seq(0,0.9,0.05),FUN=function(idx){od_test(idx,iteration = 1000)}, mc.cores = 30)
m=Get_table_from_simulation(test)
write.table(m, file = "OD_TEST_e10_l100_sensitivity_and_precision_iter1K.txt", quote = F, sep = "\t",row.names = F)

