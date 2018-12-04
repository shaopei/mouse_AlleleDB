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
#od <- c(0.01, 0, 0.01)

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
 blockList$AlleleDB_state[blockList$AlleleDB_p_value <0.05 & blockList$SimState!="S" & blockList$mat_reads/blockList$total_reads >= 0.5] <- "M"
 blockList$AlleleDB_state[blockList$AlleleDB_p_value <0.05 & blockList$SimState!="S" & blockList$mat_reads/blockList$total_reads < 0.5] <- "P" 
 
 
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
 ##test Sensitivity
 if (mat_p[2] >0.5){
   if (sum(blockList$SimState!="S") !=0) {
     SNP_Sensitivity = sum(blockList$AlleleDB_state=="M" & blockList$SimState=="M") /sum(blockList$SimState =="M") 
     AlleleHMM_Sensitivity = sum(blockList$HMM_p_value <0.05 & blockList$HMM_state=="M" & blockList$SimState=="M") /sum(blockList$SimState=="M") 
   }  else{
     SNP_Sensitivity = sum(blockList$AlleleDB_state=="M" & blockList$SimState=="M") /(sum(blockList$SimState=="M")+1) 
     AlleleHMM_Sensitivity = sum(blockList$HMM_p_value <0.05 & blockList$HMM_state=="M" & blockList$SimState=="M") /(sum(blockList$SimState=="M")+1) 
   }
 }
 else{
   if (sum(blockList$SimState!="S") !=0) {
     SNP_Sensitivity = sum(blockList$AlleleDB_state=="P" & blockList$SimState=="P") /sum(blockList$SimState =="P") 
     AlleleHMM_Sensitivity = sum(blockList$HMM_p_value <0.05 & blockList$HMM_state=="P" & blockList$SimState=="P") /sum(blockList$SimState=="P") 
   }  else{
     SNP_Sensitivity = sum(blockList$AlleleDB_state=="P" & blockList$SimState=="P") /(sum(blockList$SimState=="P")+1) 
     AlleleHMM_Sensitivity = sum(blockList$HMM_p_value <0.05 & blockList$HMM_state=="P" & blockList$SimState=="P") /(sum(blockList$SimState=="P")+1) 
   }
 }
 # specificity
    SNP_Specificity = sum(blockList$AlleleDB_p_value >=0.05 & blockList$SimState=="S") /sum(blockList$SimState=="S") 
   AlleleHMM_Specificity = sum(blockList$HMM_state=="S" & blockList$SimState=="S") /sum(blockList$SimState=="S") 
 #print(AlleleHMM_Specificity)
   
  return (list(SNP_Sensitivity=SNP_Sensitivity, 
               AlleleHMM_Sensitivity=AlleleHMM_Sensitivity,
               SNP_Specificity = SNP_Specificity,
               AlleleHMM_Specificity = AlleleHMM_Specificity
               #, data=blockList 
               ))
}

Sensitivity_betaBino_iter<- function(iteration,l, e, mat_p, od, t=1e-5){
  SNP_sen_list=c()
  HMM_sen_list=c()
  SNP_spec_list=c()
  HMM_spec_list=c()
  for (i in 1:iteration) {
    #print (i)
    s=Sensitivity_betaBino(l, e, mat_p, od, t)
    SNP_sen_list=c(SNP_sen_list, s$SNP_Sensitivity)
    HMM_sen_list=c(HMM_sen_list, s$AlleleHMM_Sensitivity)
    SNP_spec_list=c(SNP_spec_list, s$SNP_Specificity)
    HMM_spec_list=c(HMM_spec_list, s$AlleleHMM_Specificity)
  }
  return (c(mean(SNP_sen_list),mean(HMM_sen_list),
            sd(SNP_sen_list)/sqrt(length(SNP_sen_list)),
            sd(HMM_sen_list)/sqrt(length(HMM_sen_list)),
            mean(SNP_spec_list),mean(HMM_spec_list),
            sd(SNP_spec_list)/sqrt(length(SNP_spec_list)),
            sd(HMM_spec_list)/sqrt(length(HMM_spec_list))))
  #return (c(mean(SNP_sen_list), mean(HMM_sen_list)))
}
