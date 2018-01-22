# R --vanilla --slave --args $(pwd) Con_Dis_regions.bed < Con_Dis_hmm_scatterplot.R 

#arguments here
args=(commandArgs(TRUE))
setwd(args[1])
input_f = args[2] #Con_Dis_regions.bed
dummy = 1

# code start
#setwd('/Users/shaopei/Box Sync/temp')
Strand_count = read.table(input_f,header=T,sep = "\t")
original_col_number = dim(Strand_count)[2]
summary(Strand_count)
Con_col_1 = 'dark green'#'springgreen4'
#Con_col_2 = 'blue'
Dis_col_1 ='dark orange'
#Dis_col_2 = 'dark red'
filter_out_col = 'gray77'
PlusStrandState <- data.frame(do.call('rbind', strsplit(as.character(Strand_count$PlusStrandState),',',fixed=TRUE)))
names(PlusStrandState) <- c("hmm_state","pat_count","mat_count")
PlusStrandState$pat_count <- as.numeric(as.character(PlusStrandState$pat_count))
PlusStrandState$mat_count <- as.numeric(as.character(PlusStrandState$mat_count))
MinusStrandState <- data.frame(do.call('rbind', strsplit(as.character(Strand_count$MinusStrandState),',',fixed=TRUE)))
names(MinusStrandState) <- c("hmm_state","pat_count","mat_count")
MinusStrandState$pat_count <- as.numeric(as.character(MinusStrandState$pat_count))
MinusStrandState$mat_count <- as.numeric(as.character(MinusStrandState$mat_count))
Strand_count$log_mat_over_pat_plus =with(PlusStrandState,log2((mat_count+dummy)/(pat_count+dummy)))
Strand_count$log_mat_over_pat_minus =with(MinusStrandState,log2((mat_count+dummy)/(pat_count+dummy)))
# Stats about the number of regions
d = (Strand_count$log_mat_over_pat_plus *Strand_count$log_mat_over_pat_minus < 0) &(Strand_count$Con_Dis=='Dis')
c = (Strand_count$log_mat_over_pat_plus *Strand_count$log_mat_over_pat_minus > 0) & (Strand_count$Con_Dis == 'Con')
#Strand_count=Strand_count[d|c,]
# the ratio between reads count on minus and plus strand
Strand_count$log2_minus_over_plus = log2(Strand_count$log_mat_over_pat_minus / Strand_count$log_mat_over_pat_plus)

# Stats about the number of regions
# Concordant regions are filtered by log2_minus_over_plus
f = sum((Strand_count$log2_minus_over_plus[Strand_count$Con_Dis == 'Con'] < -1.5) | 
          (Strand_count$log2_minus_over_plus[Strand_count$Con_Dis == 'Con'] > 1.5), 
        na.rm = T)
cat('# of Concordant region before filter = ', sum(c),'\n')
cat('# of Concordant region after filter = ', sum(c)-f,'\n')
cat('# of Discordant region = ', sum(d),'\n')

# hist of log2_minus_over_plus
pdf(paste('Concordant','dummy',dummy,'minus_over_plus_hist.pdf', sep='_'))
#hist(Strand_count$log2_minus_over_plus[sub])
hist(Strand_count$log2_minus_over_plus[Strand_count$Con_Dis=='Con'])
dev.off()

# scatter plot 
pdf(paste(input_f,'dummy',dummy,'scatterPlot.pdf', sep='_'))
with(Strand_count, plot(log_mat_over_pat_plus, log_mat_over_pat_minus, 
                             col=ifelse(Con_Dis == 'Con', Con_col_1 , Dis_col_1 ),
                             xlim=c(-6,6),ylim=c(-6,6),
                             #xlim=c(-600,600), ylim=c(-600,600), 
                             pch=20,frame.plot=FALSE, main = paste('dummy=',dummy)))
abline(v=0)
abline(h=0)
abline(a=0, b=2^(1.5), col=Con_col_1,lty=2)
abline(a=0, b=1/2^(1.5), col=Con_col_1, lty=2)
# color the filter out points
sub = (Strand_count$log2_minus_over_plus < -1.5) | (Strand_count$log2_minus_over_plus > 1.5) & (Strand_count$Con_Dis == 'Con')
with(Strand_count, points(log_mat_over_pat_plus[sub], log_mat_over_pat_minus[sub],  
                               col= filter_out_col, pch=20))
legend("topleft", legend = c(paste("Concordant (n=",sum(c)-f,")",sep=""), paste("Discordant  (n=",sum(d),")",sep="")), col = c(Con_col_1, Dis_col_1),
       pch=20)
dev.off()
# find miss label Concordant or Discordant
#sub = (Strand_count$log_mat_over_pat_plus *Strand_count$log_mat_over_pat_minus < 0) & (Strand_count$Con_Dis == 'Con')
#Strand_count[sub,1:original_col_number]



# export the input file with region pass filter
new = Strand_count[ (Strand_count$Con_Dis == 'Dis') | (Strand_count$log2_minus_over_plus >= -1.5) & (Strand_count$log2_minus_over_plus <= 1.5),1:original_col_number]
new<- new[order(new$chrm, new$chrmStart),] 
colnames(new)[1] = paste('#',colnames(new)[1],sep='')

#View(new)
write.table(new,paste(substr(input_f, 1, nchar(input_f)-4),'_ConFiltered.bed',sep=''), row.names=FALSE, sep="\t", quote = FALSE)

