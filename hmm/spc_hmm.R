# data
setwd("~/Box Sync/Danko_lab_work/mouse_AlleleDB/hmm")
data <- read.table("counts_plus_hmm.txt",header=T,sep = "\t")
data$p <- data$mat_allele_count/data$total_reads_count
data$modified_mat <- round(30 * data$p)
#data$mu <- n*data$p
#data$vars <- data$mu*(1-data$p)
data$state_n <- 0
data$state_n[data$state == 'M'] <- 1
data$state_n[data$state == 'S'] <- 2
data$state_n[data$state == 'P'] <- 3

n=30
p1=0.9
p2=0.5
p3=0.1

data$simulated_mat <- 0
data$simulated_mat[data$state == 'M']=rnorm(sum(data$state == 'M'), mean=n*p1, sd=sqrt(n*p1*(1-p1)))
data$simulated_mat[data$state == 'S']=rnorm(sum(data$state == 'S'), mean=n*p2, sd=sqrt(n*p2*(1-p2)))
data$simulated_mat[data$state == 'P']=rnorm(sum(data$state == 'P'), mean=n*p3, sd=sqrt(n*p3*(1-p3)))


values <- data$simulated_mat
states <- data$state_n

library(rqhmm)
# hmm structure

normal.test <- function(values, states) {
#  N = length(means)
#  stopifnot(length(means) == length(vars))

  # 
  N=3
  hmm <- new.qhmm(list(1, NULL), # 2 emissions to set state path
                  rbind(c(1,2,3),c(1,2,3),c(1,2,3)), 
                  rep("discrete", N),
                  list("normal", "normal", "normal"))

  #
  # initial state
  set.initial.probs.qhmm(hmm, c(0.25, 0.5, 0.25))

  # set transition parameters
  set.transition.params.qhmm(hmm, 1, c(0.45, 0.5, 0.05))
  set.transition.params.qhmm(hmm, 2, c(0.05, 0.9, 0.05))
  set.transition.params.qhmm(hmm, 3, c(r))

  # set emission parameters
  # 1st emission track (actual normal emissions)
  set.emission.params.qhmm(hmm, 1, c(n*p1, n*p1*(1-p1)), slot = 1)
  set.emission.params.qhmm(hmm, 2, c(n*p2, n*p2*(1-p2)), slot = 1)
  set.emission.params.qhmm(hmm, 3, c(n*p3, n*p3*(1-p3)), slot = 1)

 # for (i in 1:N) {
   # 2nd emission track forces state path
 #   pars = rep(0, N)
 #   pars[i] = 1
 #   set.emission.params.qhmm(hmm, i, pars, fixed = rep(T, N), slot = 2)
 # }

 # dset = rbind(values, states)
  
  result1 = em.qhmm(hmm, list(data$simulated_mat))
  result2 = em.qhmm(hmm, list(data$modified_mat))
    
  # get means
  mus = sapply(1:N, function(i) get.emission.params.qhmm(hmm, i)[1])
  vars = sapply(1:N, function(i) get.emission.params.qhmm(hmm, i)[2])
  
  return(list(result = result, mu = mus, var = vars))
}

#
#
normal.test(data$modified_pat, data$state_n)
normal.test(data$simulated_pat, data$state_n)


#####
dset = data$simulated_mat
path = viterbi.qhmm(hmm, dset )
plot(data$snppos[1:10000], path[1:10000], xlab="SNP", ylab="states", ty='l', ylim = c(0,3)) # figure 3.6 (except the gray overlays)


# posterior decoding
fw = forward.qhmm(hmm, dset)
bk = backward.qhmm(hmm, dset)

logPx = attributes(fw)$loglik

posterior = exp(fw + bk - logPx)

lines(data$snppos[1:10000], posterior[1,1:10000], xlab="SNP", ylab="P(mat)", ty='l', col='red') # figure 3.6 (except the gray overlays)
lines(data$snppos[1:10000], posterior[3,1:10000], xlab="SNP", ylab="P(pat)", ty='l', col='blue') # figure 3.6 (except the gray overlays)
