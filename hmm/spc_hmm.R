# data
setwd("~/Box Sync/Danko_lab_work/mouse_AlleleDB/hmm")
data <- read.table("counts_plus_chrm1.txt",header=T,sep = "\t")
data$p <- data$pat_allele_count/data$total_reads_count
data$modified_pat <- round(30/data$total_reads_count * data$pat_allele_count)
#data$mu <- n*data$p
#data$vars <- data$mu*(1-data$p)
data$state_n <- 0
data$state_n[data$state == 'S'] <- 1
data$state_n[data$state == 'M'] <- 2
data$state_n[data$state == 'P'] <- 3

n=30
p1=0.5
p2=0.1
p3=0.9

data$simulated_pat <- 0
data$simulated_pat[data$state == 'S']=rnorm(sum(data$state == 'S'), mean=n*p1, sd=sqrt(n*p1*(1-p1)))
data$simulated_pat[data$state == 'M']=rnorm(sum(data$state == 'M'), mean=n*p2, sd=sqrt(n*p2*(1-p2)))
data$simulated_pat[data$state == 'P']=rnorm(sum(data$state == 'P'), mean=n*p3, sd=sqrt(n*p3*(1-p3)))


values <- data$simulated_pat
states <- data$state_n

library(rqhmm)
# hmm structure

normal.test <- function(values, states) {
#  N = length(means)
#  stopifnot(length(means) == length(vars))

  # 
  N=3
  hmm <- new.qhmm(list(c(1, 1), NULL), # 2 emissions to set state path
                  rbind(c(1,2,3),c(1,2,3),c(1,2,3)), 
                  rep("discrete", N),
                  rep(list(c("normal", "discrete")), N))

  #
  # initial state
  set.initial.probs.qhmm(hmm, c(0.8, 0.1, 0.1))

  # set transition parameters
  set.transition.params.qhmm(hmm, 1, c(0.9, 0.05, 0.05))
  set.transition.params.qhmm(hmm, 2, c(0.45, 0.5, 0.05))
  set.transition.params.qhmm(hmm, 3, c(0.45, 0.05, 0.5))

  # set emission parameters
  # 1st emission track (actual normal emissions)
  set.emission.params.qhmm(hmm, 1, c(n*p1, n*p1*(1-p1)), slot = 1)
  set.emission.params.qhmm(hmm, 2, c(n*p2, n*p2*(1-p2)), slot = 1)
  set.emission.params.qhmm(hmm, 3, c(n*p3, n*p3*(1-p3)), slot = 1)

  for (i in 1:N) {
   # 2nd emission track forces state path
    pars = rep(0, N)
    pars[i] = 1
    set.emission.params.qhmm(hmm, i, pars, fixed = rep(T, N), slot = 2)
  }

  dset = rbind(values, states)
  result = em.qhmm(hmm, list(dset))
    
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
dset = rbind(data$modified_pat, data$state_n)
path = viterbi.qhmm(hmm, dset)

# posterior decoding
fw = forward.qhmm(hmm, dset)
bk = backward.qhmm(hmm, dset)

logPx = attributes(fw)$loglik

posterior = exp(fw + bk - logPx)

plot(1:300, posterior[1,], type='l', xlab="roll", ylab="P(fair)") # figure 3.6 (except the gray overlays)

