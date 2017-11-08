import numpy as np
from math import *
import scipy.stats
#from sys import argv
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import time

#f_int = "counts_plus_hmm.txt"
#f_int = "counts_minus_hmm.txt"
f_int = "counts_hmm.txt"
#data
mat  = np.loadtxt(f_int, dtype=int ,delimiter='\t', usecols=[2], skiprows=1)
total = np.loadtxt(f_int, dtype=int ,delimiter='\t', usecols=[4], skiprows=1)
state = np.loadtxt(f_int, dtype=str ,delimiter='\t', usecols=[5], skiprows=1)
n_state = np.full(len(state), int(3), dtype=int)
n_state[state=="M"] = 0
n_state[state=="S"] = 1
n_state[state=="P"] = 2

#structure of hmm
#transition
# 0,1,2 = M, S, P
t_mm, t_ms, t_mp = 0.8, 0.1, 0.1
t_sm, t_ss, t_sp = 0.1, 0.8, 0.1
t_pm, t_ps, t_pp = 0.1, 0.1, 0.8
#t_mm, t_ms, t_mp = 0.45, 0.5, 0.05
#t_sm, t_ss, t_sp = 0.05, 0.9, 0.05
#t_pm, t_ps, t_pp = 0.05, 0.5, 0.45

T = np.log(np.array( [[t_mm, t_ms, t_mp],[t_sm, t_ss, t_sp],[t_pm, t_ps, t_pp]]))
#p_m, p_s,p_p = 0.9, 0.5, 0.1
p_m, p_s,p_p = 0.8, 0.5, 0.2
p = [p_m, p_s,p_p]

#intial prob
I_s=0.5
I_m=0.25
I_p=0.25


# emmision
binomlogpmf_dic={}
def get_emission_log_prob(x,n,p):
    if (x,n,p) not in binomlogpmf_dic:
        binomlogpmf_dic[(x,n,p)] = binomlogpmf(x, n, p)
    return binomlogpmf_dic[(x,n,p)]


def binomlogpmf(k, n, p):
    return scipy.stats.binom.logpmf(k, n, p)

def get_max_argmax(l):
    m = max(l)
    a = l.index(m)
    return m, a

def viterbi (p, T, x=mat, n=total):
    v = np.full((3, len(x)), float('-inf'))
    b = np.full((3, len(x)), int(3), dtype=int)
    #Initialization
    v[0, 0] = log(I_m) + get_emission_log_prob(x[0],n[0],p[0])
    v[1, 0] = log(I_s) + get_emission_log_prob(x[0],n[0],p[1])
    v[2, 0] = log(I_p) + get_emission_log_prob(x[0],n[0],p[2])
    b[0 ,0] = 0
    b[1 ,0] = 1
    b[2 ,0] = 2
    #Iteration
    for i in xrange(1, len(x)):
        v[0, i], b[0, i] = get_max_argmax([v[0, i-1] + T[0,0],  v[1, i-1] + T[1,0], v[2, i-1] + T[2,0]])
        v[1, i], b[1, i] = get_max_argmax([v[0, i-1] + T[0,1],  v[1, i-1] + T[1,1], v[2, i-1] + T[2,1]])
        v[2, i], b[2, i] = get_max_argmax([v[0, i-1] + T[0,2],  v[1, i-1] + T[1,2], v[2, i-1] + T[2,2]])
        v[0, i] += get_emission_log_prob(x[i],n[i],p[0])
        v[1, i] += get_emission_log_prob(x[i],n[i],p[1])
        v[2, i] += get_emission_log_prob(x[i],n[i],p[2])
    #track back pointer
    viterbi_path_backward=[]
    i = len(x) - 1
    f = np.argmax(v[:,i])
    viterbi_path_backward.append(f)
    i = i-1
    while (i>= 0):
        f = b[f,i]
        viterbi_path_backward.append(f)
        i -= 1
    viterbi_path_forward = np.array(viterbi_path_backward[::-1])
    plt.plot(np.arange(0,500), n_state[0:500], color='b')
    plt.plot(np.arange(0,500), viterbi_path_forward[0:500], color='red')
    plt.savefig('viterbi_path_forward.pdf')
    plt.close()
    #plt.show()
    return viterbi_path_forward


def sumLogProb(a, b):
    # a and b are log probabilities
    # return log(exp(a)+exp(b))
    if b is None:
        return a
    elif a is None:
        return b
    elif a > b:
        return a + np.log1p(exp(b - a))
    else:
        return b + np.log1p(exp(a - b))



def forward_probability_calculation(x=mat, n=total, p=p, T=T):
    f_p_m = np.full((3, len(x)), float('-inf'))
    # Initialization:
    f_p_m[0, 0] = log(I_m) + get_emission_log_prob(x[0],n[0],p[0])
    f_p_m[1, 0] = log(I_s) + get_emission_log_prob(x[0],n[0],p[1])
    f_p_m[2, 0] = log(I_p) + get_emission_log_prob(x[0],n[0],p[2])
    #Iteration
    for i in xrange(1, len(x)):
        f_p_m[0, i] = get_emission_log_prob(x[i],n[i],p[0]) + sumLogProb(sumLogProb(f_p_m[0, i-1] + T[0,0], f_p_m[1, i-1] + T[1,0]), f_p_m[2, i-1] + T[2,0])
        f_p_m[1, i] = get_emission_log_prob(x[i],n[i],p[1]) + sumLogProb(sumLogProb(f_p_m[0, i-1] + T[0,1], f_p_m[1, i-1] + T[1,1]), f_p_m[2, i-1] + T[2,1])
        f_p_m[2, i] = get_emission_log_prob(x[i],n[i],p[2]) + sumLogProb(sumLogProb(f_p_m[0, i-1] + T[0,2], f_p_m[1, i-1] + T[1,2]), f_p_m[2, i-1] + T[2,2])
    #Final value: 
    p_Y_f= sumLogProb(sumLogProb(f_p_m[0, len(x)-1], f_p_m[1,len(x)-1]),f_p_m[2,len(x)-1])
    return f_p_m, p_Y_f



def backward_probability_calculation(x=mat, n=total, p=p, T=T):
    b_p_m = np.full((3, len(x)), float('-inf'))
    # Initialization:
    b_p_m[0, len(x)-1] = log(1)  #???
    b_p_m[1, len(x)-1] = log(1)
    b_p_m[2, len(x)-1] = log(1)
    #Iteration
    i = len(x) - 2
    while i >= 0:
        b_p_m[0, i] = sumLogProb(sumLogProb(T[0,0] + b_p_m[0, i+1] + get_emission_log_prob(x[i+1],n[i+1],p[0]), T[0,1] + b_p_m[1, i+1] + get_emission_log_prob(x[i+1],n[i+1],p[1])),T[0,2] + b_p_m[2, i+1] + get_emission_log_prob(x[i+1],n[i+1],p[2]))
        b_p_m[1, i] = sumLogProb(sumLogProb(T[1,0] + b_p_m[0, i+1] + get_emission_log_prob(x[i+1],n[i+1],p[0]), T[1,1] + b_p_m[1, i+1] + get_emission_log_prob(x[i+1],n[i+1],p[1])),T[1,2] + b_p_m[2, i+1] + get_emission_log_prob(x[i+1],n[i+1],p[2]))
        b_p_m[2, i] = sumLogProb(sumLogProb(T[2,0] + b_p_m[0, i+1] + get_emission_log_prob(x[i+1],n[i+1],p[0]), T[2,1] + b_p_m[1, i+1] + get_emission_log_prob(x[i+1],n[i+1],p[1])),T[2,2] + b_p_m[2, i+1] + get_emission_log_prob(x[i+1],n[i+1],p[2]))
        i -=1
    #Final value:
    p_Y_b = sumLogProb(sumLogProb(log(I_m) + b_p_m[0, 0] + get_emission_log_prob(x[0],n[0],p[0]), log(I_s) + b_p_m[1, 0] + get_emission_log_prob(x[0],n[0],p[1])),log(I_p)+ b_p_m[2, 0] + get_emission_log_prob(x[0],n[0],p[2]))
    return b_p_m, p_Y_b


def em_interate(T, p, x=mat, n=total):
    t = time.time()
    f_p_m, p_Y_f = forward_probability_calculation(x, n, p, T)
    print "forward: ", t- time.time()
    b_p_m, p_Y_b = backward_probability_calculation(x, n, p, T)
    print "backward: ", t- time.time()
    
    #local P(Y)
    p_Y_l = np.full((1, len(x)), float('-inf'))
    for i in xrange(len(x)):
        p_Y_l[0,i] = sumLogProb(sumLogProb(b_p_m[0,i]+f_p_m[0,i], b_p_m[1,i]+f_p_m[1, i]),b_p_m[2, i]+f_p_m[2, i])
    print "p_Y_l ", t- time.time()
    #A = [[None, None, None], [None, None, None], [None, None, None]]
    A = np.zeros((3,3))
    new_P = [None, None, None] #P_m, P_s, P_p 
    
    # can add multiple sequence
    for i in xrange(len(x)-1):
        for k in range(3):
            for l in range(3):
                A[k,l] = A[k,l] + exp(f_p_m[k, i] + T[k,l] + get_emission_log_prob(x[i+1],n[i+1],p[l]) + b_p_m[l, i+1] - p_Y_l[0,i])
    #A = np.array(A)
    print "A : ", t- time.time()
    new_T = np.zeros((3,3))
    for k in range(3):
        new_P[k] = np.sum(np.exp(f_p_m[k,] + b_p_m[k,]-p_Y_l) * x ) / np.sum(np.exp(f_p_m[k,] + b_p_m[k,]-p_Y_l) * n ) 
        new_T[k] = A[k]/A[k].sum()
    #new_p_Y_f = forward_probability_calculation(p= new_P, T=np.log(new_T))[1]
    print "secs: ", t- time.time()
    print new_T, new_P, p_Y_f
    return np.log(new_T), new_P, p_Y_f


#def em(T ,p , threshold = 0.01, max_iter=10):
new_T, new_P, p_Y_f = em_interate(T, p, x=mat, n=total)
p_Y_f_list = [p_Y_f]
#p_Y_f += 10
max_iter = 70
for i in xrange(max_iter):
    print i
    #if p_Y_f - p_Y_f_list[i] > threshold:
    new_T, new_P, p_Y_f = em_interate(new_T, new_P, x=mat, n=total)
    p_Y_f_list.append(p_Y_f)
#print 'Start u = %.2f, theta_h = %.2f, theta_l = %.2f , the final u = %f, theta_h = %f, and theta_l = %f' %(u, theta_h, theta_l, new_u, new_theta_h, new_theta_l)
#return new_T, new_P, p_Y_f_list

def make_em_plot(em_p_Y_f_list, i, t):
    plt.plot(xrange(len(em_p_Y_f_list)),em_p_Y_f_list)
    plt.xlabel('# of iteration')
    plt.ylabel('log likelihood')
    plt.title(t)
    plt.savefig('pset3_sc2457_q3b_'+str(i)+'_plot.pdf')
    plt.close()



f_int = "counts_plus_hmm.txt"
data = np.loadtxt(f_int, dtype=str ,delimiter='\t', usecols=range(0,6), skiprows=1)
chrom = np.loadtxt(f_int, dtype=int ,delimiter='\t', usecols=[0], skiprows=1)
snppos= np.loadtxt(f_int, dtype=int ,delimiter='\t', usecols=[1], skiprows=1)
mat  = np.loadtxt(f_int, dtype=int ,delimiter='\t', usecols=[2], skiprows=1)
total = np.loadtxt(f_int, dtype=int ,delimiter='\t', usecols=[4], skiprows=1)
state = np.loadtxt(f_int, dtype=str ,delimiter='\t', usecols=[5], skiprows=1)
n_state = np.full(len(state), int(3), dtype=int)
n_state[state=="M"] = 0
n_state[state=="S"] = 1
n_state[state=="P"] = 2

v_path=[]

for i in xrange(1,20):
    t_c = total[chrom == i]
    x_c = mat[chrom == i]
    snp_c = snppos[chrom == i]

snppos_dic={}
for i in xrange(1,20):
    snppos_dic[i]=[]
    snp_c = snppos[chrom == i]
    temp=[snp_c[0]]
    for l in xrange(1,len(snp_c)):
        if snp_c[l]-snp_c[l-1] > 1000:
            snppos_dic[i].append(temp)
            temp=[snp_c[l]]
        else:
            temp.append(snp_c[l])


            
        
        
        
    v_path += (list(viterbi (x=x, n=n, p=new_P, T=new_T)))

state_map = {0:'M', 1:'S', 2:'P'}
with open(f_int[0:-4]+'_out.txt', 'w') as out:
    out.write("\t".join(['chrm','snppos','mat_allele_count','pat_allele_count','total_reads_count','state', 'hmm_state', 'hmm_post_m', 'hmm_post_s','hmm_post_p']))
    out.write("\n")
    for i in xrange(len(v_path)):
        out.write("\t".join(list(data[i])+[state_map[v_path[i]]]))
        out.write("\n")


plt.plot(snppos[0:1000], n_state[0:1000], color='b')
plt.plot(snppos[0:1000], v_path[0:1000], color='red')
plt.scatter(snppos[0:1000], v_path[0:1000], s=10)
plt.ylim(-0.5, 2.5)
plt.xlabel('snppos')
plt.savefig('viterbi_path_forward_snppos_0-1000.pdf')
plt.close()





def run():
    new_T, new_P, p_Y_f_list = em(T ,p ,threshold = 0.01, max_iter=100)

    
def main():
    
#q1_a
    v_path = viterbi (x=mat, n=total, p=p, T=T)
    plt.plot(np.arange(0,500), n_state[0:500], color='b')
    plt.plot(np.arange(0,500), v_path[0:500], color='red')
    plt.savefig('viterbi_path_forward.pdf')
    
#q1_b
    f_p_m, p_Y_f = forward_probability_calculation(x=mat, n=total, p=p, T=T)
    b_p_m, p_Y_b = backward_probability_calculation(x=mat, n=total, p=p, T=T)
    print 'logP(Y) from forward probability calculation =' , p_Y_f
    print 'logP(Y) from backward probability calculation =' , p_Y_b
    #logP(Y) from forward probability calculation = -22998.7330551
    #logP(Y) from backward probability calculation = -22998.7330551
    
    #local P(Y)
    p_Y_l = np.full((1, len(x)), float('-inf'))
    for i in xrange(len(x)):
        p_Y_l[0,i] = sumLogProb(sumLogProb(b_p_m[0,i]+f_p_m[0,i],  b_p_m[1, i]+f_p_m[1, i]),b_p_m[2, i]+f_p_m[2, i])
    
    marginal_posterior_probability_of_state_0 = np.exp(f_p_m[0,] + b_p_m[0,] -p_Y_l)
    marginal_posterior_probability_of_state_1 = np.exp(f_p_m[1,] + b_p_m[1,] -p_Y_l)
    marginal_posterior_probability_of_state_2 = np.exp(f_p_m[2,] + b_p_m[2,] -p_Y_l)
    plt.plot(np.arange(0,1000), marginal_posterior_probability_of_state_1[0,0:1000], color = 'b')
    plt.plot(np.arange(0,1000), marginal_posterior_probability_of_state_0[0,0:1000], color = 'r')
    #for u, v in index_interval:
    #    plt.plot(np.arange(u-1, v), marginal_posterior_probability_of_state_h[u-1: v], color='r', linewidth=2.0)
    #plt.savefig('pset3_sc2457_q1b_plot.pdf')
    plt.show()
    plt.close()
    

    
#q2_a
    real_interval=np.array([(66,417),(468,509),(527,728),(946,1000)])-1
    Cb = len(real_interval)*2 -1
    Cs = 1000 -1 -Cb
    Seq_h =[]
    for u, v in real_interval:
        Seq_h.append(Seq[u:v+1])
    Seq_h = ('').join(Seq_h)
    d_hG = Seq_h.count('G') + Seq_h.count('C')
    d_hA = Seq_h.count('A') + Seq_h.count('T')
    d_lG = Seq.count('G') + Seq.count('C') - d_hG
    d_lA = Seq.count('A') + Seq.count('T') - d_hA
#q3_a
    em_u, em_theta_h, em_theta_l, em_p_Y_f_list = em()
    plt.plot(xrange(len(em_p_Y_f_list)),em_p_Y_f_list)
    plt.xlabel('# of iteration')
    plt.ylabel('log likelihood')
    plt.savefig('pset3_sc2457_q3a_plot.pdf')
    #plt.show()
    plt.close()
#q3_b
    make_em_plot(em(Seq=Seq, threshold = 0.01, u = 0.5, theta_h = 0.6, theta_l = 0.4)[-1], 1, 'Start u = 0.5, theta_h = 0.6, theta_l = 0.4')
    make_em_plot(em(Seq=Seq, threshold = 0.01, u = 0.6, theta_h = 0.6, theta_l = 0.4)[-1], 2, 'Start u = 0.6, theta_h = 0.6, theta_l = 0.4')
    make_em_plot(em(Seq=Seq, threshold = 0.01, u = 0.2, theta_h = 0.8, theta_l = 0.2)[-1], 3, 'Start u = 0.2, theta_h = 0.8, theta_l = 0.2')
    make_em_plot(em(Seq=Seq, threshold = 0.01, u = 0.000001, theta_h = 0.8, theta_l = 0.2)[-1], 4, 'Start u = 0.000001, theta_h = 0.8, theta_l = 0.2')
    make_em_plot(em(Seq=Seq, threshold = 0.01, u = 0.7, theta_h = 0.8, theta_l = 0.2)[-1], 5, 'Start u = 0.7, theta_h = 0.8, theta_l = 0.2')
    make_em_plot(em(Seq=Seq, threshold = 0.01, u = 0.7, theta_h = 0.6, theta_l = 0.4)[-1], 6, 'Start u = 0.7, theta_h = 0.6, theta_l = 0.4')
    make_em_plot(em(Seq=Seq, threshold = 0.01, u = 0.7, theta_h = 0.7, theta_l = 0.3)[-1], 7, 'Start u = 0.7, theta_h = 0.7, theta_l = 0.3')
    make_em_plot(em(Seq=Seq, threshold = 0.01, u = 0.5, theta_h = 0.55, theta_l = 0.45)[-1], 9, 'Start u = 0.5, theta_h = 0.55, theta_l = 0.45')


if __name__ == '__main__':
    run()
    