import numpy as np


#infp="../counts_plus_hmm_regions_t1e-05_interestingHets_IGV.bed"
#infp="../counts_minus_hmm_regions_t1e-05_interestingHets_IGV.bed"
infp="counts_both_hmm_regions_t1e-05_interestingHets_IGV.bed"
chrom = np.loadtxt(infp, dtype=str ,delimiter='\t', usecols=[0], skiprows=0) 
bed_regions=np.loadtxt(infp, dtype=int ,delimiter='\t', usecols=[1,2], skiprows=0)
read=np.loadtxt(infp, dtype=str ,delimiter='\t', usecols=[3], skiprows=0)

def bed_overlap(bed1, bed2):
    s1, e1 = bed1
    s2, e2 = bed2
    if (s2 < e1 and e1 <= e2) or (s2 <= s1 and s1 < e2):
        return True
    else: 
        s1, e1 = bed2
        s2, e2 = bed1
        return (s2 < e1 and e1 <= e2) or (s2 <= s1 and s1 < e2)



def bed_distance(bed1, bed2):
    s1, e1 = bed1
    s2, e2 = bed2
    if bed_overlap(bed1, bed2):
        return 0
    else:
        return min(abs(s1-e2), abs(s2-e1))




dis_list=[]
dis_list_str=[]
for cn in xrange(1,23):
    c = "chr"+str(cn)
    bed_c=bed_regions[chrom== c]
    p_M = np.full((len(bed_c), len(bed_c)), 1000)
    for i in xrange(len(bed_c)):
        for j in xrange(i+1,len(bed_c)):
            #p_M[i,j] = bed_distance(bed_c[i], bed_c[j])
            #p_M[j,i] = p_M[i,j]
            d = bed_distance(bed_c[i], bed_c[j])
            #if d<100: #100Mb
            dis_list.append(d)
            dis_list_str.append(str(d))


dis=np.array(dis_list)
#m_dis=np.array(dis_list)
#p_dis=np.array(dis_list)

with open(infp+"_pairwise_distance.txt", "w") as out:
    out.write("\n".join(dis_list_str))

#import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
#n, bins, patches = plt.hist(dis[(dis>= 1000) & (dis<=10000000)], 500, cumulative=False, facecolor='g', alpha=0.75)

#b=500
for b in [50,500,5000]:
    n, bins, patches = plt.hist(dis[(dis<50000)], b, cumulative=False, facecolor='g', alpha=0.75)
    plt.title('dis[(dis<50K)]')
    plt.xlabel('AlleleHMM regions Pairewise distance (base)')
    plt.ylabel('Frequency')
    plt.savefig(str(b)+'.pdf')
    plt.close()



for b in [50,500,5000]:
    n, bins, patches = plt.hist(m_dis[(m_dis<50000)], b, cumulative=False, facecolor='blue', alpha=0.75)
    plt.title('m_dis[(m_dis<50K)]')
    plt.xlabel('AlleleHMM regions Pairewise distance (Kb)')
    plt.ylabel('Frequency')
    plt.savefig(str(b)+'.pdf')
    plt.close() 



########
infp="counts_both_hmm_regions_t1e-05_interestingHets_IGV.bed"
chrom = np.loadtxt(infp, dtype=str ,delimiter='\t', usecols=[0], skiprows=0) 
bed_regions=np.loadtxt(infp, dtype=int ,delimiter='\t', usecols=[1,2], skiprows=0)
read=np.loadtxt(infp, dtype=str ,delimiter='\t', usecols=[3], skiprows=0)
strand=np.loadtxt(infp, dtype=str ,delimiter='\t', usecols=[5], skiprows=0)

#dis_list=[]
p_M_list=[]
smallest_d_per_row=[]
smallest_d_per_row_str=[]
for cn in xrange(1,23):
    c = "chr"+str(cn)
    bed_c=bed_regions[chrom== c]
    p_M = np.full((len(bed_c), len(bed_c)), 5000000)
    for i in xrange(len(bed_c)):
        for j in xrange(i+1,len(bed_c)):
            p_M[i,j] = bed_distance(bed_c[i], bed_c[j])
            p_M[j,i] = p_M[i,j]
        smallest_d_per_row.append(min(p_M[i]))
        smallest_d_per_row_str.append(str(min(p_M[i])))
    p_M_list.append(p_M)
            #d = bed_distance(bed_c[i], bed_c[j])
            #if d<500000: #500K
            #    dis_list.append(d)


with open(infp+"_distanceToNearestRegion.txt", "w") as out:
    out.write("\n".join(smallest_d_per_row_str))

# smallest distance
smallest_d_per_row = np.array(smallest_d_per_row)
for b in [50,500,5000]:
    n, bins, patches = plt.hist(smallest_d_per_row[(smallest_d_per_row<50000)], b, cumulative=False, facecolor='g', alpha=0.75)
    plt.title("smallest_d_per_row<50Kb")
    plt.xlabel('Distance to nearesr AlleleHMM regions(base)')
    plt.ylabel('Frequency')
    plt.savefig(str(b)+'.pdf')
    plt.close() 



#make plot for EACH chromosome
Uplimit=20000
i=1
for p_M in p_M_list:
    print p_M.shape
    p_M[p_M > Uplimit] = Uplimit
    plt.imshow(p_M, cmap='hot')#, interpolation='nearest')
    plt.colorbar()
    plt.title('chr'+str(i))
    plt.savefig('chr'+str(i)+".pdf")
    plt.close()
    i += 1

#make one plot for All chromosome
p_M = np.full((len(chrom), len(chrom)), 0)
for cn in xrange(1,23):
    c = "chr"+str(cn)
    d=np.where(chrom==c)[0]
    for i in d:
        for j in d:
            p_M[i,j] = bed_distance(bed_regions[i], bed_regions[j])
    p_M_list.append(p_M)
    
p_M[p_M > 1000000] = 1000000
plt.imshow(p_M)#, cmap='hot')
plt.colorbar()
plt.savefig("chrAll.pdf")
plt.close()







# use a distance to cut bed regions into groups
# if dis(bed1, bed2) > SPECIFIC_NUMBER, put into seperate group




