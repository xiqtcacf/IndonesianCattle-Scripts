##run_loter.py prefix_ref1 prefix_ref2 prefix_mix prefix_output 30
import sys
import allel
import numpy as np

print(sys.argv[1])
print(sys.argv[2])
print(sys.argv[3])
print(sys.argv[4])
print(int(sys.argv[5]))
r1=sys.argv[1]
r2=sys.argv[2]
mix=sys.argv[3]
output=sys.argv[4]
threads=int(sys.argv[5])
chrs=["NC_037328.1","NC_037329.1","NC_037330.1","NC_037331.1","NC_037332.1","NC_037333.1",
"NC_037334.1","NC_037335.1","NC_037336.1","NC_037337.1","NC_037338.1","NC_037339.1",
"NC_037340.1","NC_037341.1","NC_037342.1","NC_037343.1","NC_037344.1","NC_037345.1",
"NC_037346.1","NC_037347.1","NC_037348.1","NC_037349.1","NC_037350.1","NC_037351.1",
"NC_037352.1","NC_037353.1","NC_037354.1","NC_037355.1","NC_037356.1"]

workdir='/home/wlk579/2.0Indonesia_Bos_project/Revision_NC/LoterNoAustralia'

#print(r1)
#print(threads)

def vcf2npy(vcfpath):
    callset = allel.read_vcf(vcfpath)
    haplotypes_1 = callset['calldata/GT'][:,:,0]
    haplotypes_2 = callset['calldata/GT'][:,:,1]
    
    m, n = haplotypes_1.shape
    mat_haplo = np.empty((2*n, m))
    mat_haplo[::2] = haplotypes_1.T
    mat_haplo[1::2] = haplotypes_2.T
    
    return mat_haplo.astype(np.uint8)
    

import os
import loter.locanc.local_ancestry as lc
import matplotlib.pyplot as plt
from joblib import Parallel, delayed

def process(ch):
  popR1 = vcf2npy(os.path.join(workdir,'data',"by_chr",r1+"_"+ch+".vcf.gz" ))
  popR2 = vcf2npy(os.path.join(workdir,'data',"by_chr",r2+"_"+ch+".vcf.gz" ))
  popMix = vcf2npy(os.path.join(workdir,'data',"by_chr",mix+"_"+ch+".vcf.gz" ))
  res_loter = lc.loter_smooth(l_H=[popR1, popR2], h_adm=popMix, num_threads=threads) ## set the number of threads
  np.savetxt(os.path.join(workdir,"results","by_chr",output+"_"+ch+"_anc.txt"), res_loter, fmt="%i")

Parallel(n_jobs=25)(delayed(process)(ch) for ch in chrs)
