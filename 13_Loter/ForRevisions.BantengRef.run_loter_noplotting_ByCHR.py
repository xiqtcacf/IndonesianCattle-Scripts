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
chrs=["1","2","3","4","5","6",
"7","8","9","10","11","12",
"13","14","15","16","17","18",
"19","20","21","22","23","24",
"25","26","27","28","29"]

workdir='/home/wlk579/2.0Indonesia_Bos_project/Revision_NC/BantengRefMap/Loter_BantengRef'

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
