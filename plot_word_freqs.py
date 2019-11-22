#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 14:31:53 2019

Plot "word" freqs for Finnish haplotypes

@author: nbaya
"""

#import pandas as pd
import matplotlib.pyplot as plt
#import operator
import ast
import numpy as np
import scipy.stats as stats

finnhaps_wd = '/Users/nbaya/Documents/lab/finnhaps/'

w=9
n_w=int(10000/w) #2000
maf = 0.01


with open(f'{finnhaps_wd}data/out.w_{w}.n_w_{n_w}.maf_{maf}.txt', 'r') as f:
    lines = f.readlines()
f.close()

sorted_ls = [ast.literal_eval(line.split('\t')[3].split(':')[1]) for line in lines] # extract the "sorted_dict" from each line
maf_ls = [line.split('\t')[2].split(':')[1].split(',') for line in lines] # extract the MAFs from each line
snp_idx_ls = [line.split('\t')[1].split(':')[1].split(',') for line in lines] # extract the SNP indices from each line
n_eff_ls = [float(line.split('\t')[4].split(':')[1].strip()) for line in lines] # extract the effective number of SNPs from each line

max_words = max([len(x) for x in sorted_ls]) # max number of words across all windows
max_ct = max([max([y[0] for y in x]) for x in sorted_ls]) # max word count across al windows

N = sum([x[0] for x in sorted_ls[0]]) # total count of haplotypes

plt_ct = 0
## plot word freq
for idx, w_sorted in enumerate(sorted_ls):
    fig, ax = plt.subplots(figsize=(6*2,4*2))
    w_sorted_r=w_sorted[::-1] #reverse order for plotting
    ax.barh(y=list(range(len(w_sorted))), height=1, align='center', 
            width=[x[0]/N for x in w_sorted_r], tick_label=[''.join(x[1]) for x in w_sorted_r])
    plt.yticks(fontsize=8)
    plt.xscale('log')
    plt.xlim([(1/N)/(10**0.05),1])
#    plt.ylim([-.5,max_words+0.5])
    plt.xlabel('Word frequency')
    plt.title(f'window #{idx} (window size= {w} SNPs, MAF>{maf})')
    plt.subplots_adjust(left=0.1)

    #plot maf
    ax2 = fig.add_axes([0.72, 0.21, 0.15, 0.2])
    ax2.bar(x=range(len(maf_ls[idx])),height=[float(x) for x in maf_ls[idx]],
            tick_label=[f'{ids}' for ids in snp_idx_ls[idx]])
    plt.xlabel('SNP idx')
    plt.ylabel('MAF')
    plt.xticks(rotation=45)
    ax2.axhline(y=maf,ls='--',c='k')
    
    plt_ct+=1
#    plt.savefig(f'{finnhaps_wd}plots/word_distr.w_{w}.window_{idx}.maf_{maf}.png',dpi=300)
    if plt_ct>1:
        plt.close()



## plot distribution of n_eff

fig, ax = plt.subplots()
ax.hist(n_eff_ls,50,density=True)
plt.title(f'Distribution of $n_e$\n(window size={w} SNPs, n windows={len(n_eff_ls)}, MAF>{maf})')
plt.xlabel('$n_e$',fontsize=15)
plt.xlim(left=1)
plt.savefig(f'{finnhaps_wd}plots/n_eff_hist.w_{w}.n_w_{len(n_eff_ls)}.maf_{maf}.png',dpi=300)
        


fit_alpha, fit_loc, fit_scale=stats.gamma.fit(n_eff_ls)
print(fit_alpha, fit_loc, fit_scale)

loc, scale = stats.expon.fit(n_eff_ls)

x = np.linspace(1, max(n_eff_ls),100)

fig, ax = plt.subplots()
ax.hist(n_eff_ls,50,density=True)
#ax.plot(x, stats.expon.pdf(x-loc)*scale)
ax.plot(x, stats.gamma.pdf(x, a=fit_alpha, loc=fit_loc, scale=fit_scale))
plt.title(f'Distribution of $n_e$\n(window size={w} SNPs, n windows={len(n_eff_ls)}, MAF>{maf})')
plt.xlabel('$n_e$',fontsize=15)
plt.xlim(left=1)



## collapse by hamming distance

def collapse(sorted_ls, hd=1):
    for window in sorted_ls:
        top = window[1] #most common word
        top_ct = window[0]
        


def hamming(a, b, ignore_ends=True):
    r'''
    From: https://colelyman.com/2018/09/19/hamming-distance-one-liner-in-python/
    '''
    return len([i for i in filter(lambda x: x[0] != x[1], 
                                  zip(a[1:-1] if ignore_ends else a, 
                                      b[1:-1] if ignore_ends else b))])
    
    
    
    