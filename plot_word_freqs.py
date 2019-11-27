#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 14:31:53 2019

Plot "word" freqs for Finnish haplotypes

@author: nbaya
"""

import pandas as pd
import matplotlib.pyplot as plt
#import operator
import ast
import math
import numpy as np
import scipy.stats as stats
from itertools import chain
finnhaps_wd = '/Users/nbaya/Documents/lab/finnhaps/'

cm_w_ls = [0.01, 0.02, 0.04, 0.08, 0.16]
chr_ls = [20,21,22]
ss = pd.DataFrame(data={'chr':list(chain(*[[x]*len(cm_w_ls) for x in chr_ls])),'cm_w':cm_w_ls*len(chr_ls), 'n_w': None, 'mean':None, 'median':None, #sumstats dataframe 
                        'std':None, 'iqr':None, 'min':None, 'max':None,
                        'snp_ct_std':None, 'max_snp_idx':None, 'prop_of_chr':None}) #, 'mode':None, })

total_cm_dict = {22: 73.93, 21:64.996, 20:111.7965}

hd_collapse = 1 # hamming distance to collapse by

def collapse(sorted_ls, hd=1):
    if hd is None or hd==0:
        return sorted_ls
    collapsed = []
    while len(sorted_ls)>0:
        top_ct, top_word = sorted_ls.pop(0)
    #    top_word_b = bin(int(word,2))
        for entry in sorted_ls:
            ct, word = entry
            if hamming(top_word, word)<=hd:
                sorted_ls.remove(entry)
                top_ct += ct
        collapsed += [(top_ct, top_word)]
    return collapsed

for chr in range(20,23):
    print(f'starting chr {chr}')
    info = 0.9
    maf = 0.01
    for cm_w in [0.01]:#cm_w_ls:
        n_w= int(total_cm_dict[chr]/cm_w) 
        
        with open(f'{finnhaps_wd}data/out.chr{chr}.n_w_{n_w}.maf_{maf}.cm_w_{cm_w}.info_{info}.txt', 'r') as f:
            lines = f.readlines()
        f.close()
        
        snp_idx_ls = [line.split('\t')[1].split(':')[1].split(',') for line in lines] # extract the SNP indices from each line
        maf_ls = [line.split('\t')[2].split(':')[1].split(',') for line in lines] # extract the MAFs from each line
        n_eff_ls = [float(line.split('\t')[3].split(':')[1].strip()) for line in lines] # extract the effective number of SNPs from each line
        w_idx_ls = [line.split('\t')[0].split(' ')[1] for line in lines] # extract the window index for each line
        sorted_ls = [ast.literal_eval(line.split('\t')[4].split(':')[1]) if len(line.split('\t'))>4 else None for line in lines ] # extract the "sorted_dict" from each line
        
        sorted_ls = [x for x in sorted_ls if x is not None]

        sorted_ls = [collapse(sorted_ls=x, hd=hd_collapse) for x in sorted_ls]
        
        n_eff_ls = []
        
        for sorted_dict in sorted_ls:
            cts = [x[0] for x in sorted_dict]
            freq = [x/sum(cts) for x in cts]
            n_eff = 1/np.linalg.norm(freq)**2
            n_eff_ls += [n_eff]
        
        max_words = max([len(x) for x in sorted_ls if x != None]) # max number of words across all windows
        max_ct = max([max([y[0] for y in x]) for x in sorted_ls if x != None]) # max word count across al windows
        
        
        ss.loc[(ss.chr==chr)&(ss.cm_w==cm_w),'n_w'] = len(snp_idx_ls)    
        ## calculate sumstats
        ss.loc[(ss.chr==chr)&(ss.cm_w==cm_w),'mean'] = np.mean(n_eff_ls)
        ss.loc[(ss.chr==chr)&(ss.cm_w==cm_w),'min'] = np.min(n_eff_ls)
        ss.loc[(ss.chr==chr)&(ss.cm_w==cm_w),'max'] = np.max(n_eff_ls)
        ss.loc[(ss.chr==chr)&(ss.cm_w==cm_w),'median'] = np.median(n_eff_ls)
        pdf = stats.kde.gaussian_kde(n_eff_ls)
        ss.loc[(ss.chr==chr)&(ss.cm_w==cm_w),'mode'] = n_eff_ls[np.argmax(pdf(n_eff_ls))]
        ss.loc[(ss.chr==chr)&(ss.cm_w==cm_w),'std'] = np.std(n_eff_ls)
        ss.loc[(ss.chr==chr)&(ss.cm_w==cm_w),'iqr'] = stats.iqr(n_eff_ls)
        ss.loc[(ss.chr==chr)&(ss.cm_w==cm_w),'snp_ct_mean'] = np.mean([len(x) for x in maf_ls])
        ss.loc[(ss.chr==chr)&(ss.cm_w==cm_w),'snp_ct_std'] = np.std([len(x) for x in maf_ls])
        ss.loc[(ss.chr==chr)&(ss.cm_w==cm_w),'max_snp_idx'] = snp_idx_ls[-1][-1]
        ss.loc[(ss.chr==chr)&(ss.cm_w==cm_w),'prop_of_chr'] = len(n_eff_ls)/n_w
        
#        print(f'cM: {cm_w}, expected windows: {n_w}, observed windows: {len(n_eff_ls)}')
        
        sorted_ls = [x for x in sorted_ls if x is not None ]
        print(f'{cm_w} {len(sorted_ls)}')
        
        ## plot distribution of n_eff
        max_n_eff = math.ceil(max(n_eff_ls))
        bins= 50 #np.linspace(1,max_n_eff,max_n_eff*2+1)
        fig, ax = plt.subplots(figsize=(6*1.2, 4*1.2))
        n_eff_ls_trunc = [n_eff for n_eff in n_eff_ls if n_eff <= max_n_eff]
        ax.hist(n_eff_ls_trunc, bins=bins)
        plt.title(f'Distribution of $n_e$ for chr {chr}, window size of {cm_w} cM\n(min: {round(min(n_eff_ls),4)}, max: {round(max(n_eff_ls),4)}, n windows={len(n_eff_ls)}{f", hd={hd_collapse}" if hd_collapse is not None else ""})')
        plt.xlabel('$n_e$',fontsize=15)
        plt.xlim(left=1)
        plt.savefig(f'{finnhaps_wd}plots/n_eff_hist.chr{chr}.{f"hd_{hd_collapse}." if hd_collapse is not None else ""}cm_w_{cm_w}.n_w_{len(n_eff_ls)}.maf_{maf}.info_{info}.png',dpi=300)
    
    for s in ['mean','median','std','iqr','n_w','snp_ct_mean','snp_ct_std','max_snp_idx','prop_of_chr']:#,'mode']:
        fig,ax=plt.subplots()
        ax.plot(ss.cm_w, ss[s],'.-')
        plt.title(f"{s.upper() if s is 'iqr' else s} of $n_e$ for chr{chr}{f', hd={hd_collapse}' if hd_collapse is not None else ''}")
        plt.xlabel('window length (cM)')
        plt.ylabel(s)
        plt.savefig(f'{finnhaps_wd}plots/n_eff.{s}.chr{chr}.{f"hd_{hd_collapse}." if hd_collapse is not None else ""}maf_{maf}.info_{info}.png',dpi=300)





plt_ct = 0
## plot word freq
for idx, w_sorted in enumerate(sorted_ls):
        if w_sorted ==None:
            continue
        fig, ax = plt.subplots(figsize=(6*2,4*2))
        w_sorted_r=w_sorted[::-1] #reverse order for plotting
        ax.barh(y=list(range(len(w_sorted))), height=1, align='center', 
                width=[x[0]/N for x in w_sorted_r], tick_label=[''.join(x[1]) for x in w_sorted_r])
        plt.yticks(fontsize=8)
        plt.xscale('log')
        plt.xlim([(1/N)/(10**0.05),1])
    #    plt.ylim([-.5,max_words+0.5])
        plt.xlabel('Word frequency')
        plt.title(f'chr {chr}, cM window: [{cm_w*int(w_idx_ls[idx])}, {cm_w*(int(w_idx_ls[idx])+1)}]\n(window size= {cm_w} cM, MAF>{maf}, INFO>{info})')
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



            

plt.bar(x = range(5), height=[122140, 117318, 114020, 104320, 71751],
        tick_label=['info≥'+str(x) for x in [0, 0.8, 0.9, 0.95, 0.99]])
plt.ylabel('SNPs')


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


%%timeit
a = 1
print(a)


'010101'


sorted_ls = [(15152, '0010000100000010000000010000'), (14060, '0010000110000100001100000000'), 
             (6127, '0000000000000000000000000000'), (5529, '0010011100010000000010000001'), 
             (5300, '0011011100000000000011000000'), (2411, '0010100110000100001100000000'), 
             (2381, '0010000100001010000000010000'), (2218, '0010000100000100001100000000'), 
             (1356, '0110001101100001000010101010'), (874, '0000000000000000000010000001'), 
             (831, '1010000100000010000000010100'), (755, '0110001101100001010010101010'), 
             (651, '0010000110000100101100000000'), (458, '0010011100010000000010000000'), 
             (241, '0010000100000000000000000000'), (211, '0000000000000000000000010100'), (39, '0010011100000000000010000000'), (28, '0010001100000000000010000000'), (22, '0010000100000000000011000000'), (8, '0010011100000010000010000000'), (8, '0010000100000010000000010100'), (7, '0010000100000010000000000000'), (6, '0010000100000000000000010000'), (5, '0010000100000010000000010001'), (5, '0010000000000010000000010000'), (3, '0010011000010000000010000001'), (3, '0010001100000010000000010000'), (3, '0010000100000110000000010000'), (2, '0011011100000000000001000000'), (2, '0010011100010110000000010000'), (2, '0010000110000100001100010000'), (2, '0010000100000100000000000000'), (2, '0010000100000000000010000001'), (1, '0110001101100001000010000000'), (1, '0011011100000000000010000000'), (1, '0010011100010000000011000001'), (1, '0010011100000000000000010000'), (1, '0010001100000000000010001010'), (1, '0010000100000000000010000000'), (1, '0010000100000000000001000000'), (1, '0010000100000000000000101010'), (1, '0010000100000000000000010100'), (1, '0000000110000100101100000000'), (1, '0000000000000000000010000000'), (1, '0000000000000000000000010000')]

## collapse by hamming distance
def collapse(sorted_ls, hd=1):
    collapsed_ls = []
    for window in sorted_ls:
        word_set = set([x[1] for x in window])
        new_window = []
        top = window[1] #most common word
        top_ct = window[0]
        
        



def hamming(a, b, ignore_ends=True):
    r'''
    From: https://colelyman.com/2018/09/19/hamming-distance-one-liner-in-python/
    '''
    return len([i for i in filter(lambda x: x[0] != x[1], 
                                  zip(a[1:-1] if ignore_ends else a, 
                                      b[1:-1] if ignore_ends else b))])
    
    
    
    
