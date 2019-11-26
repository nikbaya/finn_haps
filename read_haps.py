from datetime import datetime as dt
from sys import getsizeof
import numpy as np
import pandas as pd
import pickle
import os
import random
import math
import argparse
from nltk import word_tokenize
from sklearn.feature_extraction.text import CountVectorizer
import gzip
import nltk
# nltk.download('punkt')


parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--n_w',type=int, required=False, default=1, help='number of windows to use (use -1 to use all windows on chr)')
parser.add_argument('--maf',type=float, required=False, default=0.01, help='MAF to filter SNPs')
parser.add_argument('--info',type=float, required=False, default=0.9, help='info score to filter SNPs')
parser.add_argument('--cm_w',type=float, required=False, default=0.01, help='window width in cM')
parser.add_argument('--chr',type=int, required=False, default=22, help='which chromosome to use')
parser.add_argument('--hd',type=int, required=False, default=0, help='max Hamming distance to collapse by')
parser.add_argument('--sampling_frac',type=float, required=False, default=0, help='fraction of windows for which sorted_dict is written to file')
parser.add_argument('--print_only',action='store_true', default=False, help='whether to just print results and not write to file')
parser.add_argument('--overwrite',action='store_true', default=False, help='whether to overwrite previously written outfile')
args = parser.parse_args()

n_w = args.n_w
maf = args.maf
info = args.info
cm_w = args.cm_w
chr = args.chr
hd = args.hd
sampling_frac = args.sampling_frac
print_only = args.print_only
overwrite = args.overwrite


fname=f'/fs/projects/finn_haps/haps/FINRISK_R1_chr{chr}.haps'
#fname=f'/homes/nbaya/finn_haps/haps/FINRISK_R1_chr{chr}.haps.head10'
#fname=f'/homes/nbaya/finn_haps/haps/FINRISK_R1_chr{chr}.haps.tail1000'

snpstats = pd.read_table(f'/fs/projects/finn_haps/snp_stats/FINRISK_R1_chr{chr}.snp_stats',
			delim_whitespace=True)



def hap_patterns(maf=0.01, cm_w=0.01, info=0.9, n_w=1, snpstats=None, hd=0, sampling_frac=0,
		print_only=False, overwrite=False):
	r"""
	Look at patterns in SNPs (passing MAF of `maf`) of window width `w`, for `n_w` windows.
	Only approx. `sampling_frac`*`n_w` sorted_dicts will be written to file.
	"""	
	ngram_vec = CountVectorizer(analyzer='word', tokenizer=word_tokenize, ngram_range=(1, 1), min_df=1)
	snpstats = snpstats[(snpstats.maf>maf)&(snpstats['info']>info)]
	snpstats_dict = dict(zip(snpstats.position.tolist(), snpstats.maf.tolist()))
	snp_set = set(snpstats.position.tolist())
	missing_stats = 0
	w_lim = get_windows(chr=chr, cm_w=cm_w, write=True)
	if n_w is -1:
                print(f'... reading all windows of width {cm_w} for chr {chr} ...')
                n_w = len(w_lim)-1
	n_keep = int(sampling_frac*n_w)
	keep_sorted_dict = [1]*n_keep+[0]*(n_w-n_keep)
	random.shuffle(keep_sorted_dict)
	if not print_only:
		print(n_w)
		print(len(w_lim)-1)
		fname_out = f'/homes/nbaya/finn_haps/haps/out.chr{chr}.cm_w_{cm_w}.n_w_{"all" if n_w==len(w_lim)-1 else n_w}.maf_{maf}.info_{info}.hd_{hd}.txt'
		if os.path.isfile(fname_out):
			if overwrite:
				os.remove(fname_out)
			else:
				print(f'... appending results to existing file ...')
		print(f'... writing to {fname_out} ...')
	with open(fname,'r') as f:
		w_ct = 0
		snp_idx=-1 #negative 1 so that first SNP is indexed by 0
		consec_empty_ct = 0
		n_snps_final = 0 # total number of snps included in windows
		window_str_ls = []
		w_idx = 0
		start_bp = w_lim[w_idx] #starting limit of the current window in base pairs
		stop_bp = w_lim[w_idx+1] #stopping limit of the current window in base pairs
		print(f'... starting to read hap patterns ...')
		while w_ct < n_w:
			lines = []
			maf_ls = []
			snp_idx_ls = []
			n_snps = 0
			outside = False #whether our snp idx is outside the current window
			while w_idx < len(w_lim):	
				if not outside:
					full_line = f.readline().split()
					snp_idx += 1
					if len(full_line)==0:
						print(f'reached end of file')
						w_ct = n_w
						snp_idx -= 1
						break
						#return None
					metadata = full_line[:6]
					line = full_line[6:]
					position = int(metadata[3])
					if position not in snp_set: # if SNP is not in snpstats (and therefore has no info score)
						missing_stats +=1
						continue
				if position > stop_bp:
					#print(f'snp {position} is to the right of window {w_idx}')
					w_idx += 1
					start_bp, stop_bp = w_lim[w_idx], w_lim[w_idx+1]
					outside = True
					break
				elif position < start_bp:
					#print(f'snp {position} has not reached the start of the first window')
					continue
				outside = False
				snp_maf = snpstats_dict[position]
				# print(f'SNP {position}  MAF:{snp_maf} w_idx:{w_idx}') # NOTE: SNP is 1-indexed
				lines += [line]
				maf_ls += [snp_maf]
				snp_idx_ls += [snp_idx]
				n_snps+=1
			if len(lines)==0:
				continue
			lines_T = list(zip(*lines))
			lines_T_str = [''.join(x) for x in lines_T]
			fit = ngram_vec.fit_transform(lines_T_str)
			vocab = list(ngram_vec.get_feature_names())
			cts = fit.sum(axis=0).A1
			pattern_dict = dict(zip(vocab, cts))
			sorted_dict = sorted((val,key) for (key,val) in pattern_dict.items())[::-1]
			n_eff_ls = []
			for hd in range(0,hd+1):
				sorted_dict_tmp = sorted_dict.copy()
				sorted_dict_tmp = collapse(sorted_ls=sorted_dict_tmp, hd=hd)
				cts = [x[0] for x in sorted_dict_tmp]
				freq = [x/sum(cts) for x in cts]
				n_eff = 1/np.linalg.norm(freq)**2
				n_eff_ls += [n_eff]
			print(f'window #{w_ct} [{cm_w*w_idx} cM, {cm_w*(w_idx+1)}], len={len(sorted_dict)}, n_eff={n_eff_ls}') #\n{sorted_dict}')
			window_str = f'w_idx {w_idx}\tsnps:{",".join(str(idx) for idx in snp_idx_ls)}\t'
			#window_str += f'maf:{",".join(str(x) for x in maf_ls)}\t' #sorted_dict:{sorted_dict}\t'
			window_str += f'n_eff:{",".join(str(x) for x in n_eff_ls)}'
			window_str += f'\tsorted_dict:{sorted_dict}' if keep_sorted_dict[w_ct] else ''
			window_str += '\n'
			w_ct += 1
			if not print_only:
				with open(fname_out, 'a') as f_out:
					f_out.write(window_str)
				f_out.close()
			
			n_snps_final += n_snps
	f.close()
	print(f'read {snp_idx+1} SNPs, of which {n_snps_final} were included in words')
	print(f'SNPs missing stats: {missing_stats}')

def collapse(sorted_ls, hd=1):
	r'''
	Collapses sorted_ls, a list of tuples, by Hamming distance `hd`.
	Each tuple in sorted_ls has word count as the first element, then the word
	as the second element.
	'''
	if hd is None or hd==0:
		return sorted_ls
	collapsed = []
	while len(sorted_ls)>0:
		top_ct, top_word = sorted_ls.pop(0)
		for entry in sorted_ls:
			ct, word = entry
			if hamming(top_word, word)<=hd:
				sorted_ls.remove(entry)
				top_ct += ct
		collapsed += [(top_ct, top_word)]
	return collapsed

def hamming(a, b, ignore_ends=True):
    return len([i for i in filter(lambda x: x[0] != x[1], 
                                  zip(a[1:-1] if ignore_ends else a, 
                                      b[1:-1] if ignore_ends else b))])

def get_windows(chr, cm_w, write=False):
	r"""
	Get windows of width `cm_w` cM, returning a list of base pair positions delimiting windows.
	Assumes genmap is sorted from smallest to largest cM.
	"""
	fname = f'/homes/nbaya/finn_haps/haps/w_lim.chr{chr}.cm_w_{cm_w}.txt'
	if os.path.isfile(fname):
		print(f'... reading in w_lim from existing file {fname} ...')
		with open(fname, 'rb') as f:
			w_lim = pickle.load(f)
		f.close()
	else:
		print(f'... getting window limits for window width {cm_w} cM ...')
		genmap = pd.read_table(f'/fs/projects/finn_haps/maps/genetic_map_hg38_chr{chr}.txt',
                        delim_whitespace=True)
		genmap = genmap.rename(columns={'COMBINED_rate(cM/Mb)':'rate', 'Genetic_Map(cM)':'cm'})
		min_bp, min_cm = genmap[['position','cm']].head(1).values[0]
		max_bp, max_cm = genmap[['position','cm']].tail(1).values[0]
		genmap.drop_duplicates(subset='cm',keep='first',inplace=True)
		genmap.loc[genmap.cm==max_cm,'position'] = max_bp #this maximimizes the range of base pair positions, in case the last SNPs have duplicate values that have been dropped
		w_lim = [min_bp] #list of window delimiteres in terms of base pair position
		w_lim_cm = [min_cm] #same as w_lim, but in terms of cM
		n_w = (max_cm-min_cm)/cm_w
		n_w_lim = math.ceil(n_w) #number of window limits
		for w_idx in range(1, n_w_lim-1): #starting at 1 instead of 0 because we already have the first window limit, min_bp
			a_bp, a_cm = genmap[genmap.cm<=min_cm+w_idx*cm_w][['position','cm']].tail(1).values[0]
			b_bp, b_cm = genmap[genmap.cm>=min_cm+w_idx*cm_w][['position','cm']].head(1).values[0]
			w_lim += [a_bp+(b_bp-a_bp)*(min_cm+w_idx*cm_w-a_cm)/(b_cm-a_cm)] # keep as float
		w_lim += [max_bp] # [15295879.0, 15344263.00714182, 15359409.993022274, 15382720.99153013, 15418341.246185433, 15458239.490914984, 15492371.441934822, 15530185.25036252, 15603472.316208476, 15616423.257155754, ...] first 10 entries
		
		if write:
			with open(fname, 'wb') as f:
				pickle.dump(w_lim, f)
			f.close()
	return w_lim

if __name__=='__main__':
	start = dt.now()
	hap_patterns(cm_w = cm_w, n_w=n_w, maf=maf, sampling_frac=sampling_frac, snpstats=snpstats,
		     hd=hd, print_only=print_only, overwrite=overwrite, info=info)
	print(f'Parameters: window ("word") length: {cm_w} cM, number of SNP windows: {n_w}, maf: {maf}')
	print(f'elapsed: {round(((dt.now()-start).seconds)/60, 2)} min')























