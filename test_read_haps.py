from datetime import datetime as dt
from sys import getsizeof
import numpy as np
import pandas as pd
import pickle
import os
import argparse
from nltk import word_tokenize
from sklearn.feature_extraction.text import CountVectorizer
import gzip
import nltk
# nltk.download('punkt')


parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--w',type=int, required=False, default=10, help='width of window in SNPs')
parser.add_argument('--n_w',type=int, required=False, default=1, help='number of windows to use')
parser.add_argument('--maf',type=float, required=False, default=0.01, help='MAF to filter SNPs')
parser.add_argument('--info',type=float, required=False, default=0.9, help='info score to filter SNPs')
parser.add_argument('--cm_w',type=float, required=False, default=0.01, help='window width in cM')
parser.add_argument('--chr',type=int, required=False, default=22, help='which chromosome to use')
parser.add_argument('--print_only',action='store_true', default=False, help='whether to just print results and not write to file')
args = parser.parse_args()

w = args.w
n_w = args.n_w
maf = args.maf
info = args.info
cm_w = args.cm_w
chr = args.chr
print_only = args.print_only


fname=f'/fs/projects/finn_haps/haps/FINRISK_R1_chr{chr}.haps'
#fname=f'/homes/nbaya/finn_haps/haps/FINRISK_R1_chr{chr}.haps.head10'

snpstats = pd.read_table(f'/fs/projects/finn_haps/snp_stats/FINRISK_R1_chr{chr}.snp_stats',
			delim_whitespace=True)

def read_all_lines():
	with open('FINRISK_R1_chr22_nometa.haps','r') as f:
		line=f.readline()
		ct=1
		while line:
#			print(len(line.split()))
			line=f.readline()
			ct+=1
		print(ct)

def hap_patterns(w=10, maf=0.01, cm_w=0.01, info=0.9, n_w=None, snpstats=None):
	r"""
	Look at patterns in SNPs (passing MAF of `maf`) of window width `w`, for `n_w` windows.
	"""
	
	ngram_vec = CountVectorizer(analyzer='word', tokenizer=word_tokenize, ngram_range=(1, 1), min_df=1)
	fname_out = f'/homes/nbaya/finn_haps/haps/out.chr{chr}.w_{w}.n_w_{n_w}.maf_{maf}.cm_w_{cm_w}.info_{info}.txt'
	snpstats = snpstats[(snpstats.maf>maf)&(snpstats['info']>info)]
	snpstats_dict = dict(zip(snpstats.position.tolist(), snpstats.maf.tolist()))
	snp_set = set(snpstats.position.tolist())
	missing_stats = 0
	w_lim = get_windows(chr=chr, cm_w=cm_w, write=True)
	with open(fname,'r') as f:
		w_ct = 0
		snp_idx=-1 #negative 1 so that first SNP is indexed by 0
		consec_empty_ct = 0
		n_snps_final = 0 # total number of snps included in windows
		window_str_ls = []
		w_idx = 0
		start_bp = w_lim[w_idx] #starting limit of the current window in base pairs
		stop_bp = w_lim[w_idx+1] #stopping limit of the current window in base pairs
		#cm = 0 #current position in cm
		#cm_left = 0 #starting cm on left side of window (i.e. the start of the window)
#                cm_prev = 0 #cm of previous snp
		print(f'...starting to read hap patterns...')
		while w_ct < n_w:
			lines = []
			maf_ls = []
			snp_idx_ls = []
			n_snps = 0
			outside = False #whether our snp idx is outside the current window
			#while cm-cm_left<=cm_w:
			while w_idx < len(w_lim):	
				#if n_snps == w and cm == cm_left: #if there are already n_w snps in the window and cm hasn't changed since start of the window
				#	print(f"cM hasn't changed for {n_snps} SNPs, starting new window")
				#	break
				if not outside:
					full_line = f.readline().split()
					snp_idx += 1
					if len(full_line)==0:
						print(f'reached end of file')
						return None
						#consec_empty_ct += 1
						#if consec_empty_ct>10:
						#	print(f'Too many consecutively empty lines')
						#	return None
						#continue
					#consec_empty_ct = 0
					metadata = full_line[:6]
					line = full_line[6:]
					position = int(metadata[3])
					if position not in snp_set: # if SNP is not in snpstats (and therefore has no info score)
						missing_stats +=1
						continue
					# elif position < genmap.position.min():
					#	cm = 0
					#else:
					#	cm = genmap[genmap.position<=position]['cm'].max()
					#if cm==0:
					#	continue
					#if cm-cm_left>cm_w and n_snps>0: #if current SNP is outside the window
					#	outside = True
					#	print(f'current SNP ({position}) is outside the window ({cm}-{cm_left}>{cm_w}')
					#	break
				if position > stop_bp:
					print(f'snp {position} is to the right of window {w_idx}')
					w_idx += 1
					start_bp, stop_bp = w_lim[w_idx], w_lim[w_idx+1]
					outside = True
					break
				elif position < start_bp:
					print(f'snp {position} has not reached the start of the first window')
					continue
				outside = False
#				line_np = np.asarray(line,dtype=int)
#				alt_af = line_np.mean()
#				snp_maf = min(alt_af, 1-alt_af)
#				print(snp_maf)
				snp_maf = snpstats_dict[position]
				#line_np = np.asarray(line,dtype=int)
				#alt_af = line_np.mean()
				#snp_maf_check = min(alt_af, 1-alt_af)
				#print(f'{snp_maf} {snp_maf_check}')
				print(f'SNP {position}  MAF:{snp_maf} w_idx:{w_idx}') # NOTE: SNP is 1-indexed
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
#			print(getsizeof(lines_T))
#			lines_np_T = lines_np.T
#			lines_T = lines_np_T.tolist() 
#			pattern_dict = {x:lines_T.count(x) for x in lines_T}
			sorted_dict = sorted((val,key) for (key,val) in pattern_dict.items())[::-1]
			freq = cts/cts.sum()
			n_eff = 1/np.linalg.norm(freq)**2
			print(f'window #{w_ct}, len={len(sorted_dict)}, n_eff={n_eff}\n{sorted_dict}')
			
			window_str = f'w_idx {w_ct}\tsnps:{",".join(str(idx) for idx in snp_idx_ls)}\t'
			window_str += f'maf:{",".join(str(x) for x in maf_ls)}\tsorted_dict:{sorted_dict}\t'
			window_str += f'n_eff:{n_eff}\n'
#			window_str_ls += [window_str]
			w_ct += 1
			if not print_only:
				with open(fname_out, 'a') as f_out:
					f_out.write(window_str)
				f_out.close()
			
			n_snps_final += n_snps
	f.close()
	print(f'read {snp_idx+1} SNPs, of which {n_snps_final} were included in words')
	print(f'SNPs missing stats: {missing_stats}')
#	window_str_all = '\n'.join(window_str_ls)
#	window_str_all + '\n'
#	with open(fname_out, 'w') as f_out:
#		f_out.write(window_str_all)
#	f_out.close()

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
		if n_w%1==0: # if (max_cm-min_cm) is perfectly divisible by cm_w
			n_w_lim = int(n_w)+1 # number of window limits
		else:
			n_w_lim = int(n_w)+2 # number of window limits
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
	hap_patterns(w=w, n_w=n_w, maf=maf, snpstats=snpstats)
	print(f'Parameters: window ("word") length: {w}, number of SNP windows: {n_w}, maf: {maf}')
	print(f'elapsed: {round(((dt.now()-start).seconds)/60, 2)} min')























