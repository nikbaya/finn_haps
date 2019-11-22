from datetime import datetime as dt
from sys import getsizeof
import numpy as np
import pandas as pd
import argparse
from nltk import word_tokenize
from sklearn.feature_extraction.text import CountVectorizer
import gzip
import nltk
nltk.download('punkt')


parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--w',type=int, required=True, help='width of window in SNPs')
parser.add_argument('--n_w',type=int, required=False, default=1, help='number of windows to use')
parser.add_argument('--maf',type=float, required=False, default=0.01, help='MAF to filter SNPs')

args = parser.parse_args()

w = args.w
n_w = args.n_w
maf = args.maf

chr = 22

fname=f'/fs/projects/finn_haps/haps/FINRISK_R1_chr{chr}.haps'

genmap = pd.read_table(f'/fs/projects/finn_haps/maps/genetic_map_hg38_chr{chr}.txt',
			delim_whitespace=True)
genmap = genmap.rename(columns={'COMBINED_rate(cM/Mb)':'rate', 'Genetic_Map(cM)':'cm'})

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

def hap_patterns(w, n_w=1,maf=0.01, cm_w=0.01, info=0.9, snpstats=None):
	r"""
	Look at patterns in SNPs (passing MAF of `maf`) of window width `w`, for `n_w` windows.
	"""
	ngram_vec = CountVectorizer(analyzer='word', tokenizer=word_tokenize, ngram_range=(1, 1), min_df=1)
	fname_out = f'/homes/nbaya/finn_haps/haps/out.w_{w}.n_w_{n_w}.maf_{maf}.cm_w_{cm_w}.info_{info}.txt'
	snpstats = snpstats[(snpstats.maf>maf)&(snpstats['info']>info)]
	snpstats_dict = dict(zip(snpstats.position.tolist(), snpstats.maf.tolist()))
	missing_stats = 0
	with open(fname,'r') as f:
		snp_idx=-1 #negative 1 so that first SNP is indexed by 0
		consec_empty_ct = 0
		window_str_ls = []
		cm = 0 #current position in cm
		cm_left = 0 #starting cm on left side of window (i.e. the start of the window)
#                cm_prev = 0 #cm of previous snp
		for w_idx in range(n_w):
			lines = []
			maf_ls = []
			snp_idx_ls = []
			n_snps = 0
			outside = False #whether our snp idx is outside the current window
			while cm-cm_left<=cm_w:
				if n_snps == w and cm == cm_left: #if there are already n_w snps in the window and cm hasn't changed since start of the window
					print(f"cM hasn't changed for {n_snps} SNPs, starting new window")
					break
				if not outside:
					full_line = f.readline().split()
					snp_idx += 1
					metadata = full_line[:6]
					line = full_line[6:]
					bp = int(metadata[3])
					if bp < genmap.position.min():
						cm = 0
					else:
						cm = genmap[genmap.position<=bp]['cm'].max()
					if cm==0:
						continue
					elif cm-cm_left>cm_w:
						outside = True
						
					if len(line)==0:
						consec_empty_ct += 1
						if consec_empty_ct>10:
							print(f'Too many consecutively empty lines')
						break
					continue	
				consec_empty_ct = 0  
#				line_np = np.asarray(line,dtype=int)
#				alt_af = line_np.mean()
#				snp_maf = min(alt_af, 1-alt_af)
#				print(snp_maf)
				stats_tmp = snpstats[(snpstats.position==bp)]
				if len(stats_tmp)==0:
					#print(f'SNP {bp} missing stats!')
					missing_stats += 1
					continue
					#line_np = np.asarray(line,dtype=int)
					#alt_af = line_np.mean()
					#snp_maf = min(alt_af, 1-alt_af)
				snp_maf = stats_tmp.maf.values[0]
				print(f'SNP {bp}  MAF:{snp_maf} cm:{cm}') # NOTE: SNP is 1-indexed
				lines += [line]
				maf_ls += [snp_maf]
				snp_idx_ls += [snp_idx]
				n_snps+=1
			cm_left = cm
#			lines_np = np.asarray(lines)
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
			print(f'window #{w_idx}, len={len(sorted_dict)}, n_eff={n_eff}\n{sorted_dict}')
			
			window_str = f'w_idx:{w_idx}\tsnps:{",".join(str(idx) for idx in snp_idx_ls)}\t'
			window_str += f'maf:{",".join(str(x) for x in maf_ls)}\tsorted_dict:{sorted_dict}\t'
			window_str += f'n_eff:{n_eff}\n'
#			window_str_ls += [window_str]
			
			with open(fname_out, 'a') as f_out:
				f_out.write(window_str)
			f_out.close()
			
	f.close()
	print(f'SNPs missing stats: {missing_stats}')
#	window_str_all = '\n'.join(window_str_ls)
#	window_str_all + '\n'
#	with open(fname_out, 'w') as f_out:
#		f_out.write(window_str_all)
#	f_out.close()
if __name__=='__main__':
	start = dt.now()
	print(f'...starting to read hap patterns...')
	hap_patterns(w=w, n_w=n_w, maf=maf, snpstats=snpstats)
	print(f'window ("word") length: {w}, number of SNP windows: {n_w}, maf: {maf}')
	print(f'elapsed: {round(((dt.now()-start).seconds)/60, 2)} min')
