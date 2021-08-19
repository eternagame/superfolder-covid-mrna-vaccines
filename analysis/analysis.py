from arnie.bpps import bpps
from arnie.free_energy import free_energy
from arnie.mfe import mfe
import numpy as np
import pandas as pd
import gzip
import sys, os
from DegScore import DegScore

'''
Example usage: python analysis.py further_sequences.csv

Input: csv containing field called `sequence` that contains CDS sequences to be analyzed.
Writes csv containing metrics related to RNA stabilization against hydrolysis.

dG(MFE)
Average Unpaired Probability (AUP)
DegScore

Dependencies: Arnie (www.github.com/DasLab/arnie)
			  DeepDeg (www.github.com/Eterna/DeepDeg)

'''
from scipy.stats import gmean

C2W = {'UUU': 0.85185, 'UUC': 1.0, 'UUA': 0.2, 'UUG': 0.325, 'UCU': 0.79167,
	   'UCC': 0.91667, 'UCA': 0.625, 'UCG': 0.20833, 'UAU': 0.78571, 'UAC': 1.0,
	   'UAA': 0.6383, 'UAG': 0.51064, 'UGU': 0.85185, 'UGC': 1.0, 'UGA': 1.0,
	   'UGG': 1.0, 'CUU': 0.325, 'CUC': 0.5, 'CUA': 0.175, 'CUG': 1.0,
	   'CCU': 0.90625, 'CCC': 1.0, 'CCA': 0.875, 'CCG': 0.34375, 'CAU': 0.72414,
	   'CAC': 1.0, 'CAA': 0.36986, 'CAG': 1.0, 'CGU': 0.38095, 'CGC': 0.85714,
	   'CGA': 0.52381, 'CGG': 0.95238, 'AUU': 0.76596, 'AUC': 1.0, 'AUA': 0.3617,
	   'AUG': 1.0, 'ACU': 0.69444, 'ACC': 1.0, 'ACA': 0.77778, 'ACG': 0.30556,
	   'AAU': 0.88679, 'AAC': 1.0, 'AAA': 0.75439, 'AAG': 1.0, 'AGU': 0.625,
	   'AGC': 1.0, 'AGA': 1.0, 'AGG': 1.0, 'GUU': 0.3913, 'GUC': 0.52174,
	   'GUA': 0.26087, 'GUG': 1.0, 'GCU': 0.675, 'GCC': 1.0, 'GCA': 0.575,
	   'GCG': 0.275, 'GAU': 0.85185, 'GAC': 1.0, 'GAA': 0.72414, 'GAG': 1.0,
	   'GGU': 0.47059, 'GGC': 1.0, 'GGA': 0.73529, 'GGG': 0.73529}

def cai(seq):

	# Fixes sequence format. Convert to RNA
	seq = seq.upper().replace('T','U')
	if len(seq) % 3 != 0:
		raise ValueError("Not a valid coding sequence. Length is not a multiple of 3.")

	# Gets all the weights per codon
	w_list = []
	for i in range(0, len(seq), 3):
		codon = seq[i:i+3]
		# Do not count W or M codon since there is only one that encodes them
		if codon not in ['UGG', 'AUG']:
			w_list.append(C2W[codon])

	# return the geometric mean
	return gmean(w_list)

def write_degscore(row):
	mdl = DegScore(row['sequence'])
	struct = mdl.structure

	psu_mdl = DegScore(row['sequence'],structure=struct,mask_U=True)
	return mdl.est_half_life, psu_mdl.est_half_life

def write_ensemble_metrics(row, id_field):

	bp_filename = 'bpps_%s.npy.gz' % str(row[id_field])

	# bpps exist already
	if os.path.exists(bp_filename):
		f = gzip.GzipFile(bp_filename, "r")
		bpp_mat = np.load(f)
		print('loaded ', row[id_field])

	else: # calculate bpps and cache
		bpp_mat = bpps(row['sequence'], package='eternafold')
		print("wrote ", row[id_field])

		# cache
		f = gzip.GzipFile(bp_filename, 'w')
		np.save(file=f, arr=bpp_mat)
		f.close()

	punp_vector = 1-np.sum(bpp_mat, axis=0)

	aup = np.mean(punp_vector)
	sup_init = np.mean(punp_vector[:14])

	return aup, sup_init

def write_MFE_struct(row):
	return mfe(row['sequence'])

def write_dG_MFE(row):
	return free_energy(row['sequence'], constraint=row['_MFE_struct_vienna'].replace('.','x'))


if __name__=="__main__":

	df = pd.read_csv(sys.argv[1])
	print(df)

	df['CAI'] = df.apply(lambda row: cai(row['sequence']), axis=1)
	df[['AUP', 'AUP init 14']] = df.apply(lambda row: write_ensemble_metrics(row, 'Designer'), axis=1,result_type='expand')
	df['_MFE_struct_vienna'] = df.apply(lambda row: write_MFE_struct(row), axis=1)
	df['dG(MFE)'] = df.apply(lambda row: write_dG_MFE(row), axis=1)
	df['DegScore_half_life','DegScore_PSU_half_life'] = df.apply(lambda row: write_degscore(row), axis=1, result_type='expand')

	df.to_csv('%s_WITH_METRICS.csv' % sys.argv[1].replace('.csv',''),index=False)
