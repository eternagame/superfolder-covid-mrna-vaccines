from DeepDeg import DegScore
from arnie.bpps import bpps
from arnie.free_energy import free_energy
from arnie.mfe import mfe
import numpy as np
import sys

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

def write_degscore(row):
    mdl = DegScore(row['sequence'])
    return mdl.degscore

def write_AUP(row):
    mat = bpps(row['sequence'],package='eternafold')
    return np.mean(1-np.sum(mat, axis=0))

def write_MFE_struct(row):
    return mfe(row['sequence'])

def write_dG_MFE(row):
    return free_energy(row['sequence'], constraint=row['_MFE_struct_vienna'].replace('.','x'))


if __name__=="__main__":

	df = pd.read_csv(sys.argv[1])

	df['CAI'] = df.apply(lambda row: write_cai(row), axis=1)
	df['AUP'] = df.apply(lambda row: write_AUP(row), axis=1)
	df['_MFE_struct_vienna'] = df.apply(lambda row: write_MFE_struct(row), axis=1)
	df['dG(MFE)'] = df.apply(lambda row: write_dG_MFE(row), axis=1)
	df['DegScore'] = df.apply(lambda row: write_degscore(row), axis=1)

	df.to_csv('further_sequences_processed.csv',index=False)
