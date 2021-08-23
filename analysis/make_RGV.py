from RiboGraphViz import RGV
import sys
import matplotlib.pyplot as plt
import pandas as pd
from DegScore import DegScore


def rgv_local(row):
    deg_mdl = DegScore(row['sequence'], structure = row['_MFE_struct_vienna'])
    deg_vec = deg_mdl.degscore_by_position
    
    deg_mdl_PSU = DegScore(row['sequence'], structure = row['_MFE_struct_vienna'],mask_U=True)
    deg_vec_PSU = deg_mdl_PSU.degscore_by_position
    
    rg = RGV(row['_MFE_struct_vienna'])
    plt.figure(figsize=(20,20))
    rg.draw(c=deg_vec)
    plt.savefig('%s.pdf' % row['Designer'],bbox_inches='tight')

    plt.figure(figsize=(20,20))
    rg.draw(c=deg_vec_PSU)
    plt.savefig('%s_PSU.pdf' % row['Designer'],bbox_inches='tight')

df = pd.read_csv(sys.argv[1])

df.apply(lambda row: rgv_local(row), axis=1)
