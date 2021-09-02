import pandas as pd
import matplotlib.pyplot as plt
from RiboGraphViz import RGV
from DegScore import DegScore

color_vec=[]
struct_string = []
labels=[]

df_o = pd.read_csv('superfolders_ALL.csv')

df = df_o.loc[df_o.Variant=='S-2P, Delta'][df_o.Nucleotide=='Unmod']

for _, row in df.iterrows():
    deg_mdl = DegScore(row['sequence'], structure=row['_MFE_struct_vienna'])
    color_vec.extend(deg_mdl.degscore_by_position)
    struct_string.append(row['_MFE_struct_vienna'])
    labels.append(row['Designer'])

plt.figure(figsize=(50,50))
mdl = RGV(' '.join(struct_string))
mdl.draw(c=color_vec, cmap='plasma', struct_label=labels)
plt.savefig("Delta_unmod_designs.pdf", bbox_inches='tight')


df = df_o.loc[df_o.Variant=='S-2P, Delta'][df_o.Nucleotide=='PSU']
print(len(df))
for _, row in df.iterrows():
    deg_mdl = DegScore(row['sequence'], structure=row['_MFE_struct_vienna'],mask_U=True)
    color_vec.extend(deg_mdl.degscore_by_position)
    struct_string.append(row['_MFE_struct_vienna'])
    labels.append(row['Designer'])

plt.figure(figsize=(50,50))
mdl = RGV(' '.join(struct_string))
mdl.draw(c=color_vec, cmap='plasma', struct_label=labels)
plt.savefig("Delta_PSU_designs.pdf", bbox_inches='tight')
