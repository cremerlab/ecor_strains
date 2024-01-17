#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import size.viz 
cor, pal = size.viz.matplotlib_style()

data = pd.read_csv('../../data/mapped_relative_protein_abundances.csv') 
rp_data = pd.read_csv('../../data/ecor_rp_measurements.csv')
rib_genes = ['rrsA', 'rpsA', 'rpsB', 'rpsC', 'rpsD', 'rpsE', 'rpsF', 'rpsG',
'rpsH', 'rpsI', 'rpsJ', 'rpsK', 'rpsL', 'rpsM', 'rpsN', 'rpsO', 'rpsP', 'rpsQ',
'rpsR', 'rpsS', 'rpsT', 'rpsU', 'rrlA', 'rrfA', 'rplA', 'rpsB', 'rplC', 'rplD',
'rplE', 'rplF', 'rplJ', 'RplL', 'rplI', 'rplK', 'rplM', 'rplN', 'rplO', 'rplP',
'rplQ', 'rplR', 'rplS', 'rplT', 'rplU', 'rplV', 'rplW', 'rplX', 'rplY', 'rpmA',
'rpmB', 'rpmC', 'rpmD', 'rplE', 'rpmF', 'rplmG', 'rpmH', 'rpmI','rpmJ', 'rrsA',
'rpsA', 'rpsB', 'rpsC', 'rpsD', 'rpsE', 'rpsF', 'rpsG', 'rpsH','rpsI', 'rpsJ',
'rpsK', 'rplsL', 'rpsM', 'rpsN', 'rpsO', 'rpsP', 'rpsQ', 'rpsR', 'rpsS', 'rpsT',
'rpsU']

rib = data[data['gene_name'].isin(rib_genes)].groupby(['strain', 'carbon_source', 'growth_rate']).sum().reset_index()

lit_data = pd.read_csv('../../data/Fig4_ecoli_ribosomal_mass_fractions.csv')


fig, ax = plt.subplots(1,1, figsize=(2,2))

ax.plot(lit_data['growth_rate_hr'], lit_data['mass_fraction'], 'o', ms=4, color=cor['light_black'],
        markeredgewidth=0, alpha=0.35, label='literature data')

for g, d in rp_data.groupby('strain'):
    ax.plot(d['growth_rate_hr'], d['mass_frac'], 'v', label=f'{g} - RNA/Protein', ms=4, alpha=0.65,
            markeredgewidth=0)

for g, d in rib.groupby('strain'):
    ax.plot(d['growth_rate'], d['mass_frac'], 'o', label=f'{g} - mass spec', ms=4,
            markeredgecolor='k', markeredgewidth=0.75)

ax.set_xlabel('growth rate [hr$^{-1}$]', fontsize=6)
ax.set_ylabel('ribosomal mass fraction', fontsize=6)
ax.legend(fontsize=4, bbox_to_anchor=(1,1))
plt.savefig('./rline_trends_mass_spec.pdf', bbox_inches='tight')
