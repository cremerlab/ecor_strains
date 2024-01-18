#%%
import numpy as np 
import pandas as pd 

standard = pd.read_csv('../../../data/compiled_mass_fractions.csv')
standard = standard[(standard['strain']=='NCM3722') & 
                    (standard['growth_rate_hr'] >= 0.90) & 
                    (standard['growth_rate_hr'] <= 1.2)]
ref_cond = standard.groupby(['gene_name', 'b_number', 'cog_class', 'cog_letter', 'go_terms'])['mass_frac'].agg(('mean')).reset_index()
ref_cond
data = pd.read_csv('../../../data/relative_protein_abundances.csv')

#%% 
# data['gn_lower'] = data['Gene names'].str.lower()
data['ref_mass_frac'] = None
data['gene_name'] = None
data['b_nubmer'] = None
data['cog_letter'] = None
data['cog_class'] = None
names = data['Gene names'].str.lower()
cols = [k for k in data.keys() if ('log2' not in k) & ('Relative_concentration' in k)]
ref_col = data['Relative_concentration_NCM_Glu'].values
mapped = pd.DataFrame([])
err = 0

growth_rates = {'ECOR2': {'glucose': 0.8861, 'acetate': 0.2738, 'strain':'ECOR02'},
                'ECOR51': {'glucose': 0.3208, 'acetate':0.2756, 'strain':'ECOR51'},
                'ECOR63': {'glucose': 0.5048, 'acetate': 0.5526, 'strain':'ECOR63'},
                'AC': {'glucose':0.7299, 'acetate': 0.2021, 'strain':'AC1'},
                'NCM': {'glucose':0.9823, 'acetate':0.4069, 'strain':'NCM3722'}}
#%%
for g, d in ref_cond.groupby(['gene_name', 'b_number', 'cog_class', 'cog_letter', 'go_terms']):
    ind = None
    for i, n in enumerate(names):
            if g[0].lower() in n:
                ind = i
                break
            elif g[1].lower() in n:
                 ind = i 
                 break 
    
    if ind is not None:
        _df = data.iloc[ind]
        for c in cols:
            cond = c.split('_')[-1]
            if cond == 'Ace':
                cond = 'acetate'
            else:
                cond = 'glucose'
  
            df = pd.DataFrame({'gene_name': g[0],
                               'b_number': g[1],
                               'cog_letter': g[3],
                               'cog_class': g[2],
                               'GO_terms': d['go_terms'].values[0],
                               'relative_conc': _df[c] / ref_col[ind],
                               'ref_mass_frac': d['mass_frac'].values[0],
                               'strain': growth_rates[c.split('_')[2]]['strain'],
                               'growth_rate': growth_rates[c.split('_')[2]][cond],
                               'carbon_source': cond}, index =[0]) 
            mapped = pd.concat([mapped, df], sort=False)
    else:
        err += 1
#%%
mapped['mass_frac'] = mapped['relative_conc'] * mapped['ref_mass_frac']
mapped.to_csv('../../../data/mapped_relative_protein_abundances.csv', index=False)