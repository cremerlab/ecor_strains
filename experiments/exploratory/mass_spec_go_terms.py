#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import size.viz
cor, pal = size.viz.matplotlib_style()
data = pd.read_csv('../../data/mapped_relative_protein_abundances.csv')

glycolysis_enzymes = ['pgi', 'kduL', 'yggF',
                      'ybhA', 'gplX', 'pfkA', 'pfkB',
                      'fbaB', 'fbaA', 'tpiA',
                      'gapA', 'pgk', 'gpmA',  
                      'gpmM', 'eno', 'ppsA', 
                      'pykF', 'pykA']

glycolysis = data[data['gene_name'].isin(glycolysis_enzymes)]
glycolysis = glycolysis.groupby(['strain', 'carbon_source', 'growth_rate'])['mass_frac'].sum().reset_index()
glycolysis['sector'] = 'glycolysis'
glycolysis['discriminator'] = f"genes: {', '.join(glycolysis_enzymes)}"

carb_tport_GO = ['GO:0008643; ', 'GO:0006861; ', 
                 'GO:0008644']

carb_tport = data[data['GO_terms'].str.contains(carb_tport_GO[0]) |
                  data['GO_terms'].str.contains(carb_tport_GO[1]) |
                  data['GO_terms'].str.contains(carb_tport_GO[2])]
carb_tport = carb_tport.groupby(['strain', 'carbon_source', 'growth_rate'])['mass_frac'].sum().reset_index()
carb_tport['sector'] = 'carbohydrate transport'
carb_tport['discriminator'] = f"GO Terms: {''.join(carb_tport_GO)}"


tca = data[data['GO_terms'].str.contains('GO:0006099; ')]
tca = tca.groupby(['strain', 'carbon_source', 'growth_rate'])['mass_frac'].sum().reset_index()
tca['sector'] = 'TCA cycle'
tca['discriminator'] = 'GO Terms: GO:0006099'

txn_reg_GO = ['GO:0006355; ',  'GO:0032583; ',
              'GO:0045449; ', 'GO:0061019; ']

txn_reg = data[data['GO_terms'].str.contains(txn_reg_GO[0]) |
                  data['GO_terms'].str.contains(txn_reg_GO[1]) |
                  data['GO_terms'].str.contains(txn_reg_GO[2]) | 
                  data['GO_terms'].str.contains(txn_reg_GO[3])]
txn_reg = txn_reg.groupby(['strain', 'carbon_source', 'growth_rate'])['mass_frac'].sum().reset_index()
txn_reg['sector'] = 'transcriptional regulators'
txn_reg['discriminator'] = f"GO Terms: {''.join(txn_reg_GO)}"

merged = pd.concat([glycolysis, carb_tport, tca, txn_reg], sort=False)

_dfs = []
for g, d in merged.groupby(['sector', 'carbon_source']):
    NCM = d[d['strain']=='NCM3722']
    d['NCM_fold_change'] = d['mass_frac'] / NCM['mass_frac'].values[0]
    _dfs.append(d)

df = pd.concat(_dfs, sort=False)

df.to_csv('../../data/filtered_relative_mass_fracs.csv', index=False)
