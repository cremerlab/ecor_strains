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


#%% Make a figure

fig, ax = plt.subplots(1, 2, figsize=(2,1), sharey=True)
ax[0].set_yticks([0, 1, 2, 3])
ax[0].set_yticklabels(['\ntranscription\nregulation', '\ncarb.\ntransport', '\nTCA cycle', '\nglycolysis'],
                      fontsize=5)

fig.text(0.2, -0.1, 'expression change relative to NCM3722', fontsize=5)

idx = {'glycolysis':3, 'TCA cycle':2, 'carbohydrate transport':1, 'transcriptional regulators':0}

nudges = {'AC1':0, 'ECOR02':0.2, 'ECOR51':0.4, 'ECOR63':0.6}
colors = {'AC1':cor['primary_blue'], 'ECOR02':cor['primary_green'], 'ECOR51':cor['primary_purple'], 'ECOR63':cor['primary_black']}
for k, v in colors.items():
    ax[1].plot([], [], 's', ms=4, alpha=0.75, color=v, label=k)
leg = ax[1].legend(fontsize=4, loc='upper right', handletextpad=0.05, bbox_to_anchor=(1.1, 1))

ax[0].set_title('acetate growth', fontsize=5)
ax[1].set_title('glucose growth', fontsize=5)
for g, d in df[df['strain']!='NCM3722'].groupby(['strain', 'sector', 'carbon_source']):
    if g[-1] == 'acetate':
        a = ax[0]
    else:
        a = ax[1]
    a.barh(idx[g[1]]-nudges[g[0]], d['NCM_fold_change'], height=0.2,

            color=colors[g[0]], alpha=0.75)
plt.savefig('./mass_spec_change.pdf', bbox_inches='tight')