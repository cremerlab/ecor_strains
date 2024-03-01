#%%
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import ecor.viz 
import ecor.model
const = ecor.model.load_constants()
cor, pal = ecor.viz.matplotlib_style()

# Load data and select interesting carbon source
data = pd.read_csv('../../data/ecor_rp_measurements.csv')

carb_source = data[data['carbon_source']=='Glucose']
carb_source = carb_source.groupby('strain')[['growth_rate_hr', 'rna_to_protein']].mean().reset_index()
carb_source = carb_source[carb_source['strain'].isin(['ECOR02', 'NCM3722', 'ECOR59', 'ECOR63'])]

# Evaluate the optimal allocation model
nu_range = np.linspace(0.05, 35, 300)



def estimate_nu(phiRb, lam):
    pref = lam / (1 - const['phi_O'] - phiRb)
    numer = const['Kd_cpc'] * lam
    denom = const['gamma_max'] * phiRb * (1 - lam/phiRb)
    return pref * (1 + numer/denom)
data['estimated_nu_inv_hr'] =  estimate_nu(data['rna_to_protein'].values * 0.4558,
                                    data['growth_rate_hr'])
#%%
opt_phiRb = ecor.model.phiRb_optimal_allocation(const['gamma_max'], nu_range,
                                                const['Kd_cpc'], const['phi_O'])
opt_lam = ecor.model.steady_state_growth_rate(const['gamma_max'], opt_phiRb, 
                                              nu_range, const['Kd_cpc'], const['phi_O'])

fig, ax = plt.subplots(1, 2, figsize=(4, 2))

strains = {'ECOR02': 'X', 'ECOR39':'D', 'ECOR55': '*', 'ECOR59': 'v', 'ECOR63':'s', 'NCM3722':'o',
           'Isolate_1': '^'}

c = {'ECOR63':cor['light_red'], 'ECOR59':cor['primary_red'],  'NCM3722': cor['red'], 'ECOR02':cor['dark_red']}
for g, d in data[data['carbon_source'] != 'Glucose'].groupby('strain'):
    if g == 'NCM3722':
        color = cor['primary_black']
        alpha = 0.75
    else:
        color = cor['light_black']
        alpha=0.5
    ax[0].plot(d['growth_rate_hr'], d['rna_to_protein'], strains[g],
               color=color, alpha=alpha, ms=5, markeredgecolor=cor['primary_black'])

phiRb_range = np.linspace(0, 0.44, 200)
for g, d in carb_source.groupby('strain'):
    nu_est = estimate_nu(d['rna_to_protein']*0.4558, d['growth_rate_hr']).values[0]
    ax[0].plot(d['growth_rate_hr'], d['rna_to_protein'], strains[g], ms=5, 
               color=c[g], markeredgecolor=cor['primary_black'], label=g)
    ax[1].plot(d['rna_to_protein'], d['growth_rate_hr'], strains[g], ms=5,
               markeredgecolor=cor['primary_black'], color=c[g], label=f'{g}\n' +r'$\nu_{max} \approx$' + f'{int(nu_est)}' + ' hr$^{-1}$' +'\n')
    nu_est = estimate_nu(d['rna_to_protein']*0.4558, d['growth_rate_hr']).values[0]
    lam = ecor.model.steady_state_growth_rate(const['gamma_max'], phiRb_range,
                                              nu_est, const['Kd_cpc'], const['phi_O'])
    ax[1].plot(phiRb_range/0.4558, lam, color=c[g], ls='-', lw=1)

ax[0].plot(opt_lam, opt_phiRb/0.4558, '-', color=cor['primary_blue'], lw=1, label='optimal behavior')
ax[0].set_ylim([0.05, 0.8])
ax[0].set_xlim([0.1, 3])

ax[0].set_xlabel('growth rate [hr$^{-1}$]\n$\lambda$', fontsize=6)
ax[0].set_ylabel('RNA-to-protein ratio', fontsize=6)

ax[1].set_xlim([phiRb_range.min()/0.4558, phiRb_range.max() / 0.4558])
ax[1].set_ylim([0, 1.1])
ax[1].set_xlabel('RNA-to-protein ratio', fontsize=6)
ax[1].set_ylabel('$\lambda$\ngrowth rate [hr$^{-1}$]', fontsize=6)
ax[1].legend(bbox_to_anchor=(1,1))

plt.savefig('./metabolic_rate_example.pdf', bbox_inches='tight')


#%% use a linear fit to the r-line to estimate metabolic rates
import scipy.stats
data = pd.read_csv('../../data/results_summary.csv', skiprows=1)
r_data = pd.read_csv('../../data/Fig4_ecoli_ribosomal_mass_fractions.csv')
popt = scipy.stats.linregress(r_data['growth_rate_hr'], r_data['mass_fraction'])
data = data[['Strain', 'Buffer', 'Experiment', 'carbon_source', 'growth_rate']]
data['estimated_phiRb'] = popt[1] + popt[0] * data['growth_rate']
data['estimated_nu_inv_hr'] = estimate_nu(data['estimated_phiRb'], data['growth_rate'])
data.to_csv('../../data/estimated_metabolic_rate_all_strains.csv', index=False)
