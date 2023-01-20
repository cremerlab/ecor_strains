# %%
import numpy as np
import pandas as pd
import glob
import ecor.io

# Load each of the abundance file paths.
files = glob.glob(
    '../../../data/transcriptomics/2022-12-17_diauxic_shifts/*/abundance.tsv')

# Assemble the total dataframe
rnaseq_df = pd.DataFrame([])
for i, f in enumerate(files):
    date, strain, phase = f.split('/')[-2].split('_')
    df = pd.read_csv(f, delimiter='\t')
    # Get the gene names
    names = [v.split('|')[0] for v in df['target_id'].values]
    # b_num = [v.split('|')[-2] for v in df['target_id'].values]
    std_names = ecor.io.standardize_genes(names)
    df['gene_name'] = std_names['names']
    df['cog_letter'] = std_names['cog_letters']
    df['cog_class'] = std_names['cog_classes']
    df['cog_desc'] = std_names['cog_descs']
    df['allo_sector'] = std_names['sectors']
    df['fraction'] = df['tpm'].values / 1E6
    df['date'] = date
    df['strain'] = strain
    df['growth_phase'] = phase
    df = df[['date', 'strain', 'growth_phase', 'gene_name',
             'cog_letter', 'cog_class', 'cog_desc',
             'allo_sector', 'est_counts', 'tpm',
             'fraction']]
    rnaseq_df = pd.concat([rnaseq_df, df], sort=False)
rnaseq_df.to_csv('output/2022-12-17_RNASEQ_quantification.csv', index=False)
