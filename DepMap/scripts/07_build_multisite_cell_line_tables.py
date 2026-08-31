from pathlib import Path
import pandas as pd

ROOT = Path('/Users/jingruimu/Documents/ChatGPT/federated AI/DepMap')
model = pd.read_csv(ROOT/'Model.csv', low_memory=False)
maf = pd.read_csv(ROOT/'OmicsSomaticMutationsMAF.maf', sep='\t', usecols=['DepMap_ID','ModelID','Hugo_Symbol','Variant_Classification','Protein_Change'], low_memory=False)
prism = pd.read_csv(ROOT/'Drug_sensitivity_AUC_(PRISM_Repurposing_Secondary_Screen)_subsetted.csv', usecols=[0,1,2,3,4,5,6], low_memory=False)
prism.columns = ['DepMap_ID','cell_line_display_name','lineage_1','lineage_2','lineage_3','lineage_6','lineage_4']
prism_ids = set(prism.DepMap_ID.dropna())

pairs = [('SKCM','BRAF'), ('BRCA','ERBB2'), ('COAD','BRAF'), ('GB','EGFR')]
non_syn = ~maf.Variant_Classification.fillna('').isin(['Silent','Synonymous_SNP','RNA','Intron','3\'UTR','5\'UTR','IGR'])
maf = maf[non_syn]

for disease, gene in pairs:
    m = model[(model.OncotreeCode == disease) & model.ModelID.isin(prism_ids)].copy()
    g = maf[maf.Hugo_Symbol.eq(gene)].copy()
    # The MAF contains both a sequencing identifier (DepMap_ID, often CDS-*)
    # and the analysis model identifier (ModelID, ACH-*); ModelID joins Model.csv/PRISM.
    by_id = g.groupby('ModelID').agg(
        variants=('Protein_Change', lambda x: ';'.join(sorted(set(str(v) for v in x.dropna())))),
        mutation_rows=('DepMap_ID','size')
    ).reset_index()
    out = m[['ModelID','CellLineName','OncotreeCode','OncotreeLineage','OncotreePrimaryDisease']].rename(columns={'ModelID':'DepMap_ID'})
    out = out.merge(by_id, left_on='DepMap_ID', right_on='ModelID', how='left').drop(columns=['ModelID'], errors='ignore')
    out['mutated'] = out.mutation_rows.notna().astype(int)
    out['variants'] = out.variants.fillna('')
    out['prism_available'] = 1
    out = out.sort_values(['mutated','DepMap_ID'], ascending=[False,True])
    path = ROOT/'processed'/f'{disease.lower()}_{gene.lower()}_cell_line_table.tsv'
    out.to_csv(path, sep='\t', index=False)
    print(f'{disease}-{gene}: {len(out)} lines; mutated={out.mutated.sum()}; wild_type={(out.mutated==0).sum()}; output={path}')

