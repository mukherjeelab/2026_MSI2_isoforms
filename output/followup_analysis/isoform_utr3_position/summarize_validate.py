from pathlib import Path
import json
import numpy as np
import pandas as pd
p = Path('output/followup_analysis/isoform_utr3_position')
a = pd.read_csv(p/'assigned_read_centers.tsv.gz', sep='\t')
c = pd.read_csv(p/'cohort_memberships.tsv.gz', sep='\t')
pr = pd.read_csv(p/'read_position_profiles.tsv', sep='\t')
r = pd.read_csv(p/'endogenous_assigned_summits.tsv.gz', sep='\t')
f = pd.read_csv(p/'candidate_utr3.tsv.gz', sep='\t').set_index('feature_id')
# Independently reconstruct spliced transcript offsets from the JSON segments.
lookup = {fid: json.loads(row.segments_0based_half_open) for fid, row in f.iterrows()}
for fid, group in a.groupby('feature_id'):
    segs = lookup[fid]
    centers = group.center0.to_numpy()
    offsets = np.full(len(group), -1, dtype=int)
    preceding = 0
    for start, end in segs:
        keep = (centers >= start) & (centers < end)
        assert (offsets[keep] == -1).all()
        offsets[keep] = preceding + centers[keep] - start
        preceding += end - start
    assert (offsets >= 0).all()
    if f.loc[fid, 'strand'] == '-': offsets = preceding - 1 - offsets
    assert np.array_equal(offsets, group.offset0.to_numpy())
    assert np.max(np.abs((offsets+.5)/preceding-group.fraction.to_numpy())) < 1e-12
counts=[]
for cohort in c.cohort.unique():
    z = a[a.feature_id.isin(c.loc[c.cohort==cohort,'feature_id'])].copy()
    for (iso, rep, sid), g in z.groupby(['construct','replicate','sample_id']):
        counts.append(dict(cohort=cohort,construct=str(iso),replicate=rep,sample_id=sid,reads=len(g),features=g.feature_id.nunique()))
        curve = pr[(pr.cohort==cohort)&(pr.construct==iso)].sort_values('position')
        assert len(curve)==100
        weights=1/g.groupby('feature_id').feature_id.transform('size')/g.feature_id.nunique()
        observed=weights.groupby((g.fraction*100).astype(int)).sum().reindex(range(100),fill_value=0).to_numpy()
        assert np.max(np.abs(observed-curve['replicate_'+rep].to_numpy())) < 1e-12
        assert np.max(np.abs((curve.replicate_A+curve.replicate_B)/2-curve.observed)) < 1e-12
    for iso,g in z.groupby('construct'):
        counts.append(dict(cohort=cohort,construct=str(iso),replicate='A+B union',sample_id='combined',reads=len(g),features=g.feature_id.nunique()))
    counts.append(dict(cohort=cohort,construct='324+328',replicate='all four union',sample_id='combined',reads=len(z),features=z.feature_id.nunique()))
counts=pd.DataFrame(counts)
counts.to_csv(p/'analysis_counts.tsv',sep='\t',index=False)
pd.DataFrame([dict(summits=len(r),features=r.feature_id.nunique(),minimum_length_nt=100,weighting='equal feature')]).to_csv(p/'endogenous_counts.tsv',sep='\t',index=False)
w=1/r.groupby('feature_id').feature_id.transform('size')/r.feature_id.nunique()
ref_curve=w.groupby((r.fraction*100).astype(int)).sum().reindex(range(100),fill_value=0).to_numpy()
assert np.max(np.abs(ref_curve-pd.read_csv(p/'endogenous_plot_inputs.tsv',sep='\t').observed.to_numpy()))<1e-12
print(counts.to_string(index=False))
print('Endogenous:',len(r),'summits;',r.feature_id.nunique(),'3UTR features')
(p/'validation.txt').write_text('PASS: Independent reconstruction of all read spliced offsets, 12 observed replicate profiles, six replicate means, and the endogenous equal-feature curve. Agreement tolerance 1e-12. Counts use unique feature IDs and sum distinct library records.\n')
