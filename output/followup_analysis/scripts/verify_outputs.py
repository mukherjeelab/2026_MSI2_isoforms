import csv,gzip,hashlib,json
from pathlib import Path
root=Path('.');out=Path('output/followup_analysis')
def rows(p):
 with (gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)) as h:return list(csv.DictReader(h,delimiter='\t'))
checks=[]
def check(k,v):
 checks.append({'check':k,'passed':bool(v)})
 if not v:raise AssertionError(k)
a=rows(out/'tables/read_annotation_summary.tsv')
for sid in sorted(set(x['sample_id'] for x in a)):
 z=[x for x in a if x['sample_id']==sid];check(sid+' fractional mass=1',abs(sum(float(x['fractional_fraction']) for x in z)-1)<1e-8)
for name in ['raw_count_quantile/continuous_delta_delta_peaks.tsv','tables/downsampling_peak_summary.tsv']:
 z=rows(out/name);check(name+' fixed16046',len(z)==16046 and len(set(x['peak_id'] for x in z))==16046)
z=rows(out/'tables/metaintron_equal_shape_profiles.tsv')
for sid in set(x['sample_id'] for x in z):
 for q in range(5):
  r=[x for x in z if x['sample_id']==sid and int(x['quartile'])==q];check(sid+' shape Q'+str(q),len(r)==100 and abs(sum(float(x['equal_intron_fraction']) for x in r)-1)<1e-8)
z=rows('output/external_rbp_splice_comparison/liftover_qc.tsv')
for r in z:check(r['RBP']+' liftover mass balance',int(r['original_N'])==sum(int(r[k]) for k in ['malformed','unmapped','split_or_multiple','width_QC_rejected','retained_N']))
z=rows('output/external_rbp_splice_comparison/splice_bin_profiles.tsv');check('4 RBP x2 boundaries x200bins',len(z)==1600)
with (out/'tables/validation_checks.tsv').open('w') as h:
 w=csv.DictWriter(h,fieldnames=['check','passed'],delimiter='\t');w.writeheader();w.writerows(checks)
print(len(checks),'checks passed')
