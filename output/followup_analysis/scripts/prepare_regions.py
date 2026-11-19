"""GENCODE v26 transcript-derived segments, 0-based half-open; same coding/UTR rule as round2."""
import gzip,re,csv
from collections import defaultdict
from pathlib import Path
tr={}
with gzip.open('data/reference_followup/gencode.v26.primary_assembly.annotation.gtf.gz','rt') as h:
 for line in h:
  if line.startswith('#'):continue
  f=line.rstrip().split('\t')
  if f[2] not in ('exon','CDS','stop_codon'):continue
  a=dict(re.findall(r'(\w+) "([^"]*)"',f[8])); t=a.get('transcript_id')
  if not t:continue
  d=tr.setdefault(t,{'chr':f[0],'strand':f[6],'gene':a['gene_id'],'exon':[],'coding':[]})
  d['exon' if f[2]=='exon' else 'coding'].append((int(f[3])-1,int(f[4])))
rows=set()
for t,d in tr.items():
 e=sorted(set(d['exon'])); c=d['coding']; prefix=(d['chr'],d['strand'],d['gene'])
 for (a,b),(x,y) in zip(e,e[1:]):
  if x>b:rows.add(prefix+(b,x,'intron'))
 if c:
  lo=min(x[0] for x in c); hi=max(x[1] for x in c)
  for a,b in c:rows.add(prefix+(a,b,'CDS'))
  for a,b in e:
   if a<lo: rows.add(prefix+(a,min(b,lo),'5UTR' if d['strand']=='+' else '3UTR'))
   if b>hi: rows.add(prefix+(max(a,hi),b,'3UTR' if d['strand']=='+' else '5UTR'))
 else:
  for a,b in e:rows.add(prefix+(a,b,'noncoding_exon'))
with gzip.open('output/followup_analysis/cache/gencode_region_segments.tsv.gz','wt') as h:
 w=csv.writer(h,delimiter='\t');w.writerow(['chrom','strand','gene_id','start','end','category']);w.writerows(sorted(rows))
print(len(rows),'segments')
