import csv,hashlib
from pathlib import Path
out=Path('output/followup_analysis')
inputs=list(Path('data/reference_followup').glob('*'))+list(Path('data/external_rbp_controls/original').glob('*'))+[Path(p) for p in ['data/MSI2_peaks_omniCLIP.bed','data/results/peak_analysis/peak_count_matrix.tsv','data/results/peak_analysis/bam_library_sizes.tsv','code/clip_peak_analysis/deseq_samples.tsv','data/results/peak_analysis/figures/all_nonzero_peaks/transcript_region_bias/peak_transcript_region_annotations.tsv','data/results/landmark_sites.tsv.gz','data/results_round2/features.tsv.gz','data/results_round2/sample_metadata.tsv','data/results_round2/analysis_config.yaml','data/results_round2/reference_peaks/peaks.tsv']]+list(Path('data/results_round2').glob('NM*/read_annotations.tsv.gz'))+list(Path('data/results_round2').glob('NM*/feature_bin_counts.tsv.gz'))
def save(paths,name):
 rows=[]
 for p in sorted(set(paths)):
  if not p.is_file():continue
  h=hashlib.sha256()
  with p.open('rb') as f:
   for c in iter(lambda:f.read(1024*1024),b''):h.update(c)
  rows.append(dict(path=str(p),bytes=p.stat().st_size,sha256=h.hexdigest()))
 with (out/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=['path','bytes','sha256'],delimiter='\t');w.writeheader();w.writerows(rows)
 print(name,len(rows))
save(inputs,'input_manifest.tsv')
save(list((out/'scripts').glob('*'))+list((out/'figures').glob('*'))+list((out/'tables').glob('*'))+list((out/'raw_count_quantile').glob('*'))+list(Path('output/external_rbp_splice_comparison').glob('*')),'output_manifest.tsv')
