"""Append completed follow-up configurations without rewriting historical entries."""
from pathlib import Path
import csv,shutil,html,json,hashlib,re
ROOT=Path(globals().get('SOURCE_ROOT','.')).resolve();OUT=Path(globals().get('REGISTRY_OUT',ROOT/'output/analysis_comparison_registry')).resolve();FU=ROOT/'output/followup_analysis'
def read(p):
 with open(p) as h:return list(csv.DictReader(h,delimiter='\t'))
def write(p,rs,fields):
 with open(p,'w') as h:
  w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rs)
config=read(FU/'new_analysis_configurations.tsv');old=read(OUT/'analysis_registry.tsv');fields=list(old[0]);baseline=[x for x in old if not x['analysis_id'].startswith('FU-')]
snap=OUT/'pre_followup_snapshot';snap.mkdir(exist_ok=True)
for name in ['analysis_registry.tsv','figure_registry.tsv','report.html','build_metadata.json']:
 if not (snap/name).exists():shutil.copy2(OUT/name,snap/name)
new=[]
for c in config:
 r={k:'' for k in fields};r.update(analysis_id=c['analysis_id'],analysis_family='Follow-up sensitivities and published controls',short_name=c['Analysis'],status='COMPLETED',approx_date='2026-09-21',biological_question='Which conclusions survive normalization, weighting and matched backgrounds?',data_universe=c['Universe'],n_peaks_or_sites=c['N'],inclusion_filter=c['Universe'],library_normalization=c['Normalization'],transformation=c['Normalization'],pseudocount='See exact method in follow-up report',clip_smi_definition=c['Normalization'],replicate_aggregation='Biological replicates retained when available; external archived peak sets pooled/unresolved',construct_comparison=c['Construct_direction'],positional_aggregation=c['Weighting'],annotation_or_landmark='hg38 GENCODEv26; fixed reference peaks and original landmark annotation',statistical_tests='See per-method tables; resampling is technical, not biological',key_effect_size=c['Main_result'],overall_direction=c['Construct_direction'],main_result=c['Main_result'],major_limitation='Sparse reads; transcript ambiguity; endogenous source provenance provisional; methods have distinct denominators',role='Sensitivity / diagnostic',script_path=c['script'],result_table_path=c['result_table'],primary_figure_path=c['figure'],what_changed_from_related_analysis='See followup_analysis_status.tsv; exact old duplicates reused',held_constant='Original reference identifiers and source genome annotation',question_ids='Follow-up tasks1-10',provenance_note='New completed follow-up; source scripts and inputs identified; seeds saved',portable_script_path='sources/'+c['script'],portable_result_table_path='sources/'+c['result_table'],portable_primary_figure_path='sources/'+c['figure']);new.append(r)
write(OUT/'analysis_registry.tsv',baseline+new,fields)
fr=read(OUT/'figure_registry.tsv');ff=list(fr[0]);fr=[r for r in fr if not r['figure_id'].startswith('FU-FIG-')];captions={r['path']:r['caption'] for r in read(FU/'figure_metadata.tsv')}
paths=sorted(list((FU/'figures').glob('*.png'))+list((FU/'raw_count_quantile').glob('*.png'))+list((ROOT/'output/external_rbp_splice_comparison').glob('*.png')))
for i,path in enumerate(paths,1):
 rel=str(path.relative_to(ROOT));cs=[c for c in config if c['figure']==rel]
 if not cs:
  stem=path.stem
  key='FU-EXT-MSI2' if 'external_rbp_splice_comparison' in rel else ('FU-FRACTIONAL' if 'annotation' in stem else ('FU-UTR' if 'UTR' in stem else ('FU-QNRAW' if any(x in stem for x in ['region','delta_delta']) else 'FU-META-DEPTH')))
  cs=[c for c in config if c['analysis_id']==key]
 c=cs[0];r={k:'' for k in ff};r.update(figure_id='FU-FIG-'+str(i).zfill(3),analysis_ids=';'.join(x['analysis_id'] for x in cs),analysis_family='Follow-up sensitivities and published controls',figure_name=path.stem,png_path=rel,pdf_path=str(path.with_suffix('.pdf').relative_to(ROOT)),portable_png_path='sources/'+rel,portable_pdf_path='sources/'+str(path.with_suffix('.pdf').relative_to(ROOT)),universe=c['Universe'],normalization=c['Normalization'],main_quantity=c['Unit'],comparison=c['Construct_direction'],main_result=c['Main_result'],caveat='See follow-up report for denominator, inclusion rule and provenance',script_path=c['script'],result_table_path=c['result_table'],gallery='Yes',caption=captions.get(rel,c['Analysis']+': '+c['Main_result']));fr.append(r)
write(OUT/'figure_registry.tsv',fr,ff)
# Reviewable output snapshots; huge genome/BAM exports and per-read caches remain in the workspace.
for folder in [FU/'scripts',FU/'tables',FU/'figures',FU/'raw_count_quantile',ROOT/'output/external_rbp_splice_comparison']:
 for p in folder.rglob('*'):
  if p.is_file() and p.suffix not in ['.log']:
   dst=OUT/'sources'/p.relative_to(ROOT);dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dst)
for p in FU.glob('*'):
 if p.is_file() and p.suffix in ['.tsv','.html','.md','.json']:
  dst=OUT/'sources'/p.relative_to(ROOT);dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dst)
# Preserve navigable paths from the nested review copy without recursively copying this registry.
copyreport=OUT/'sources/output/followup_analysis/report.html'
if copyreport.exists():
    copyreport.write_text(copyreport.read_text().replace('../analysis_comparison_registry/report.html','../../../report.html'))
for rel in ['output/splicing_first_pass/README.md','output/round2_comparisons/report.html','output/round2_comparisons/README.md']:
    src=ROOT/rel; dst=OUT/'sources'/rel
    if src.exists():
        dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)

# Link to the bundled follow-up report and show the new registry rows in the existing HTML.
p=OUT/'report.html';s=p.read_text();s=re.sub(r'<!-- FOLLOWUP-START -->.*?<!-- FOLLOWUP-END -->','',s,flags=re.S)
block='<!-- FOLLOWUP-START --><section style="padding:24px;background:#e7f1f2;border:2px solid #337b85"><h2>Completed follow-up extension: '+str(len(new))+' new configurations</h2><p>The original registry below remains the historical baseline. The combined TSV now contains '+str(sum(r['status']=='COMPLETED' for r in baseline+new))+' completed configurations plus proposed work. <a href="sources/output/followup_analysis/report.html">Open the integrated follow-up report</a>.</p><table><tr><th>ID</th><th>Analysis</th><th>Result</th></tr>'+''.join('<tr><td>'+html.escape(c['analysis_id'])+'</td><td>'+html.escape(c['Analysis'])+'</td><td>'+html.escape(c['Main_result'])+'</td></tr>' for c in config)+'</table></section><!-- FOLLOWUP-END -->'
s=re.sub(r'(<body[^>]*>)',lambda m:m.group(1)+block,s,count=1);p.write_text(s)
meta_path=OUT/'build_metadata.json';m=json.loads(meta_path.read_text());m['followup_completed_variants']=len(new);m['combined_completed_variants']=sum(r['status']=='COMPLETED' for r in baseline+new);m['followup_figure_records']=len(paths);m['followup_report']='sources/output/followup_analysis/report.html';meta_path.write_text(json.dumps(m,indent=2))
print('Registry extended:',len(new),'new configurations,',len(paths),'new figure pairs')
