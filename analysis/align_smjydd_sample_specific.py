"""Map each sample only to the references specified by the user."""
from collections import Counter
import csv
import json
from align_smjydd_references import ROOT, REFS, LEFT, RIGHT, mp, mapper, classify, locus_call, write_table

OUT = ROOT/'outputs/SMJYDD_sample_specific_alignment'
PANELS = {'L2': ['328_ref', '346_MSI2'], 'L4': ['324_altUTR'], 'L8': ['251_Exon1']}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    summaries=[]
    with (OUT/'read_assignments.tsv').open('w',newline='') as f, (OUT/'alignments.paf').open('w') as paf:
        writer=None
        for path in sorted((ROOT/'data/SMJYDD_fastq').glob('*.fastq')):
            group=path.stem.split('_')[2]
            names=PANELS[group]
            maps={n:mapper(REFS[n]) for n in names}
            c=Counter()
            for name,seq,qual in mp.fastx_read(str(path)):
                hits={n:list(m.map(seq,cs=True)) for n,m in maps.items()}
                for n,hs in hits.items():
                    for h in hs:
                        cols=str(h).split('\t');cols[3]=n
                        paf.write('\t'.join([path.stem+':'+name,str(len(seq))]+cols)+'\n')
                call=classify(hits,names)
                c['total_reads']+=1;c['mapping_'+call[0]]+=1
                if group in ('L2','L4'):
                    basis='328_ref' if group=='L2' else '324_altUTR'
                    allele='346_MSI2' if group=='L2' else '324_altUTR'
                    args=dict(reference=REFS[basis],insertion_references=[REFS[allele]])
                    locus,net,err=locus_call(seq,hits[basis],**args)
                    strict=locus_call(seq,hits[basis],flank_identity=0.90,**args)[0]
                else:
                    locus=strict='not_applicable';net=err=''
                c[locus]+=1;c['strict_'+strict]+=1
                row=dict(sample=path.stem,read_id=name,read_length=len(seq),allowed_references='|'.join(names),mapping_call=call[0],compatible_references=call[1],best_score=call[2],alignment_identity=call[3],reference_fraction=call[4],locus_call=locus,net_window_insertion_length=net,window_edit_distance=err,strict_locus_call=strict)
                if writer is None:
                    writer=csv.DictWriter(f,fieldnames=list(row),delimiter='\t');writer.writeheader()
                writer.writerow(row)
            row=dict(sample=path.stem,allowed_references='|'.join(names),total_reads=c['total_reads'])
            for key in ['mapping_'+n for n in REFS]+['mapping_ambiguous','mapping_no_qualifying_alignment','target_present_spanning','target_absent_spanning','other_locus_sequence','conflicting_locus_calls','not_callable','not_applicable','strict_target_present_spanning','strict_target_absent_spanning']:
                row[key]=c[key]
            denom=c['target_present_spanning']+c['target_absent_spanning']
            row['target_present_percent_of_resolved']=round(100*c['target_present_spanning']/denom,3) if denom else 'NA'
            assert sum(c['mapping_'+n] for n in names)+c['mapping_ambiguous']+c['mapping_no_qualifying_alignment']==c['total_reads']
            assert sum(c[k] for k in ['target_present_spanning','target_absent_spanning','other_locus_sequence','conflicting_locus_calls','not_callable','not_applicable'])==c['total_reads']
            summaries.append(row)
            print(path.stem,dict(c),flush=True)
    write_table(OUT/'sample_summary.tsv',summaries)
    (OUT/'reference_panels.json').write_text(json.dumps(PANELS,indent=2)+'\n')


if __name__=='__main__':
    # Confirm that both L4 alleles can be evaluated on the actual L4 reference.
    s=REFS['324_altUTR'];m=mapper(s)
    absent=s[:s.index(LEFT)+40]+s[s.index(RIGHT):]
    for seq,expect in [(s,'target_present_spanning'),(absent,'target_absent_spanning')]:
        for q in [seq,mp.revcomp(seq)]:
            assert locus_call(q,list(m.map(q)),reference=s,insertion_references=[s])[0]==expect
    main()
