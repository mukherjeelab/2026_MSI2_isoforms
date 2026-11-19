"""Sample-specific analysis; all downstream denominators use reads >=500 nt."""
import csv
from collections import Counter, defaultdict
from statistics import mean, median
from align_smjydd_references import ROOT, REFS, mp, mapper, classify, score, locus_call, write_table
from align_smjydd_sample_specific import PANELS

MIN_LENGTH = 500
OUT = ROOT/'outputs/SMJYDD_min500'


def normalize(ref, pos, deleted, inserted):
    allele = deleted or inserted
    while pos > 0 and allele and ref[pos-1] == allele[-1]:
        allele = ref[pos-1] + allele[:-1]
        pos -= 1
    return pos, allele if deleted else '', allele if inserted else ''


def events(seq, qual, hit, ref):
    s=seq if hit.strand==1 else mp.revcomp(seq)
    qs=qual if hit.strand==1 else qual[::-1]
    q=hit.q_st if hit.strand==1 else len(seq)-hit.q_en
    r=hit.r_st
    found=set(); covered=set()
    for length,op in hit.cigar:
        if op in (0,7,8):
            for j in range(length):
                if ord(qs[q+j])-33 >= 20:
                    covered.add(r+j)
                    if s[q+j] in 'ACGT' and s[q+j]!=ref[r+j]:
                        found.add((r+j,ref[r+j],s[q+j]))
            q+=length;r+=length
        elif op==1:
            found.add(normalize(ref,r,'',s[q:q+length]));q+=length
        elif op==2:
            found.add(normalize(ref,r,ref[r:r+length],''));r+=length
    return found, covered


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    summaries=[]; variants=[]; readrows=[]
    for path in sorted((ROOT/'data/SMJYDD_fastq').glob('*.fastq')):
        group=path.stem.split('_')[2];names=PANELS[group]
        maps={n:mapper(REFS[n]) for n in names}
        counts=Counter();lengths=[]
        evcounts=defaultdict(Counter);evstrands=defaultdict(lambda:defaultdict(Counter))
        coverage=defaultdict(Counter);spans=defaultdict(list)
        for name,seq,qual in mp.fastx_read(str(path)):
            if len(seq)<MIN_LENGTH:
                continue
            lengths.append(len(seq))
            hits={n:list(m.map(seq,cs=True)) for n,m in maps.items()}
            call=classify(hits,names);counts['mapping_'+call[0]]+=1
            locus='not_applicable';net=err=''
            if group in ('L2','L4'):
                basis='328_ref' if group=='L2' else '324_altUTR'
                allele='346_MSI2' if group=='L2' else '324_altUTR'
                locus,net,err=locus_call(seq,hits[basis],reference=REFS[basis],insertion_references=[REFS[allele]])
            counts[locus]+=1
            row=dict(sample=path.stem,read_id=name,length=len(seq),mapping_call=call[0],compatible_references=call[1],identity=call[3],reference_fraction=call[4],locus_call=locus,net_window_insertion_length=net,window_edit_distance=err,variant_screen_reference='',variant_events='')
            readrows.append(row)
            if call[0] not in names:
                continue
            refname=call[0];ref=REFS[refname]
            eligible=[h for h in hits[refname] if h.mlen/h.blen>=.85 and h.q_en-h.q_st>=80]
            h=max(eligible,key=score)
            if h.mlen/h.blen<.90 or (h.r_en-h.r_st)/len(ref)<.80:
                continue
            counts['variant_screen_'+refname]+=1
            ev,covered=events(seq,qual,h,ref)
            row['variant_screen_reference']=refname
            row['variant_events']='|'.join(f'{pos+1}:{deleted}>{inserted}' for pos,deleted,inserted in sorted(ev))
            coverage[refname].update(covered)
            spans[refname].append((h.r_st,h.r_en))
            for e in ev:
                pos,deleted,inserted=e
                if pos<10 or pos+len(deleted)>len(ref)-10:
                    continue
                if not (len(deleted)==len(inserted)==1) and not (h.r_st<=pos-5 and h.r_en>=pos+max(1,len(deleted))+5):
                    continue
                evcounts[refname][e]+=1
                evstrands[refname][e][h.strand]+=1
        summary=dict(sample=path.stem,retained_reads=len(lengths),mean_length=round(mean(lengths),2),median_length=median(lengths))
        for key in ['mapping_'+n for n in REFS]+['mapping_ambiguous','mapping_no_qualifying_alignment','target_present_spanning','target_absent_spanning','other_locus_sequence','conflicting_locus_calls','not_callable','not_applicable']:
            summary[key]=counts[key]
        for n in REFS:
            summary['variant_screen_'+n]=counts['variant_screen_'+n]
        assert sum(counts['mapping_'+n] for n in names)+counts['mapping_ambiguous']+counts['mapping_no_qualifying_alignment']==len(lengths)
        summaries.append(summary)
        for refname,evs in evcounts.items():
            for (pos,deleted,inserted),n in evs.items():
                if len(deleted)==len(inserted)==1:
                    denom=coverage[refname][pos]
                    kind='substitution'
                else:
                    denom=sum(start<=pos-5 and end>=pos+max(1,len(deleted))+5 for start,end in spans[refname])
                    kind='deletion' if deleted else 'insertion'
                if n<3:continue
                strands=evstrands[refname][(pos,deleted,inserted)]
                variants.append(dict(sample=path.stem,reference=refname,kind=kind,position_1based=pos+1,deleted=deleted,inserted=inserted,support_reads=n,coverage_reads=denom,percent=round(100*n/denom,2) if denom else '',forward_reads=strands[1],reverse_reads=strands[-1]))
        print(summary,flush=True)
    write_table(OUT/'sample_summary.tsv',summaries)
    write_table(OUT/'read_assignments.tsv',readrows)
    variants.sort(key=lambda r:(r['sample'],-r['support_reads']))
    write_table(OUT/'recurring_variants.tsv',variants)


if __name__=='__main__':
    main()
