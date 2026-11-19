"""Prepare independent per-reference IGV views from >=500-nt reads.

Run with system Python (mappy dependencies installed for Python 3.8).
Then run convert_smjydd_igv_bam.py with the bundled Python 3.12.
"""
import csv
from pathlib import Path
from xml.etree import ElementTree as ET
from align_smjydd_references import ROOT, REFS, mp, mapper, score
from align_smjydd_sample_specific import PANELS

OUT=ROOT/'outputs/SMJYDD_IGV'
TEMP=ROOT/'tmp/SMJYDD_IGV_SAM'


def main():
    OUT.mkdir(parents=True,exist_ok=True);TEMP.mkdir(parents=True,exist_ok=True)
    (OUT/'references.fasta').write_text(''.join(f'>{n}\n{s}\n' for n,s in REFS.items()))
    calls={(r['sample'],r['read_id']):r for r in csv.DictReader((ROOT/'outputs/SMJYDD_min500/read_assignments.tsv').open(),delimiter='\t')}
    summary=[];groups={}
    for path in sorted((ROOT/'data/SMJYDD_fastq').glob('*.fastq')):
        group=path.stem.split('_')[2]
        for refname in PANELS[group]:
            base=path.stem+'__on_'+refname
            sam=TEMP/(base+'.sam');m=mapper(REFS[refname])
            mapped=unmapped=0
            with sam.open('w') as f:
                f.write('@HD\tVN:1.6\tSO:unsorted\n')
                for n,s in REFS.items():f.write(f'@SQ\tSN:{n}\tLN:{len(s)}\n')
                f.write(f'@RG\tID:{path.stem}\tSM:{path.stem}\tPL:ONT\n')
                f.write('@CO\tIndependent single-reference visualization; >=500 nt; one best qualifying alignment per read. Do not sum across reference views.\n')
                for name,seq,qual in mp.fastx_read(str(path)):
                    if len(seq)<500:continue
                    call=calls[(path.stem,name)]
                    tags=[f'RG:Z:{path.stem}',f'ZC:Z:{call["mapping_call"]}',f'ZL:Z:{call["locus_call"]}']
                    hs=[h for h in m.map(seq,MD=True) if h.mlen/h.blen>=.85 and h.q_en-h.q_st>=80]
                    if not hs:
                        unmapped+=1
                        fields=[name,4,'*',0,0,'*','*',0,0,seq,qual]
                    else:
                        mapped+=1;h=max(hs,key=score)
                        oriented=seq if h.strand==1 else mp.revcomp(seq)
                        oq=qual if h.strand==1 else qual[::-1]
                        left=h.q_st if h.strand==1 else len(seq)-h.q_en
                        right=len(seq)-h.q_en if h.strand==1 else h.q_st
                        cigar=(f'{left}S' if left else '')+h.cigar_str+(f'{right}S' if right else '')
                        assert sum(n for n,op in h.cigar if op in (0,1,7,8))+left+right==len(seq)
                        assert sum(n for n,op in h.cigar if op in (0,2,3,7,8))==h.r_en-h.r_st
                        fields=[name,16 if h.strand==-1 else 0,refname,h.r_st+1,h.mapq,cigar,'*',0,0,oriented,oq]
                        tags += [f'NM:i:{h.NM}',f'AS:i:{score(h)}',f'MD:Z:{h.MD}']
                    f.write('\t'.join(map(str,fields+tags))+'\n')
            summary.append(dict(sample=path.stem,view_reference=refname,mapped=mapped,unmapped=unmapped,total=mapped+unmapped,bam=base+'.bam'))
            groups.setdefault((group,refname),[]).append(base+'.bam')
    with (OUT/'view_counts.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(summary[0]),delimiter='\t');w.writeheader();w.writerows(summary)
    bed=[('328_ref',688,751,'Candidate_63bp_deletion'),('346_MSI2',749,752,'Candidate_3bp_junction_deletion'),('346_MSI2',753,810,'57bp_insertion_vs_328'),('324_altUTR',849,906,'Original_57bp_target'),('324_altUTR',852,906,'Candidate_54bp_deletion_equivalent_breakpoints'),('324_altUTR',123,125,'Candidate_2bp_deletion'),('251_Exon1',403,452,'Candidate_49bp_deletion')]
    (OUT/'regions_of_interest.bed').write_text(''.join('\t'.join(map(str,x))+'\n' for x in bed))
    loci={'328_ref':'328_ref:630-810','346_MSI2':'346_MSI2:680-860','324_altUTR':'324_altUTR:790-970','251_Exon1':'251_Exon1:1-744'}
    for (group,refname),bams in groups.items():
        root=ET.Element('Session',genome=str(OUT/'references.fasta'),locus=loci[refname],version='3')
        resources=ET.SubElement(root,'Resources')
        ET.SubElement(resources,'Resource',path=str(OUT/'regions_of_interest.bed'),name='Candidate regions (not confirmed biological variants)')
        for bam in bams:ET.SubElement(resources,'Resource',path=str(OUT/bam),index=str(OUT/(bam+'.bai')),name=bam.replace('.bam',''))
        ET.ElementTree(root).write(OUT/f'{group}_on_{refname}.xml',encoding='UTF-8',xml_declaration=True)
    print(summary)


if __name__=='__main__':main()
