"""Alternative semiglobal L2/346 view; original read calls stay unchanged."""
import csv
import re
from align_smjydd_references import ROOT, REFS, Align

TEMP=ROOT/'tmp/SMJYDD_IGV_SAM'
OUT=ROOT/'outputs/SMJYDD_IGV'
REF=REFS['346_MSI2']
ALIGN=Align.PairwiseAligner(mode='global',match_score=2,mismatch_score=-4,open_gap_score=-12,extend_gap_score=-0.2)
ALIGN.target_end_gap_score=0
ALIGN.query_end_gap_score=0


def realign(seq):
    a=ALIGN.align(REF,seq)[0]
    coords=a.coordinates
    blocks=[]
    for j in range(coords.shape[1]-1):
        r,q=map(int,coords[:,j]);end_r,end_q=map(int,coords[:,j+1])
        blocks.append((r,q,end_r,end_q))
    diagonal=[i for i,(r,q,er,eq) in enumerate(blocks) if er>r and eq>q]
    first,last=diagonal[0],diagonal[-1]
    start=blocks[first][0];leading=blocks[first][1];trailing=len(seq)-blocks[last][3]
    ops=[];nm=matches=aligned=run=0;md=[]
    if leading:ops.append((leading,'S'))
    for r,q,er,eq in blocks[first:last+1]:
        if er>r and eq>q:
            assert er-r==eq-q
            ops.append((er-r,'M'));aligned+=er-r
            for x,y in zip(REF[r:er],seq[q:eq]):
                if x==y:run+=1;matches+=1
                else:md.extend([str(run),x]);run=0;nm+=1
        elif er>r:
            ops.append((er-r,'D'));nm+=er-r
            md.extend([str(run),'^'+REF[r:er]]);run=0
        else:
            ops.append((eq-q,'I'));nm+=eq-q
    md.append(str(run))
    if trailing:ops.append((trailing,'S'))
    assert sum(n for n,op in ops if op in 'MIS')==len(seq)
    return start,''.join(f'{n}{op}' for n,op in ops),nm,''.join(md),matches,aligned,leading+trailing,a.score


def main():
    rows=[]
    for sample in ['SMJYDD_1_L2_HEK','SMJYDD_2_L2_H295R']:
        src=TEMP/f'{sample}__on_346_MSI2.sam'
        dst=TEMP/f'{sample}__on_346_MSI2_semiglobal.sam'
        with src.open() as inp,dst.open('w') as out:
            for line in inp:
                if line.startswith('@'):
                    out.write(line);continue
                f=line.rstrip().split('\t')
                if int(f[1])&4:
                    out.write(line);continue
                oldcigar=f[5];seq=f[9]
                assert len(seq)>=500
                start,cigar,nm,md,matches,aligned,clipped,sc=realign(seq)
                oldops=[(int(n),op) for n,op in re.findall(r'(\d+)([MIDNSHP=X])',oldcigar)]
                oldnm=int(next(t[5:] for t in f[11:] if t.startswith('NM:i:')))
                oldmatches=sum(n for n,op in oldops if op=='M')-(oldnm-sum(n for n,op in oldops if op in 'ID'))
                # Do not replace a stronger original alignment with a poor fit.
                accept=matches>=0.95*oldmatches and matches/max(1,aligned)>=.85
                rows.append(dict(sample=sample,read_id=f[0],original_cigar=oldcigar,alternative_cigar=cigar,accepted=int(accept),original_softclip=sum(n for n,op in oldops if op=='S'),alternative_softclip=clipped))
                if accept:
                    f[3]=str(start+1);f[4]='255';f[5]=cigar
                    f=f[:11]+[t for t in f[11:] if not t.startswith(('NM:','MD:','AS:'))]
                    f += [f'NM:i:{nm}',f'MD:Z:{md}',f'AS:i:{round(sc*10)}','ZA:Z:semiglobal_gap_extension_0.2_score_scaled10']
                    out.write('\t'.join(f)+'\n')
                else:out.write(line)
        print(sample,'done',flush=True)
    with (OUT/'semiglobal_changes.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)


if __name__=='__main__':main()
