"""Independent exact-junction checks on >=500-nt reads; no mapping required."""
from summarize_smjydd_min500 import ROOT, OUT, REFS, MIN_LENGTH, mp, write_table

CANDIDATES=[('L2','328_ref',688,63),('L2','346_MSI2',749,3),
            ('L4','324_altUTR',850,54),('L4','324_altUTR',123,2),
            ('L8','251_Exon1',403,49)]


def main():
    rows=[]; fasta=[]
    for group,name,start,length in CANDIDATES:
        ref=REFS[name]
        alt=ref[:start]+ref[start+length:]
        label=f'{name}_del{start+1}_{start+length}'
        fasta.append(f'>{label}\n{alt}\n')
        left=ref[start-20:start];right=ref[start+length:start+length+20]
        junction=left+right;full=left+ref[start:start+length]+right
        for path in sorted((ROOT/'data/SMJYDD_fastq').glob(f'*_{group}_*.fastq')):
            seqs=[s for _,s,_ in mp.fastx_read(str(path)) if len(s)>=MIN_LENGTH]
            rows.append(dict(sample=path.stem,candidate=label,
                deletion_junction_exact_reads=sum(junction in s or mp.revcomp(junction) in s for s in seqs),
                reference_window_exact_reads=sum(full in s or mp.revcomp(full) in s for s in seqs),
                retained_reads=len(seqs)))
    write_table(OUT/'exact_junction_validation.tsv',rows)
    (OUT/'candidate_alternatives.fasta').write_text(''.join(fasta))


if __name__=='__main__':
    main()
