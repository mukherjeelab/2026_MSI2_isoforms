"""Run with bundled Python 3.12; convert, sort, index, and validate IGV files."""
import sys
import csv
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tmp/sequence_browser_py312'))
import pysam
OUT=ROOT/'outputs/SMJYDD_IGV'
TEMP=ROOT/'tmp/SMJYDD_IGV_SAM'


def main():
    pysam.faidx(str(OUT/'references.fasta'))
    results=[]
    for row in csv.DictReader((OUT/'view_counts.tsv').open(),delimiter='\t'):
        name=row['bam'];sam=TEMP/name.replace('.bam','.sam');bam=OUT/name
        pysam.sort('-o',str(bam),str(sam))
        pysam.index(str(bam))
        pysam.quickcheck(str(bam))
        with pysam.AlignmentFile(bam,'rb') as f:
            records=list(f.fetch(until_eof=True))
            assert len(records)==int(row['total'])
            assert len({r.query_name for r in records})==len(records)
            assert sum(not r.is_unmapped for r in records)==int(row['mapped'])
            assert all(r.query_length>=500 for r in records)
            for r in records:
                if r.is_unmapped:continue
                assert r.reference_name==row['view_reference']
                assert r.infer_query_length()==r.query_length
        with pysam.AlignmentFile(bam,'rb') as f:
            assert len(list(f.fetch(row['view_reference'])))==int(row['mapped'])
        results.append(f'{name}: validated {len(records)} records, {row["mapped"]} mapped, indexed fetch OK')
    (OUT/'validation.txt').write_text('\n'.join(results)+f'\npysam {pysam.__version__}\n')
    print('\n'.join(results))


if __name__=='__main__':main()
