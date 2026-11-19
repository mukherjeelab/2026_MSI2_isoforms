"""Nanopore reference comparison and anchored target-locus calls.

Dependencies: biopython and mappy (installed in tmp/sequence_mapping_deps).
Run: python3 analysis/align_smjydd_references.py
"""
import csv
import json
import sys
from collections import Counter
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tmp/sequence_mapping_deps'))
import Bio
from Bio import Align, SeqIO
import mappy as mp

OUT = ROOT / 'outputs/SMJYDD_reference_alignment'
TARGET = 'TCAGTCCTGAATAGCTACAGTGCTCAACCGAATTTTGGCGCGCCCGCTTCCCCGGCA'
REFS = {r.id: str(r.seq) for r in SeqIO.parse(ROOT / 'data/SMJYDD_expected_references.fasta', 'fasta')}
BASE = REFS['328_ref']
LEFT, RIGHT = BASE[713:753], BASE[753:793]
EDIT = Align.PairwiseAligner(mode='global', match_score=0, mismatch_score=-1, open_gap_score=-1, extend_gap_score=-1)


def distance(a, b):
    return int(-EDIT.score(a, b)) if a and b else max(len(a), len(b))


def mapper(seq):
    # Standard Nanopore preset; each reference is indexed independently so
    # secondary-hit suppression cannot arbitrarily break cross-reference ties.
    return mp.Aligner(seq=seq, preset='map-ont', best_n=10)


def score(hit):
    # Reconstruct the map-ont dual-affine alignment score from its CIGAR.
    gaps = [(n, op) for n, op in hit.cigar if op in (1, 2)]
    mismatches = hit.NM - sum(n for n, _ in gaps)
    return 2 * hit.mlen - 4 * mismatches - sum(min(4 + 2*n, 24 + n) for n, _ in gaps)


def locus_call(seq, hits, flank_identity=0.85, reference=BASE, insertion_references=None):
    calls = []
    insertion_references = insertion_references or (REFS['324_altUTR'], REFS['346_MSI2'])
    inserts = [s[s.index(LEFT)+40:s.index(RIGHT)] for s in insertion_references]
    lstart, rstart = reference.index(LEFT), reference.index(RIGHT)
    boundaries = (lstart, lstart+29, rstart+10, rstart+39)
    for hit in hits:
        if hit.r_st > lstart or hit.r_en < rstart+40:
            continue
        oriented = seq if hit.strand == 1 else mp.revcomp(seq)
        q = hit.q_st if hit.strand == 1 else len(seq)-hit.q_en
        r = hit.r_st
        positions = {}
        for length, op in hit.cigar:
            if op in (0, 7, 8):
                for pos in boundaries:
                    if r <= pos < r+length:
                        positions[pos] = q + pos-r
                q += length
                r += length
            elif op == 1:
                q += length
            elif op in (2, 3):
                r += length
        if len(positions) != 4:
            continue
        left = oriented[positions[boundaries[0]]:positions[boundaries[1]]+1]
        right = oriented[positions[boundaries[2]]:positions[boundaries[3]]+1]
        if any(1-distance(a,b)/max(len(a),len(b)) < flank_identity for a,b in [(left,LEFT[:30]),(right,RIGHT[-30:])]):
            continue
        window = oriented[positions[boundaries[0]]:positions[boundaries[3]]+1]
        present_err = min(distance(window, LEFT+allele+RIGHT) for allele in inserts)
        absent_err = distance(window, LEFT+RIGHT)
        err = min(present_err, absent_err)
        # Compare the entire anchored window: this avoids indel placement
        # ambiguity caused by repeats immediately beside the insertion.
        if absent_err <= 8 and present_err-absent_err >= 10:
            call = 'target_absent_spanning'
        elif present_err <= 14 and absent_err-present_err >= 10:
            call = 'target_present_spanning'
        else:
            call = 'other_locus_sequence'
        calls.append((call, len(window)-80, err))
    if not calls:
        return 'not_callable', '', ''
    if len({c[0] for c in calls}) > 1:
        return 'conflicting_locus_calls', '', ''
    return min(calls, key=lambda c: c[2])


def classify(alignments, names):
    best = sorted([(score(h), name, h) for name in names for h in alignments[name]
                   if h.mlen/h.blen >= 0.85 and h.q_en-h.q_st >= 80], reverse=True, key=lambda x:x[0])
    if not best:
        return 'no_qualifying_alignment', '', '', '', ''
    s, name, h = best[0]
    compatible = sorted({n for sc,n,_ in best if sc > s-20})
    label = name if len(compatible)==1 else 'ambiguous'
    return label, '|'.join(compatible), s, round(h.mlen/h.blen,4), round((h.r_en-h.r_st)/h.ctg_len,4)


def write_table(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # Add the observed L4 H295R consensus as an explicitly labeled diagnostic
    # alternative, not as a user-supplied expected construct.
    consensus = {p.stem: str(next(SeqIO.parse(p, 'fasta')).seq) for p in sorted((ROOT/'data/SMJYDD_results/SMJYDD_fasta-files').glob('*.fasta'))}
    observed = mp.revcomp(consensus['SMJYDD_4_L4_H295R'])
    all_refs = dict(REFS, observed_L4_H295R_consensus=observed)
    (OUT/'references_with_observed_variant.fasta').write_text(''.join('>'+n+'\n'+s+'\n' for n,s in all_refs.items()))
    maps = {n:mapper(s) for n,s in all_refs.items()}
    # Unit checks guard reverse-strand and insertion CIGAR coordinate handling.
    for name, expected in [('324_altUTR','target_present_spanning'),('346_MSI2','target_present_spanning'),('328_ref','target_absent_spanning'),('observed_L4_H295R_consensus','target_absent_spanning')]:
        for seq in (all_refs[name],mp.revcomp(all_refs[name])):
            assert locus_call(seq,list(maps['328_ref'].map(seq)))[0] == expected, name
    assert locus_call(REFS['251_Exon1'],list(maps['328_ref'].map(REFS['251_Exon1'])))[0]=='not_callable'
    summary, comparisons = [], []
    GLOBAL = Align.PairwiseAligner(mode='global',match_score=2,mismatch_score=-4,open_gap_score=-6,extend_gap_score=-1)
    with (OUT/'consensus_alignments.txt').open('w') as f:
        for sample,s in consensus.items():
            choices=[(GLOBAL.score(ref,q),name,q,strand) for name,ref in REFS.items() for q,strand in [(s,'+'),(mp.revcomp(s),'-')]]
            _,name,q,strand=max(choices)
            aln=GLOBAL.align(REFS[name],q)[0]
            f.write(f'\n{sample}; reference={name}; consensus_strand={strand}\n{aln}\n')
            comparisons.append(dict(sample=sample,closest_reference=name,consensus_strand=strand,consensus_length=len(s),edit_distance=distance(REFS[name],q)))
    write_table(OUT/'consensus_comparison.tsv',comparisons)
    with (OUT/'read_assignments.tsv').open('w',newline='') as reads, (OUT/'alignments.paf').open('w') as paf:
        writer = None
        for path in sorted((ROOT/'data/SMJYDD_fastq').glob('*.fastq')):
            counts, expected_counts, augmented_counts = Counter(), Counter(), Counter()
            lengths=[]
            for name,seq,qual in mp.fastx_read(str(path)):
                lengths.append(len(seq))
                alignments = {n:list(a.map(seq,cs=True)) for n,a in maps.items()}
                for ref,hits in alignments.items():
                    for hit in hits:
                        # mappy's string is PAF columns 3 onward; replace N/A
                        # from single-sequence indexes with the actual ref name.
                        cols=str(hit).split('\t'); cols[3]=ref
                        paf.write('\t'.join([path.stem+':'+name,str(len(seq))]+cols)+'\n')
                exp=classify(alignments,REFS)
                aug=classify(alignments,all_refs)
                locus,inslen,err=locus_call(seq,alignments['328_ref'])
                strict=locus_call(seq,alignments['328_ref'],0.90)[0]
                exact=TARGET in seq or mp.revcomp(TARGET) in seq
                counts[locus]+=1
                counts['strict_'+strict]+=1
                counts['exact_target']+=exact
                counts['exact_target_not_callable']+=exact and locus=='not_callable'
                expected_counts[exp[0]]+=1; augmented_counts[aug[0]]+=1
                row=dict(sample=path.stem,read_id=name,read_length=len(seq),expected_reference_call=exp[0],compatible_expected_references=exp[1],best_score=exp[2],best_alignment_identity=exp[3],best_reference_fraction=exp[4],augmented_reference_call=aug[0],compatible_augmented_references=aug[1],locus_call=locus,net_window_insertion_length=inslen,window_edit_distance=err,strict_flank_locus_call=strict,exact_target=int(exact))
                if writer is None:
                    writer=csv.DictWriter(reads,fieldnames=list(row),delimiter='\t');writer.writeheader()
                writer.writerow(row)
            row=dict(sample=path.stem,total_reads=len(lengths),median_read_length=median(lengths),reads_at_least_700nt=sum(n>=700 for n in lengths))
            for key in ['target_present_spanning','target_absent_spanning','other_locus_sequence','conflicting_locus_calls','not_callable','exact_target','exact_target_not_callable']:
                row[key]=counts[key]
            denom=counts['target_present_spanning']+counts['target_absent_spanning']
            row['present_percent_of_present_or_absent']=round(100*counts['target_present_spanning']/denom,3) if denom else 'NA'
            for key in ['target_present_spanning','target_absent_spanning','not_callable','other_locus_sequence','conflicting_locus_calls']:
                row['strict_'+key]=counts['strict_'+key]
            for key in list(REFS)+['ambiguous','no_qualifying_alignment']:
                row['expected_'+key]=expected_counts[key]
            for key in list(all_refs)+['ambiguous','no_qualifying_alignment']:
                row['augmented_'+key]=augmented_counts[key]
            summary.append(row)
            print(path.stem,dict(counts),flush=True)
    write_table(OUT/'sample_summary.tsv',summary)
    (OUT/'versions.json').write_text(json.dumps(dict(biopython=Bio.__version__,mappy=mp.__version__,preset='map-ont',reference_score_margin=20,minimum_alignment_identity=0.85,minimum_aligned_query_bases=80,flank_identity=0.85,strict_flank_identity=0.90),indent=2)+'\n')


if __name__=='__main__':
    main()
