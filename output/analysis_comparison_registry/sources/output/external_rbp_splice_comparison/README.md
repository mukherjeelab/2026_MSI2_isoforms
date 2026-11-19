# Published RBP splice-boundary comparison

This is a called-site comparison, not a read-coverage comparison. The primary figure uses MSI2's 43,794 fixed, score>=1, same-strand merged reference peaks, represented by the summit of the highest-score source interval. Their endogenous origin remains provisional: the user believes this is the endogenous set, but no accession/publication has been supplied. No construct read-support filter is applied to these reference peaks.

QKI and PUM2 are the exact requested doRiNA Hafner2010 datasets. Their archived JSON describes PCP contiguous coverage clusters from pooled PAR-CLIP libraries, not summits. ELAVL1 is the exact archived PARCLIP_ELAVL1A_hg19 set; JSON identifies Kishore2011, but does not identify the contributing nuclease/library subset. It must not be equated to GSM714639 alone. All external intervals use floor((start+end-1)/2) as midpoint. Individual external biological replicates cannot be reconstructed from these archived BEDs. Untouched BEDs and JSON are in data/external_rbp_controls/original/.

## Coordinate conversion and annotation

The UCSC hg19ToHg38.over.chain is imported with rtracklayer. Retain exactly one mapped interval and 95–105% original width; reject unmapped, split/multiple, malformed, and substantially truncated/expanded intervals. Strand follows the chain transformation. Counts and retention are in liftover_qc.tsv. The multiple category includes split fragments, not just alternate genomic mappings. Converted coordinates remain in processed/.

All RBPs use the exact original MSI2 `data/results/landmark_sites.tsv.gz`: 287,085 donors and 279,145 acceptors, including unsupported sites. Same-strand nearest donor and nearest acceptor are assigned separately to each peak. A peak therefore contributes at most once per boundary-type profile. Nearest ties choose the lower genomic coordinate. The 10-nt bins are half-open, [-1000,1000); +/-250 figures crop the same bins. Donor negative=exonic side, positive=intronic side. Acceptor negative=intronic side, positive=exonic side. The inherited zero is the first intronic donor nucleotide and LAST intronic acceptor nucleotide; thus acceptor [0,50) contains that one anchor nucleotide and 49 exon-side positions. There is no silent one-nucleotide coordinate shift.

## Matched null and denominator

For all four RBPs, intersect the center/summit with same-strand GENCODEv26 transcript-derived segments. Resolve overlapping categories by CDS > 5UTR > 3UTR > intron > noncoding exon, then longest containing segment, gene ID, segment index. The selected gene, category and contiguous segment stay fixed. Draw a uniform integer center within that segment for each peak in each of 500 iterations, seed260923. Preserve peak identity, number, chromosome, strand, gene context, center category, and original width; only the center is required to lie inside the segment, so the conceptual width can extend outside it. No sequence/mappability or GC matching is imposed. This tests localization conditional on genic context, not global intronic preference. Centers without a same-strand modeled genic segment are excluded from observed and null profiles and counted explicitly in background_eligibility.tsv. No profile eligibility depends on proximity to a splice site.

Observed/expected = observed eligible-peak count in a bin / mean randomized count. Fractions divide by ALL null-eligible peaks of that RBP, including those farther than1kb from either boundary. Zero expected counts yield NA, not infinity. Grey bands are the 2.5–97.5% null count quantiles divided by null mean; they are not confidence intervals on biological replicates. Empirical enrichment p=(1+#null>=observed)/501. Bin BH correction spans all1600 bin tests. Saved null matrices reproduce every profile.

Proximity uses ALL post-conversion peaks, including null-ineligible peaks, and ALL annotated splice sites. +/-50 and +/-250 coverage checks are inclusive of both endpoints; profile bins are half-open. These definitions explain small differences between proximity and binned totals. Multiple nearby sites count for site occupancy, but each peak counts once for peak proximity.

## Interpretation

MSI2's exact donor[-10,0) peak enrichment is about1.10-fold and its broad exon-side donor[-50,0) enrichment0.92-fold. Acceptor[0,50) is0.84-fold. Thus the called reference peaks do not reproduce a strong exon-side spike at the read-profile extrema. QKI shows some broader shape similarity, not a compelling sharp boundary-positive control. The descriptive MSI2/QKI profile correlations over+/-250 are0.43 donor and0.48 acceptor; these are not independent tests or evidence of splicing regulation.

Only1.52% of all donors and2.15% of all acceptors have an MSI2 summit within250nt. QKI is also sparse (0.36%/0.41%); low annotation-wide occupancy alone does not invalidate CLIP signal. Peak-count-normalized proximity is more comparable:9.64% of MSI2 peaks are donor-proximal versus7.78% of QKI peaks.

The profiles are not identical across RBPs: ELAVL1 has a distinct acceptor intron-side feature, whereas PUM2 has weak immediate-boundary enrichment. Several sets show immediate-boundary depletion. This can reflect peak-width/cluster construction, feature annotation, or protocol effects, and the null does not model all of them. MSI2's read-based exon-side feature may reflect coverage/background structure; these peak results do not establish its biological specificity. This analysis calibrates plausibility and does not prove or disprove splicing regulation.

## Reproduce

From the repository root, run `Rscript output/external_rbp_splice_comparison/analyze_external.R`, then `Rscript output/external_rbp_splice_comparison/finish_external.R`. First prepare GENCODE segments using the follow-up prepare_regions.py script. Downloaded originals and reference chain are required. The main follow-up report embeds these results alongside existing MSI2 read-based profiles.

## Primary sources

- [Hafner et al.2010](https://pubmed.ncbi.nlm.nih.gov/20371350/)
- [Kishore et al.2011](https://www.nature.com/articles/nmeth.1608)
- [Mukherjee et al.2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4053807/)
- [doRiNA source archive](https://dorina.mdc-berlin.de/regulators/h_sapiens/hg19/)
- [UCSC hg19-to-hg38 chain](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz)
