# Round-two comparisons

Open `report.html` for a visual index of every plot. All figures are available as PNG and vector PDF in `figures/`; underlying values and comparisons are in `tables/`.

Run from the repository root: `Rscript analysis/compare_round2.R`.

## Main observations

- Equal-genome sampling preserves the broad normalized metaprofile shapes: across all four group means and both motif sets, correlation with the matching full-genome profile is at least 0.999 for donor/acceptor windows, 0.998 for transcript ends, 0.966 for 3′UTRs, and 0.916 for introns. Relative RMS changes remain below 4%. Peak-only selection changes the question and substantially changes the profiles.
- In genome-wide standardized profiles, 328 CLIP is higher across the 3′UTR, whereas 324 CLIP is higher across the intron. SMI differences are visible too: CLIP coverage alone is not a measure of construct-specific enrichment. The four profiles largely converge after peak selection, with greater replicate variation.
- Expanded donor/acceptor windows show sharp coverage rises on the exonic side of the splice boundary in both CLIP and SMI. UAG filtering reduces these boundary spikes in the genome profiles. The transcript-end profiles have strong upstream SMI peaks, especially for 328; these are not CLIP-specific binding peaks.
- Continuous peak preference is strongly sensitive to sparse support and denominator choice. For all reads, median 3′UTR delta-delta is −2.293 in full genome, −1.251 in equal genome, −5.572 in full peaks, −0.148 in peaks after equal genome, and 0 in equal peaks. Negative values favor 328 relative CLIP/SMI enrichment. These shifts must not be interpreted as a biological effect of normalization.
- Full-genome and full-peak variants have **identical per-peak counts**; only the selected-library CPM denominator changes. The same holds for equal genome versus peaks after equal genome. These identities are checked by the script. Thus the corresponding continuous-plot changes are normalization effects, even though the metaprofiles can also change because non-peak reads are removed.
- Of 16,711 peaks supported in the full all-read libraries, 4,032 have no support after equal-genome sampling, and 6,357 have no support after UAG equal-genome sampling. Zero-support peaks are retained as delta-delta zero in the fixed-universe plot. A separate supported-only plot shows the consequence of excluding them, but then peak membership differs by variant.
- Even in the full all-read variant, 80.5% of peaks have zero 324 CLIP counts in both replicates, versus 42.7% for 328 CLIP. The 0.5-CPM pseudocount therefore materially affects preference magnitudes. A 0.1-CPM sensitivity column is included in the tables. Broad genome metaprofiles and sparse peak preference distributions weight observations differently and need not have the same apparent direction.

## Plot definitions

Every main panel grid has five rows (full genome, equal genome, full peaks, peaks after equal genome, equal peaks) and two columns (all reads and UAG reads). Axes are shared across the ten variants within each figure.

**Continuous preference:** arithmetic mean replicate CLIP and SMI CPM, then `log2((CLIP + 0.5)/(SMI + 0.5))` for each construct; delta-delta is 324 minus 328. CPM denominators are the unique selected-read totals, never summed per-peak counts. The fixed universe is 16,711 peaks with at least one count across all eight all-read full-genome libraries. Round-two exclusive peak annotations use CDS > 5UTR > 3UTR > intron > other. This is not an exact recreation of the older 16,046-peak universe or its annotation method. Violin widths are normalized independently and do not encode region size; region sizes are 13,031 3′UTR, 2,290 CDS, 1,222 intron, 106 other, and 62 5′UTR.

**Expanded windows:** splice donor and acceptor −1,000 to +1,000 nt; annotated transcript end −500 to +500 nt, using the supplied 10-nt bins and transcript orientation. Donor zero is the first intronic base; acceptor zero is the last intronic base; transcript-end zero is the last annotated transcribed base, not a validated cleavage site.

**Standardized frames:** intron and spliced 3′UTR, 100 bins spanning 0–100% from 5′ to 3′. 3′UTR introns are excluded. Profiles retain the fixed annotation universe, including zero-coverage features, and exclude zero-width bins from the relevant denominator. They are not scaled to unit height or unit area.

**Normalized metaprofiles:** supplied mean per-feature depth per million selected reads. Thin lines show biological replicates A/B; thick lines are their arithmetic means. CLIP is solid and SMI dashed. Raw mean-depth versions are provided separately. Peak-selected reads contribute their full aligned blocks, not only bases within a peak interval. Motif presence is in the read, not necessarily at the plotted position.

**Construct-difference profiles:** `(324 CLIP − 324 SMI) − (328 CLIP − 328 SMI)` on the normalized mean-depth scale. Thin lines show A/B labels across constructs; they are descriptive combinations, not an assertion of experimental pairing across constructs. These additive differences use no pseudocount and are different from the continuous peak log-ratio statistic.

## Validation and limits

The script checks all 80 selected-library totals, exclusive-category sum conservation, 56,000 unique sample/variant/feature/bin rows, invariant feature/bin denominators, complete peak matrices, finite preference values, complete region annotations, and paired genome/peak count identities. Source BAM extraction is not rerun or independently audited here. Figures were visually inspected. This is a descriptive comparison of two biological replicates per construct/assay and one sampling seed; no new differential tests, confidence intervals, or sampling-variability estimates are claimed. Changes in the number of unsupported peaks and zero-to-nonzero transitions are included in the sign-change statistic.
