# First report: common-CPM landmark reanalysis

## Scope and direction

All eight libraries were placed on one scale using
`raw_count / primary_mapped_reads * 1e6`. Enrichment was calculated within each
matched CLIP/SMI biological replicate as
`log2((CLIP_CPM + p) / (SMI_CPM + p))`, then averaged across the two replicates
per construct. The primary pseudocount was 0.1 CPM; 0.5 CPM was retained as a
sensitivity analysis. Delta-delta is enrichment 324 minus enrichment 328, so a
negative value indicates stronger relative enrichment for 328.

The previous bundle in `output/splicing_first_pass/` was read only for QC
comparison and was not overwritten.

## 1. Common-scale normalization

| Landmark | Previous mean DD | Common-CPM mean DD | Change | DD sign changes | Profile Pearson r | Strongest negative bin |
|---|---:|---:|---:|---:|---:|---:|
| Donor | -0.100 | -0.176 | -0.076 | 0/50 | 0.99999 | -5 nt, retained |
| Acceptor | -0.108 | -0.188 | -0.080 | 0/50 | 0.99996 | +15 nt, retained |
| Annotated transcript 3' end | -0.357 | -0.595 | -0.238 | 0/10 | 0.99965 | -45 nt, unchanged |

There are also zero sign changes in either construct's enrichment profile. The
common CPM scale changes the magnitude materially—means are approximately
0.08 more negative at splice sites and 0.24 more negative at annotated 3' ends—but
not the qualitative result. No bin changes sign, profile shapes are essentially
unchanged, and the donor -5 and acceptor +15 extrema remain.

With the 0.5 CPM sensitivity pseudocount, mean delta-delta is -0.091 donor,
-0.099 acceptor, and -0.331 at annotated 3' ends. No bin changes sign relative
to the 0.1 CPM analysis, and Pearson correlations are at least 0.9993. Thus the
pseudocount affects amplitude, as expected for these sparse counts, but not
direction or shape.

## 2. Replicate-level regional estimates

Every matched A and B construct contrast is negative in all nine regions. The
construct direction is therefore visible independently in both biological
replicates and is not created only by the averaged metaprofile.

| Landmark region | A: 324 - 328 | B: 324 - 328 |
|---|---:|---:|
| Donor exon side | -0.515 | -0.359 |
| Donor intron side | -0.096 | -0.061 |
| Donor full immediate | -0.306 | -0.210 |
| Acceptor intron side | -0.067 | -0.044 |
| Acceptor exon side | -0.540 | -0.383 |
| Acceptor full immediate | -0.303 | -0.214 |
| Transcript 3' end upstream | -0.774 | -0.557 |
| Transcript 3' end downstream | -0.610 | -0.440 |
| Transcript 3' end full region | -0.692 | -0.498 |

The complete table reports the four underlying enrichment estimates—324 A,
324 B, 328 A, and 328 B—rather than only these construct contrasts.

## 3. Side of splice

Side difference is exonic delta-delta minus intronic delta-delta. Bootstrap
intervals are 95% intervals for the mean; the site-level bootstrap resamples
genes as clusters. Paired Wilcoxon tests compare the two sides within the same
site or gene. `d_z` is the paired standardized mean difference, and `r_rb` is
the paired rank-biserial effect.

| Landmark/unit | N sites | N genes | Mean side difference | Median | Bootstrap 95% CI | Paired Wilcoxon p | d_z | r_rb |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Donor/site | 101,162 | 13,128 | -0.359 | 0.000 | [-0.370, -0.346] | <2.2e-308 | -0.205 | -0.299 |
| Donor/gene median | 101,162 | 13,128 | -0.224 | 0.000 | [-0.238, -0.208] | 1.98e-276 | -0.278 | -0.536 |
| Acceptor/site | 100,801 | 12,992 | -0.406 | 0.000 | [-0.419, -0.393] | <2.2e-308 | -0.227 | -0.328 |
| Acceptor/gene median | 100,801 | 12,992 | -0.250 | 0.000 | [-0.264, -0.237] | 3.17e-303 | -0.301 | -0.564 |

Both donor and acceptor effects are concentrated on the exonic side. The zero
medians reflect the large number of sparse sites or genes with zero paired
difference; the negative means, cluster-bootstrap intervals, Wilcoxon tests,
and rank-based effects consistently show an exonic shift. Individual 10-nt
bins were not treated as replicates.

## 4. Coverage dependence and the weighting discrepancy

Coverage strata use total normalized whole-window coverage summed across all
eight samples. The statistic is condition-independent and does not use
construct direction or CLIP enrichment. Balanced quartiles were defined within
each landmark type, with site ID used only to break exact coverage ties.

Site-weighted full-region delta-delta remains negative in every quartile:

| Landmark | Q1 lowest | Q2 | Q3 | Q4 highest |
|---|---:|---:|---:|---:|
| Donor | -0.206 | -0.214 | -0.218 | -0.395 |
| Acceptor | -0.189 | -0.206 | -0.224 | -0.415 |
| Annotated transcript 3' end | -0.472 | -0.528 | -0.424 | -0.957 |

Thus ordinary high-coverage sites do not reverse the site-weighted result; Q4
is actually the most negative stratum. The reversal appears only after pooling
reads within a stratum. Pooled-read full-region delta-delta is negative in Q1,
Q2, and Q3, but positive in Q4: +0.656 donor, +0.453 acceptor, and +0.468 at
annotated transcript 3' ends.

Coverage is strongly concentrated. Q4 contains 25% of sites but 75.5% of donor,
75.7% of acceptor, and 82.9% of transcript-end coverage. The top 1% of sites
alone contributes 17.2%, 17.2%, and 38.4%, respectively. Consequently, the
read-weighted calculation is governed by a small, extreme upper tail, whereas
the site-weighted calculation gives each landmark equal influence.

The two summaries also apply the log ratio at different stages: site weighting
averages site-level log ratios with many exact-zero combinations, while the
aggregate method sums reads first and then takes a log ratio. Because log ratios
and pseudocounts are nonlinear, these estimands need not agree. The data localize
the directional reversal to pooled reads in Q4, not to a general positive
coverage trend across sites.

## Conclusion

Common primary-mapped CPM strengthens the magnitude but preserves the original
site-weighted direction and landmark shapes. The 328-over-324 effect is
reproduced in both replicate contrasts, is strongly exonic at donor and
acceptor sites, and remains negative across all site-weighted coverage strata.
The opposing read-weighted summary is best interpreted as an upper-tail pooled-
read effect rather than the typical landmark response. No transcript-level
differential splicing or PSI analysis was performed.
