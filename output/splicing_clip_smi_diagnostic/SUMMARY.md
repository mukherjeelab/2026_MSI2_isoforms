# CLIP and SMI component diagnostic

## Result in brief

The splice-site effect is not produced solely by SMI normalization. Higher 328
CLIP coverage is already visible on the exonic side of both donors and
acceptors, while intronic CLIP coverage is nearly identical or slightly higher
for 324. SMI differs as well and strengthens the site-weighted CLIP/SMI contrast.
At the donor -5 and acceptor +15 delta-delta minima, both the CLIP numerator and
SMI background contribute, with the exact SMI log-component accounting for
about 60% of the negative delta.

All calculations use the unchanged retained landmark universe, implicit zeros,
common primary-mapped CPM normalization, transcript-oriented 10-nt bins, and
the original replicate pairing. The existing primary enrichment and
delta-delta values were read from `output/splicing_common_cpm_investigation/`
and were not recalculated or changed.

## Donors

- **CLIP-only:** 328 CLIP is higher on the exon side. Across -50 to 0 nt, mean
  CPM per retained site is 0.928 for 328 and 0.756 for 324. At -5 nt it is
  0.957 versus 0.783, and 4.41% versus 1.80% of sites have a nonzero CLIP count.
  Intronic CLIP is small and slightly higher for 324 (0.122 versus 0.099 CPM).
- **SMI-only:** 328 also has higher arithmetic-mean SMI CPM on the exon side
  (1.098 versus 0.813; at -5 nt, 1.126 versus 0.841). However, 328 SMI is more
  concentrated: at -5 nt it is nonzero at 10.42% of sites, compared with
  18.38% for 324.
- **Effect of normalization:** On the exon side, the raw CLIP separation and
  CLIP/SMI delta point in the same direction, and SMI normalization strengthens
  the difference. On the intron side, normalization reverses the small raw
  CLIP difference and produces a modest 328-favoring enrichment contrast.
- **The -5 nt minimum:** delta-delta is -0.459. Its exact decomposition is
  -0.175 from the CLIP log-component and -0.285 from the signed SMI component,
  approximately 38% and 62% of the total magnitude. The minimum is therefore
  driven by both components, with the larger contribution from SMI.
- **Exonic asymmetry:** The exonic-side localization is already apparent in
  CLIP-only coverage. It is not created by SMI normalization, although the
  normalization makes it stronger.

## Acceptors

- **CLIP-only:** 328 CLIP is higher on the exon side. Across 0 to +50 nt, mean
  CPM is 1.021 for 328 and 0.805 for 324. At +15 nt it is 1.043 versus 0.832,
  with nonzero CLIP at 4.66% versus 1.87% of sites. Intron-side CLIP is almost
  identical (0.079 for 328 and 0.081 for 324).
- **SMI-only:** 328 arithmetic-mean SMI CPM is also higher on the exon side
  (1.191 versus 0.875; at +15 nt, 1.222 versus 0.899). As at donors, it is more
  concentrated: +15-nt SMI is nonzero at 11.05% of 328 sites and 19.11% of 324
  sites.
- **Effect of normalization:** CLIP/SMI normalization strengthens the exonic
  raw CLIP difference and reverses the very small intronic raw CLIP difference.
- **The +15 nt minimum:** delta-delta is -0.473, comprising -0.188 from the
  CLIP log-component and -0.284 from the signed SMI component, approximately
  40% and 60% of the total magnitude. Both contribute, with SMI larger.
- **Exonic asymmetry:** The pattern is visible in CLIP-only coverage and becomes
  more pronounced after SMI normalization.

## Annotated transcript 3' ends

- **CLIP-only:** 328 CLIP is higher throughout the window. Full-window mean CPM
  is 1.747 for 328 and 1.239 for 324.
- **SMI-only:** 328 SMI/background is also higher throughout the window and is
  especially elevated upstream. Full-window means are 2.874 for 328 and 1.671
  for 324.
- **Effect of normalization:** The site-weighted CLIP/SMI result remains in the
  same direction and is stronger than the CLIP log-component alone. At the
  -45-nt minimum, delta-delta is -0.692, with nearly balanced contributions:
  -0.356 from CLIP and -0.336 from SMI.
- **Shape:** Unlike the splice landmarks, this is a broad difference across the
  annotated +/-50-nt window. These remain annotated transcript ends, not
  validated cleavage or PAS coordinates.

## Why higher mean 328 SMI can still strengthen 328 enrichment

Arithmetic mean CPM and the primary site-weighted enrichment summarize the
distribution differently. The primary analysis averages site-level log ratios
with a 0.1 CPM pseudocount. At the splice minima, 328 SMI has higher arithmetic
mean CPM but occurs at a smaller fraction of sites, indicating more concentrated
counts. The exact mean-log SMI term therefore favors 328 relative enrichment,
even though the arithmetic mean SMI CPM is higher for 328. This is a consequence
of sparsity, concentration, and the nonlinear log/pseudocount calculation, not
an algebraic contradiction.

The component columns in the diagnostic TSV make the distinction explicit:

```text
delta_delta = clip_contribution_to_delta + smi_contribution_to_delta
```

Raw CLIP-only and SMI-only CPM differences describe observed read distributions;
they are not interpreted as binding differences on their own.
