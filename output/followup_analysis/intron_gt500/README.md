# Intron summit analysis: >500 nt follow-up

Reproduce from the project root: `Rscript output/followup_analysis/scripts/intron_gt500.R`.

The length histogram counts each unique intron represented in the original summit-quartile figure once (not every annotated intron and not once per summit). The table includes all original represented introns and identifies those retained by length >500 nt. The x axis is logarithmic.

The new profiles filter the existing longest-matching, strand-aware intron assignments; they do not reassign summits. Position 0% is the transcript-oriented intron start and 100% the end. All 100 bins are retained, including zero-count bins. Site weighting gives every summit equal weight. Equal-intron weighting gives each represented intron total weight one.

The background reproduces the prior 3′ UTR adaptation in site_spatial.R: independent positions on the half-integer grid 12.5:(L-12.5), 500 draws, preserving each assigned intron length and summit count. Both weighting schemes use identical randomized positions. The panel-D-style plot shows the background median in black, the 25th-75th percentiles in dark blue, the 10th-90th percentiles in light blue, and the remaining range to the minimum/maximum in gray; observed MSI2 is red. These are pointwise null percentiles, not p-value thresholds or observed-curve confidence intervals. The earlier 2.5th/97.5th percentiles remain in the data table. Color interpretation follows Figure 1 and the spatial-pattern methods of https://doi.org/10.1186/gb-2014-15-1-r12; gray is rendered as an outer envelope here.

This is a PARPIPE-style spatial-background adaptation, not execution of the full historical PARPIPE script. In particular, the historical simulate() function also contains inter-site spacing logic that the previous UTR adaptation did not use; this follow-up preserves that previous approach for comparability. It does not model local mappability or sequence composition.

                 population introns summits min_nt q25_nt median_nt q75_nt
                     <char>   <int>   <int>  <int>  <num>     <num>  <num>
1: Original plotted introns    9732   20542    101 5329.5     16289  47402
2:          Introns >500 nt    9565   20364    504 5660.0     16857  48286
    max_nt
     <int>
1: 1240120
2: 1240120
