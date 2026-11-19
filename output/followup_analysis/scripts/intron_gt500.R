suppressPackageStartupMessages({library(data.table); library(ggplot2)})

# Run from the project root. Reuse the exact assignments behind the quartile plot.
out <- 'output/followup_analysis/intron_gt500'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
save_plot <- function(p, name, w=10, h=5.5) {
  for (ext in c('png', 'pdf')) ggsave(file.path(out, paste0(name, '.', ext)), p,
                                    width=w, height=h, dpi=220)
}
hits <- fread('output/followup_analysis/tables/summit_normalized_positions.tsv.gz')[kind=='intron']
introns <- unique(hits[, .(feature_id, length, sites_in_feature)])
stopifnot(!anyDuplicated(introns$feature_id),
          all(hits[, .N, by=feature_id]$N == introns$sites_in_feature))
introns[, retained_gt500 := length > 500]
fwrite(introns, file.path(out, 'included_intron_lengths.tsv'), sep='\t')
summary <- rbindlist(lapply(c(FALSE, TRUE), function(filtered) {
  z <- if (filtered) introns[length>500] else introns
  data.table(population=if(filtered) 'Introns >500 nt' else 'Original plotted introns',
             introns=nrow(z), summits=sum(z$sites_in_feature), min_nt=min(z$length),
             q25_nt=quantile(z$length,.25), median_nt=median(z$length),
             q75_nt=quantile(z$length,.75), max_nt=max(z$length))
}))
fwrite(summary, file.path(out, 'intron_length_summary.tsv'), sep='\t')
save_plot(ggplot(introns, aes(length)) +
  geom_histogram(bins=55, fill='#236d83', color='white', linewidth=.2) +
  geom_vline(xintercept=500, linetype=2, color='#b44d36') +
  scale_x_log10(labels=scales::label_number(big.mark=',')) + theme_bw(base_size=12) +
  labs(title='Lengths of introns represented in the summit analysis',
       subtitle=sprintf('%s unique introns; one count per intron; dashed line = 500 nt',
                        format(nrow(introns),big.mark=',')),
       x='Intron length (nt; logarithmic scale)', y='Number of introns'),
  'intron_length_distribution')

z <- hits[length>500]
n <- nrow(z); nf <- uniqueN(z$feature_id)
stopifnot(n>0, all(z$fraction>=0 & z$fraction<1))
weights <- 1/z$sites_in_feature
profile <- function(bin, w) {
  ans <- numeric(100)
  v <- rowsum(w, bin+1L, reorder=FALSE)
  ans[as.integer(rownames(v))] <- v[,1]
  ans/sum(w)
}
bins <- floor(z$fraction*100)
observed <- cbind(site_weighted=profile(bins,rep(1,n)),
                  intron_weighted=profile(bins,weights))
stopifnot(all(abs(colSums(observed)-1)<1e-10), abs(sum(weights)-nf)<1e-8)
seed <- 260922L; iterations <- 500L; set.seed(seed)
null <- array(0, c(100,iterations,2))
for (i in seq_len(iterations)) {
  # Same independent half-integer grid adaptation as site_spatial.R's UTR null.
  # Does not reproduce the historical simulate() inter-site spacing constraint.
  random_position <- 12.5 + floor(runif(n)*(floor(z$length-25)+1))
  b <- pmin(99L, floor(100*random_position/z$length))
  null[,i,1] <- profile(b,rep(1,n))
  null[,i,2] <- profile(b,weights)
}
pr <- rbindlist(lapply(1:2, function(j) data.table(
  position=(0:99)+.5, weighting=colnames(observed)[j], observed=observed[,j],
  background_median=apply(null[,,j],1,median),
  background_low=apply(null[,,j],1,quantile,.025),
  background_high=apply(null[,,j],1,quantile,.975),
  summits=n, introns=nf, iterations=iterations, seed=seed)))
fwrite(pr,file.path(out,'intron_gt500_profiles.tsv'),sep='\t')
saveRDS(null,file.path(out,'PARpipe_intron_gt500_nulls.rds'))
subtitle <- sprintf('%s summits in %s introns; strictly >500 nt; original assignments retained',
                    format(n,big.mark=','),format(nf,big.mark=','))
labels <- c(site_weighted='Site weighted', intron_weighted='Equal intron weight')
save_plot(ggplot(pr,aes(position,observed,color=weighting)) + geom_line(linewidth=.6) +
  scale_color_manual(values=c(site_weighted='#F8766D',intron_weighted='#00BFC4'),labels=labels) +
  theme_bw(base_size=12) + theme(legend.position='bottom') +
  labs(title='Intron summit density: introns longer than 500 nt',subtitle=subtitle,
       x='Intron position (%)',y='Fraction of summits per 1% bin',color=NULL),
  'metaintron_site_weighting_gt500')
# Panel-D-style background percentiles; same null draws and observed profiles.
for (j in 1:2) {
  ix <- which(pr$weighting == colnames(observed)[j])
  for (prob in c(0,.10,.25,.75,.90,1)) {
    column <- paste0('null_p',sprintf('%02d',round(prob*100)))
    set(pr, i=ix, j=column, value=apply(null[,,j],1,quantile,prob))
  }
}
fwrite(pr,file.path(out,'intron_gt500_profiles.tsv'),sep='\t')
save_plot(ggplot(pr,aes(position,observed)) +
  geom_ribbon(aes(ymin=null_p00,ymax=null_p100,fill='Outer range (min-max)')) +
  geom_ribbon(aes(ymin=null_p10,ymax=null_p90,fill='10th-90th percentiles')) +
  geom_ribbon(aes(ymin=null_p25,ymax=null_p75,fill='25th-75th percentiles')) +
  geom_line(aes(y=null_p25),color='grey45',linewidth=.3) +
  geom_line(aes(y=null_p75),color='grey45',linewidth=.3) +
  geom_line(aes(y=background_median,color='Background median'),linewidth=.5) +
  geom_line(aes(color='Observed MSI2'),linewidth=.7) +
  scale_fill_manual(name='Randomized background',
    breaks=c('25th-75th percentiles','10th-90th percentiles','Outer range (min-max)'),
    values=c('25th-75th percentiles'='#4682B4','10th-90th percentiles'='#ADD8E6',
             'Outer range (min-max)'='#DEDEDE')) +
  scale_color_manual(name=NULL,values=c('Background median'='black','Observed MSI2'='#D62728')) +
  facet_wrap(~weighting,labeller=as_labeller(labels)) + theme_bw(base_size=12) +
  theme(legend.position='bottom',legend.box='vertical',panel.grid.minor=element_blank()) +
  labs(title='Intron summit positions versus PARPIPE-style background',subtitle=subtitle,
       x='Intron position (%)',y='Weighted fraction per 1% bin',
       caption='Pointwise percentiles of 500 randomizations; colors describe the null distribution, not p-value thresholds.\nSame 12.5 nt end-margin adaptation used for the previous 3-prime UTR analysis.'),
  'PARpipe_intron_background_gt500',12,6.5)
writeLines(c('# Intron summit analysis: >500 nt follow-up', '',
  'Reproduce from the project root: `Rscript output/followup_analysis/scripts/intron_gt500.R`.', '',
  'The length histogram counts each unique intron represented in the original summit-quartile figure once (not every annotated intron and not once per summit). The table includes all original represented introns and identifies those retained by length >500 nt. The x axis is logarithmic.', '',
  'The new profiles filter the existing longest-matching, strand-aware intron assignments; they do not reassign summits. Position 0% is the transcript-oriented intron start and 100% the end. All 100 bins are retained, including zero-count bins. Site weighting gives every summit equal weight. Equal-intron weighting gives each represented intron total weight one.', '',
  'The background reproduces the prior 3′ UTR adaptation in site_spatial.R: independent positions on the half-integer grid 12.5:(L-12.5), 500 draws, preserving each assigned intron length and summit count. Both weighting schemes use identical randomized positions. The panel-D-style plot shows the background median in black, the 25th-75th percentiles in dark blue, the 10th-90th percentiles in light blue, and the remaining range to the minimum/maximum in gray; observed MSI2 is red. These are pointwise null percentiles, not p-value thresholds or observed-curve confidence intervals. The earlier 2.5th/97.5th percentiles remain in the data table. Color interpretation follows Figure 1 and the spatial-pattern methods of https://doi.org/10.1186/gb-2014-15-1-r12; gray is rendered as an outer envelope here.', '',
  'This is a PARPIPE-style spatial-background adaptation, not execution of the full historical PARPIPE script. In particular, the historical simulate() function also contains inter-site spacing logic that the previous UTR adaptation did not use; this follow-up preserves that previous approach for comparability. It does not model local mappability or sequence composition.', '',
  capture.output(print(summary))),file.path(out,'README.md'))
print(summary)
