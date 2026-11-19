suppressPackageStartupMessages({library(data.table); library(ggplot2)})

# Run from the project root. This script reads cached results only.
out <- 'output/followup_analysis/figure_restyles'
observed_path <- 'output/followup_analysis/intron_gt500/intron_gt500_profiles.tsv'
null_path <- 'output/followup_analysis/intron_gt500/PARpipe_intron_gt500_nulls.rds'
analysis_script <- 'output/followup_analysis/scripts/intron_gt500.R'
script_path <- file.path(out, 'plot_intron_null_cloud_percentiles.R')
sources <- c(observed_path, null_path, analysis_script)
hash_before <- tools::md5sum(sources)
dir.create(out, recursive=TRUE, showWarnings=FALSE)
profiles <- fread(observed_path)
null <- readRDS(null_path)
# The generating script stores site_weighted in slice 1, intron_weighted in slice 2.
weight_order <- c('site_weighted', 'intron_weighted')
stopifnot(identical(dim(null), c(100L,500L,2L)), all(is.finite(null)),
          all(null>=0), all(profiles$iterations==500),
          all(profiles$summits==20364), all(profiles$introns==9565),
          all(profiles$seed==260922))
cloud <- rbindlist(lapply(seq_along(weight_order), function(j) {
  p <- profiles[weighting==weight_order[j]][order(position)]
  stopifnot(nrow(p)==100, identical(p$position, (0:99)+.5),
            max(abs(p$background_median-apply(null[,,j],1,median)))<1e-12,
            abs(sum(p$observed)-1)<1e-12,
            max(abs(colSums(null[,,j])-1))<1e-12)
  data.table(position=rep(p$position,500), draw=rep(1:500,each=100),
             fraction=as.vector(null[,,j]), weighting=weight_order[j])
}))
panel_levels <- c('intron_weighted', 'site_weighted')
panel_names <- c('Equal intron weight', 'Site weighted')
profiles[, panel:=factor(weighting,levels=panel_levels,labels=panel_names)]
cloud[, panel:=factor(weighting,levels=panel_levels,labels=panel_names)]

# Verify saved percentile columns against the exact cached draws.
for (j in seq_along(weight_order)) {
  p <- profiles[weighting==weight_order[j]][order(position)]
  for (q in c(10,25,75,90)) {
    stopifnot(max(abs(p[[paste0('null_p',q)]] -
      apply(null[,,j],1,quantile,q/100))) < 1e-12)
  }
}
fig <- ggplot() +
  geom_ribbon(data=profiles,aes(position,ymin=null_p10,ymax=null_p90,
    fill='10th-90th percentiles')) +
  geom_ribbon(data=profiles,aes(position,ymin=null_p25,ymax=null_p75,
    fill='25th-75th percentiles')) +
  scale_fill_manual(name='Null bands',
    breaks=c('25th-75th percentiles','10th-90th percentiles'),
    values=c('25th-75th percentiles'='#7CAACD','10th-90th percentiles'='#D6EBF5')) +
  geom_line(data=cloud,aes(position,fraction,group=draw,color='Randomized null (500 draws)'),
            linewidth=.18,alpha=.035) +
  geom_line(data=profiles,aes(position,background_median,color='Null median'),linewidth=.55) +
  geom_line(data=profiles,aes(position,observed,color='Observed MSI2'),linewidth=.8) +
  facet_wrap(~panel,nrow=1) +
  scale_color_manual(name=NULL,
    breaks=c('Observed MSI2','Randomized null (500 draws)','Null median'),
    values=c('Observed MSI2'='#D62728','Randomized null (500 draws)'='#46B8D5',
             'Null median'='black')) +
  guides(color=guide_legend(order=1,override.aes=list(alpha=1,linewidth=c(.8,.5,.55))),
         fill=guide_legend(order=2)) +
  scale_x_continuous(breaks=c(0,25,50,75,100),limits=c(0,100),expand=expansion(mult=.015)) +
  scale_y_continuous(labels=scales::label_number(accuracy=.001)) +
  labs(title='Intron summit positions versus PARPIPE-style background',
       subtitle='20,364 summits in 9,565 introns >500 nt; original assignments retained',
       x='Intron position (%)', y='Weighted fraction per 1% bin',
       caption=paste('All 500 cached null draws; raw 1% bins, no smoothing. Seed: 260922.',
         'Left: each represented intron has total weight one. Right: each summit has equal weight.',sep='\n')) +
  theme_classic(base_size=12) +
  theme(strip.background=element_rect(fill='#F2F2F2',color='#B0B0B0',linewidth=.4),
        strip.text=element_text(size=12), legend.position='bottom', legend.box='vertical',
        legend.key.width=grid::unit(1.4,'cm'),
        plot.title=element_text(size=16),plot.subtitle=element_text(size=11),
        plot.caption=element_text(size=10,hjust=0),
        panel.spacing=grid::unit(1.0,'cm'),plot.margin=margin(12,14,12,12))
stem <- file.path(out,'intron_summit_positions_PARpipe_style_cloud_percentiles')
ggsave(paste0(stem,'.png'),fig,width=12,height=6.8,dpi=300,bg='white')
ggsave(paste0(stem,'.pdf'),fig,width=12,height=6.8,bg='white',useDingbats=FALSE)

stopifnot(identical(hash_before,tools::md5sum(sources)))
provenance <- data.table(
  field=c('observed_table','observed_table_md5','null_cache','null_cache_md5',
          'original_analysis_script','original_analysis_script_md5','plotting_script',
          'plotting_script_md5','null_array_axes','null_slice_1','null_slice_2',
          'null_draws_available','null_draws_displayed','smoothing','randomizations_regenerated',
          'original_seed','length_filter','summits','introns','null_line_alpha',
          'source_files_unchanged','validation','R_version','ggplot2_version'),
  value=c(observed_path,unname(hash_before[1]),null_path,unname(hash_before[2]),
          analysis_script,unname(hash_before[3]),script_path,unname(tools::md5sum(script_path)),
          '100 bins x 500 draws x 2 weighting schemes','site_weighted','intron_weighted',
          '500','500','None; raw 1% bins connected with straight segments','No',
          '260922','Strictly >500 nt','20364','9565','0.035','Yes; MD5 before/after identical',
          'Cached medians and 10/25/75/90 percentiles match observed table; observed and all null profiles sum to one (tolerance 1e-12)',
          R.version.string,as.character(packageVersion('ggplot2'))))
fwrite(provenance,file.path(out,'plot_provenance_cloud_percentiles.tsv'),sep='\t')
caption <- paste('Reference MSI2 summit positions within introns longer than 500 nt, shown in',
  'transcript-oriented normalized 0-100% intron coordinates. Red: observed MSI2 profile.',
  'Transparent blue curves: all 500 cached randomized profiles from the existing PARPIPE-style',
  'background procedure. Light blue band: 10th-90th percentiles; darker blue: 25th-75th percentiles.',
  'Black line: null median. No gray outer envelope is drawn. Left, each represented intron contributes',
  'total weight one; right, each summit contributes equally. Profiles use unsmoothed 1% bins.',
  'The observed excess toward the intron beginning is consistent with donor-proximal enrichment',
  'across multiple relative-position bins; these bins do not specify a fixed nucleotide distance',
  'from the splice site. This is a restyle of the existing analysis, with no new randomizations.')
writeLines(c('# Intron summit combined null-cloud and percentile restyle','',
  '## Figure caption','',caption,'',
  '## Reproduction and provenance','',
  paste0('Run from the project root: `Rscript ',script_path,'`.'),'',
  paste0('- Observed profiles: `',observed_path,'`.'),
  paste0('- Cached null array: `',null_path,'`.'),
  paste0('- Original generating code: `',analysis_script,'`.'),
  '- All 500 draws are displayed in each panel, with no subset selection or smoothing.',
  '- Existing tables, caches, analytical figures, assignments, >500 nt filter, weighting, and numerical results are unchanged.',
  '- Blue overlap shows the cached null cloud; it is not a p-value color scale.',
  '- The null remains the prior independent half-integer positional-grid adaptation with 12.5 nt end margins, not a new implementation of the full historical PARPIPE script.',
  '- See plot_provenance_cloud_percentiles.tsv for source MD5 hashes, array-slice mapping, seed, software versions, and validation.',
  '', '## Checks','',
  'Verified 100 bins and 500 draws per weighting scheme; cached medians and 10/25/75/90 percentiles match the saved profile table; every observed/null profile sums to one; source MD5 hashes are identical before and after plotting.'),
  file.path(out,'README_cloud_percentiles.md'))
cat('Restyle saved; all 500 cached draws used; source hashes unchanged.\n')
