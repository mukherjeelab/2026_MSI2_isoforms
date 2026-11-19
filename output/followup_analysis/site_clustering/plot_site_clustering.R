suppressPackageStartupMessages({library(data.table);library(ggplot2)})
out <- 'output/followup_analysis/site_clustering'
cache_paths <- file.path(out,'cache',paste0(c('intron','3UTR'),'_site_distances.rds'))
hash_before <- tools::md5sum(cache_paths)
caches <- lapply(cache_paths,readRDS)
names(caches) <- c('intron','3UTR')
stats <- fread(file.path(out,'MSI2_site_clustering_all_statistics.tsv'))
save_plot <- function(p,name,w=12,h=6.5) {
  ggsave(file.path(out,paste0(name,'.png')),p,width=w,height=h,dpi=250,bg='white')
  ggsave(file.path(out,paste0(name,'.pdf')),p,width=w,height=h,bg='white',useDingbats=FALSE)
}
# A shared broad view covers >=99% of observed and null distances in each type.
# Densities remain normalized to ALL distances; tails are never renormalized away.
full_limit <- ceiling(max(vapply(caches,function(d) max(
  quantile(d$observed_nearest,.99),quantile(d$nearest_neighbor_null,.99)),numeric(1)))/1000)*1000
bw <- 25
kde <- function(x,limit,nout) {
  # Reflection boundary correction on nonnegative distances, shared bandwidth for all curves.
  # Fine FFT grid <=5 nt, then interpolate only for the displayed line coordinates.
  grid_n <- 2^ceiling(log2(max(2048,limit/5+1)))
  d <- density(c(-x,x),bw=bw,from=0,to=limit,n=grid_n)
  xx <- seq(0,limit,length.out=nout)
  data.table(distance_nt=xx,density=pmax(0,2*approx(d$x,d$y,xout=xx)$y))
}
tail_rows <- list()
for(view in c('broad','zoom')) {
  limit <- if(view=='broad') full_limit else 1000
  nout <- if(view=='broad') 2048 else 501
  rows <- list()
  for(kind_name in names(caches)) {
    d <- caches[[kind_name]]
    obs <- kde(d$observed_nearest,limit,nout)
    obs[,`:=`(kind=kind_name,curve='observed',iteration=0L)]
    nul <- rbindlist(lapply(1:500,function(it) {
      v <- kde(d$nearest_neighbor_null[,it],limit,nout)
      v[,`:=`(kind=kind_name,curve='null',iteration=it)]
      v
    }))
    med <- nul[,.(density=median(density)),by=.(kind,distance_nt)]
    med[,`:=`(curve='null_median',iteration=0L)]
    rows[[kind_name]] <- rbindlist(list(obs,nul,med),use.names=TRUE)
    tail_rows[[paste(view,kind_name)]] <- data.table(view=view,kind=kind_name,x_max_nt=limit,
      observed_fraction_beyond=mean(d$observed_nearest>limit),
      null_fraction_beyond=mean(d$nearest_neighbor_null>limit),bandwidth_nt=bw)
  }
  curves <- rbindlist(rows)
  fwrite(curves,file.path(out,paste0('MSI2_site_distance_density_inputs_',view,'.tsv.gz')),sep='\t')
  curves[,panel:=factor(kind,levels=c('intron','3UTR'),labels=c(
    'Introns >500 nt | 3,799 features; 14,598 summits',
    "3'UTRs >=100 nt | 5,317 features; 22,952 summits"))]
  p <- ggplot() +
    geom_line(data=curves[curve=='null'],aes(distance_nt,density,group=iteration,
      color='Randomized null (500 draws)'),alpha=.035,linewidth=.17) +
    geom_line(data=curves[curve=='null_median'],aes(distance_nt,density,color='Null median'),linewidth=.5) +
    geom_line(data=curves[curve=='observed'],aes(distance_nt,density,color='Observed MSI2'),linewidth=.8) +
    facet_wrap(~panel,nrow=1,scales='free_y') +
    scale_color_manual(name=NULL,breaks=c('Observed MSI2','Randomized null (500 draws)','Null median'),
      values=c('Observed MSI2'='#D62728','Randomized null (500 draws)'='#46B8D5','Null median'='black')) +
    guides(color=guide_legend(override.aes=list(alpha=1,linewidth=c(.8,.5,.5)))) +
    scale_x_continuous(labels=scales::label_number(big.mark=','),expand=expansion(mult=c(0,.015))) +
    theme_classic(base_size=11) + theme(legend.position='bottom',
      strip.background=element_rect(fill='#F2F2F2',color='grey70'),plot.caption=element_text(hjust=0)) +
    labs(title='MSI2 nearest-neighbor summit distances within features',
      subtitle=sprintf('Same feature assignments and site counts; 500 within-feature randomizations | 0-%s nt view',
                        format(limit,big.mark=',')),
      x='Nearest-neighbor summit distance (nt)',y='Density (per nt)',
      caption=paste('One distance per summit; only features with >=2 sites. Panel y scales differ.',
        'Display only: reflected Gaussian density, bandwidth 25 nt for every curve. Numerical tests use raw distances.',
        'Same x limits in both panels; tails outside the view are retained in density normalization.',sep='\n'))
  save_plot(p,if(view=='broad') 'MSI2_site_distance_density' else 'MSI2_site_distance_density_zoom')
  cat(view,'density plot complete; upper x limit',limit,'nt\n')
}
fwrite(rbindlist(tail_rows),file.path(out,'MSI2_density_display_parameters.tsv'),sep='\t')
t <- stats[distance_type=='nearest_neighbor' & metric!='median_distance_nt']
t[,threshold_nt:=as.numeric(sub('within_([0-9]+)_nt','\\1',metric))]
t[,panel:=factor(kind,levels=c('intron','3UTR'),labels=c('Introns >500 nt',"3'UTRs >=100 nt"))]
t[,stratum:=factor(stratum,levels=c('All','2 sites','3-4 sites','>=5 sites'))]
fwrite(t,file.path(out,'MSI2_site_clustering_threshold_plot_inputs.tsv'),sep='\t')
save_plot(ggplot(t,aes(threshold_nt)) +
  geom_ribbon(aes(ymin=null_low_95,ymax=null_high_95),fill='#CBE5F0') +
  geom_line(aes(y=null_reference,color='Null expectation'),linewidth=.55) +
  geom_point(aes(y=null_reference,color='Null expectation'),size=1.4) +
  geom_line(aes(y=observed,color='Observed MSI2'),linewidth=.7) +
  geom_point(aes(y=observed,color='Observed MSI2'),size=1.7) +
  facet_grid(panel~stratum) +
  scale_color_manual(name=NULL,values=c('Observed MSI2'='#D62728','Null expectation'='#18394B')) +
  scale_x_continuous(breaks=c(50,100,250,500)) +
  scale_y_continuous(labels=scales::label_percent(),limits=c(0,1)) +
  theme_bw(base_size=11) + theme(legend.position='bottom',panel.grid.minor=element_blank(),
                              plot.caption=element_text(hjust=0)) +
  labs(title='MSI2 local clustering: observed fractions versus feature-matched null',
    subtitle='Each summit contributes once; site-count strata are preserved in every randomization',
    x='Neighbor-distance threshold (nt)',y='Fraction of summits with a neighbor at or below threshold',
    caption='Dark line: mean across 500 null iterations. Blue band: pointwise 95% null interval. Red: observed.'),
  'MSI2_site_clustering_thresholds',14,7)
stopifnot(identical(hash_before,tools::md5sum(cache_paths)))
cat('All plotting complete; distance caches unchanged.\n')
