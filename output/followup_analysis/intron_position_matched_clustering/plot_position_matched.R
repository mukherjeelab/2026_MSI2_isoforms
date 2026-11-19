suppressPackageStartupMessages({library(data.table);library(ggplot2)})
out <- 'output/followup_analysis/intron_position_matched_clustering'
models <- c('original_uniform','position_matched','first20_uniform')
labels <- c('Uniform, full introns','Position-matched, full introns','Uniform, first 20% only')
cache_paths <- c('output/followup_analysis/site_clustering/cache/intron_site_distances.rds',
  file.path(out,'cache',paste0(models[-1],'_distances.rds')))
before <- tools::md5sum(cache_paths)
stats <- fread(file.path(out,'uniform_vs_position_matched_vs_first20.tsv'))
stats[,model_label:=factor(model,levels=models,labels=labels)]
stats[,stratum:=factor(stratum,levels=c('All','2 sites','3-4 sites','>=5 sites'))]
save_plot <- function(p,name,w=13,h=6) {
  ggsave(file.path(out,paste0(name,'.png')),p,width=w,height=h,dpi=250,bg='white')
  ggsave(file.path(out,paste0(name,'.pdf')),p,width=w,height=h,bg='white',useDingbats=FALSE)
}
t <- stats[metric!='median_distance_nt']
t[,threshold_nt:=as.numeric(sub('within_([0-9]+)_nt','\\1',metric))]
fwrite(t,file.path(out,'threshold_plot_inputs.tsv'),sep='\t')
save_plot(ggplot(t,aes(threshold_nt))+
  geom_ribbon(aes(ymin=null_low_95,ymax=null_high_95),fill='#CBE5F0')+
  geom_line(aes(y=null_reference,color='Null expectation'),linewidth=.5)+
  geom_point(aes(y=null_reference,color='Null expectation'),size=1.3)+
  geom_line(aes(y=observed,color='Observed MSI2'),linewidth=.7)+
  geom_point(aes(y=observed,color='Observed MSI2'),size=1.5)+
  facet_grid(model_label~stratum)+
  scale_color_manual(name=NULL,values=c('Observed MSI2'='#D62728','Null expectation'='#18394B'))+
  scale_x_continuous(breaks=c(50,100,250,500))+
  scale_y_continuous(labels=scales::label_percent(),limits=c(0,1))+
  theme_bw(base_size=10)+theme(legend.position='bottom',panel.grid.minor=element_blank(),plot.caption=element_text(hjust=0))+
  labs(title='Intron clustering after accounting for relative-position preference',
    subtitle='Full introns: 3,799 features / 14,598 summits. First 20%: 1,089 features / 3,475 summits.',
    x='Neighbor-distance threshold (nt)',y='Fraction of summits with a neighbor at or below threshold',
    caption=paste('Blue band: 95% null interval from 500 draws. Each exact feature and its analyzed site count are preserved.',
      'First-20% strata count sites inside that zone, and nearest neighbors must also lie inside the zone.',sep='\n')),
  'intron_position_matched_thresholds',14,9)

curves <- list()
for(j in seq_along(models)) {
  d <- readRDS(cache_paths[j])
  kde <- function(x) {
    q <- density(c(-x,x),bw=25,from=0,to=1000,n=2048)
    data.table(distance_nt=seq(0,1000,2),density=pmax(0,2*approx(q$x,q$y,xout=seq(0,1000,2))$y))
  }
  obs <- kde(d$observed_nearest);obs[,`:=`(curve='observed',iteration=0L)]
  nul <- rbindlist(lapply(1:500,function(it) {
    r<-kde(d$nearest_neighbor_null[,it]);r[,`:=`(curve='null',iteration=it)];r
  }))
  med <- nul[,.(density=median(density)),by=distance_nt];med[,`:=`(curve='null_median',iteration=0L)]
  curves[[j]] <- rbindlist(list(obs,nul,med),use.names=TRUE)[,model:=models[j]]
}
curves <- rbindlist(curves)
fwrite(curves,file.path(out,'density_zoom_plot_inputs.tsv.gz'),sep='\t')
curves[,model_label:=factor(model,levels=models,labels=labels)]
save_plot(ggplot()+
  geom_line(data=curves[curve=='null'],aes(distance_nt,density,group=iteration,
    color='Null (500 draws)'),alpha=.035,linewidth=.17)+
  geom_line(data=curves[curve=='null_median'],aes(distance_nt,density,color='Null median'),linewidth=.5)+
  geom_line(data=curves[curve=='observed'],aes(distance_nt,density,color='Observed MSI2'),linewidth=.75)+
  facet_wrap(~model_label,nrow=1)+
  scale_color_manual(name=NULL,breaks=c('Observed MSI2','Null (500 draws)','Null median'),
    values=c('Observed MSI2'='#D62728','Null (500 draws)'='#46B8D5','Null median'='black'))+
  guides(color=guide_legend(override.aes=list(alpha=1,linewidth=c(.75,.5,.5))))+
  theme_classic(base_size=11)+theme(legend.position='bottom',strip.background=element_rect(fill='#F2F2F2',color='grey70'),
    plot.caption=element_text(hjust=0))+
  labs(title='MSI2 intron nearest-neighbor distances: uniform and position-matched nulls',
    subtitle='Full introns: 14,598 summits in 3,799 introns. First 20%: 3,475 summits in 1,089 introns.',
    x='Nearest-neighbor summit distance (nt)',y='Density (per nt)',
    caption=paste('All introns >500 nt. Shared axes; densities retain all tail mass beyond the 0-1,000 nt view.',
      'Display only: reflected Gaussian KDE, bandwidth 25 nt. Tests use raw distances. First-20% panel uses a different subset.',sep='\n')),
  'intron_position_matched_density_zoom',15,6)

p <- fread(file.path(out,'positional_null_diagnostics.tsv'))[model=='position_matched']
v <- fread(file.path(out,'position_matching_validation.tsv'))[model=='position_matched']
p <- merge(p,v[,.(bin,analytic_expected_probability)],by='bin')
save_plot(ggplot(p,aes(position_percent))+
  geom_ribbon(aes(ymin=null_low,ymax=null_high),fill='#D6EBF5')+
  geom_line(aes(y=global_empirical_probability,color='Global empirical preference'),linewidth=.7)+
  geom_line(aes(y=null_mean,color='Simulated mean'),linewidth=.55)+
  geom_line(aes(y=analytic_expected_probability,color='Expected after grid/margins'),linetype=2,linewidth=.4)+
  scale_color_manual(name=NULL,values=c('Global empirical preference'='#D62728',
    'Simulated mean'='#2487A5','Expected after grid/margins'='black'))+
  theme_classic(base_size=11)+theme(legend.position='bottom',plot.caption=element_text(hjust=0))+
  labs(title='Position-matched null retains the global intron-position preference',
    subtitle='Preference fitted to all 20,364 reference summits in 9,565 eligible introns >500 nt',
    x='Intron position (%)',y='Fraction per 1% bin',
    caption=paste('Blue shading: 95% interval across 500 position-matched draws for the multi-site intron cohort.',
      'Small end-bin differences reflect the original half-integer grid and 12.5 nt margin conditioning.',sep='\n')),
  'intron_position_preference_check',11,5.5)
stopifnot(identical(before,tools::md5sum(cache_paths)))
cat('Plots saved; null caches unchanged.\n')
