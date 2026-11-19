suppressPackageStartupMessages({library(data.table);library(ggplot2)})
setDTthreads(2)
out <- 'output/publication_figures'
ref <- fread('output/followup_analysis/intron_gt500/intron_gt500_profiles.tsv')[weighting=='intron_weighted'][order(position)]
rn <- readRDS('output/followup_analysis/intron_gt500/PARpipe_intron_gt500_nulls.rds')[,,2]
stopifnot(identical(dim(rn),c(100L,500L)),nrow(ref)==100,
          max(abs(ref$background_median-apply(rn,1,median)))<1e-12,
          abs(sum(ref$observed)-1)<1e-12)
for(q in c(10,25,75,90)) stopifnot(max(abs(ref[[paste0('null_p',q)]]-
  apply(rn,1,quantile,q/100)))<1e-12)
cloud <- data.table(position=rep(ref$position,500),draw=rep(1:500,each=100),value=as.vector(rn))
p <- ggplot() +
  geom_ribbon(data=ref,aes(position,ymin=null_p10,ymax=null_p90,fill='10th-90th')) +
  geom_ribbon(data=ref,aes(position,ymin=null_p25,ymax=null_p75,fill='25th-75th')) +
  geom_line(data=cloud,aes(position,value,group=draw),color='#46B8D5',alpha=.035,linewidth=.15) +
  geom_line(data=ref,aes(position,background_median),color='black',linewidth=.4) +
  geom_line(data=ref,aes(position,observed),color='#D62728',linewidth=.65) +
  scale_fill_manual(name='Null percentiles',values=c('10th-90th'='#D6EBF5','25th-75th'='#7CAACD')) +
  scale_x_continuous(breaks=c(0,25,50,75,100),limits=c(0,100),expand=expansion(mult=.015)) +
  scale_y_continuous(breaks=seq(.005,.020,.005),labels=scales::label_number(accuracy=.001)) +
  coord_cartesian(ylim=c(.005,.0235),expand=FALSE) +
  labs(x='Intron position (%)',y='Weighted fraction per 1% bin') +
  theme_minimal(base_size=10,base_family='Helvetica') +
  theme(panel.grid.major=element_line(color='#E5E5E5',linewidth=.4),
        panel.grid.minor=element_blank(),axis.text=element_text(color='black'),
        axis.title.y=element_text(margin=margin(r=6)),
        legend.position='bottom',legend.text=element_text(size=8),legend.title=element_text(size=8),
        legend.key.size=grid::unit(3,'mm'),legend.spacing.x=grid::unit(1.2,'mm'),
        legend.margin=margin(2,0,0,0),plot.margin=margin(5,7,5,5))
stem <- file.path(out,'endogenous_intron_equal_weight')
ggsave(paste0(stem,'.png'),p,width=3.8,height=3.2,dpi=600,bg='white')
ggsave(paste0(stem,'.pdf'),p,width=3.8,height=3.2,device='pdf',useDingbats=FALSE,bg='white')
fwrite(ref,paste0(stem,'_source_data.tsv'),sep='\t')
