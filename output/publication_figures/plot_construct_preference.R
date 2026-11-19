suppressPackageStartupMessages({library(data.table); library(ggplot2)})
# Run from project root. Use the exact retained peak set in the reference figure.
input <- 'data/results/peak_analysis/figures/all_nonzero_peaks/transcript_region_bias/relaxed_filter_continuous_delta_delta_peaks.tsv'
out <- 'output/publication_figures'
d <- fread(input)
regions <- c("3' UTR", 'CDS', 'Intron', 'Other')
d[, original_region := transcript_region]
d[transcript_region == "5' UTR", transcript_region := 'Other']
d[, transcript_region := factor(transcript_region, levels=regions)]
stopifnot(nrow(d)==9369, all(is.finite(d$delta_delta_324_minus_328)),
          identical(as.integer(table(d$transcript_region)),c(8271L,398L,290L,410L)))
# Monte Carlo conditional KS test preserves ties in the pooled peak values.
# Each region is compared to 3' UTR, not to zero. Peaks are the observations.
set.seed(20260928)
reference <- d[transcript_region == "3' UTR",delta_delta_324_minus_328]
stats <- rbindlist(lapply(setdiff(regions,"3' UTR"),function(region) {
  x <- d[transcript_region == region,delta_delta_324_minus_328]
  test <- ks.test(x,reference,alternative='two.sided',exact=FALSE,
                  simulate.p.value=TRUE,B=99999)
  data.table(transcript_region=region,reference="3' UTR",n=length(x),n_reference=length(reference),
             KS_D=unname(test$statistic),p_value=test$p.value,permutations=99999L)
}))
stats[,q_BH := p.adjust(p_value,method='BH')]
fwrite(stats,file.path(out,'construct_preference_vs_3UTR_statistics.tsv'),sep='\t')
fwrite(d[transcript_region=='Other',.N,by=transcript_region_detail],
       file.path(out,'construct_preference_other_annotations.tsv'),sep='\t')
fwrite(d[,.(n=.N,median=median(delta_delta_324_minus_328)),by=transcript_region],
       file.path(out,'construct_preference_region_summary.tsv'),sep='\t')
stats[, transcript_region := factor(transcript_region, levels=regions)]
stats[, label := vapply(q_BH, function(v) {
  if(v < .001) sprintf('q = %.1e',v) else sprintf('q = %.3f',v)
}, character(1))]
p <- ggplot(d,aes(transcript_region,delta_delta_324_minus_328)) +
  geom_hline(yintercept=0,color='grey45',linewidth=.4) +
  geom_violin(fill='#BFBFBF',color='black',linewidth=.4,scale='width',trim=TRUE) +
  geom_boxplot(width=.14,outlier.shape=NA,fill='white',color='black',linewidth=.4) +
  geom_text(data=stats,aes(x=transcript_region,y=13.2,label=label),
            inherit.aes=FALSE,size=2.35,color='black') +
  annotate('text',x=1,y=13.2,label='Reference',size=2.35) +
  scale_y_continuous(breaks=seq(-10,10,5),minor_breaks=seq(-12.5,12.5,2.5),
                     expand=expansion(mult=c(.025,.025))) +
  labs(x=NULL,y='LFC difference (MSI2-324 - MSI2-328)') +
  theme_minimal(base_size=10,base_family='Helvetica') +
  theme(panel.grid.major=element_line(color='#E5E5E5',linewidth=.45),
        panel.grid.minor.y=element_line(color='#EEEEEE',linewidth=.3),
        panel.grid.minor.x=element_blank(),
        axis.text=element_text(color='black'),
        axis.title.y=element_text(margin=margin(r=6)),
        plot.margin=margin(7,7,5,5))
stem <- file.path(out,'construct_preference_by_transcript_region')
ggsave(paste0(stem,'.png'),p,width=3.8,height=3.5,dpi=600,bg='white')
ggsave(paste0(stem,'.pdf'),p,width=3.8,height=3.5,device='pdf',useDingbats=FALSE,bg='white')
fwrite(d,paste0(stem,'_source_data.tsv'),sep='\t')
