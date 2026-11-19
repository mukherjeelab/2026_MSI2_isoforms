#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ragg)})
input <- 'data/results_round2'
out <- 'output/round2_comparisons'
dir.create(file.path(out, 'figures'), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(out, 'tables'), recursive=TRUE, showWarnings=FALSE)
write_tab <- function(x, name) fwrite(x, file.path(out,'tables',paste0(name,'.tsv')), sep='\t')
save_plot <- function(p, name, w=12, h=12) {
  ggsave(file.path(out,'figures',paste0(name,'.png')),p,width=w,height=h,dpi=180,device=ragg::agg_png,bg='white')
  ggsave(file.path(out,'figures',paste0(name,'.pdf')),p,width=w,height=h,device=grDevices::pdf,bg='white')
}
theme_set(theme_bw(base_size=11)+theme(legend.position='bottom',panel.grid.minor=element_blank(),
  strip.background=element_rect(fill='#EDF1F4'),plot.title=element_text(face='bold'),plot.caption=element_text(hjust=0)))
modes <- c('genome_full','genome_equal','peaks_full','peaks_after_genome_equal','peaks_equal')
mode_labels <- c('Full genome','Equal genome','Full peaks','Peaks after equal genome','Equal peaks')
decorate <- function(d) {
  d[, motif:=factor(sub('__.*','',variant),levels=c('all','uag'),labels=c('All reads','UAG reads'))]
  d[, mode:=factor(sub('.*__','',variant),levels=modes,labels=mode_labels)]
  d
}
meta <- fread(file.path(input,'sample_metadata.tsv'))
read_all <- function(name) rbindlist(lapply(meta$sample_id,function(s) fread(file.path(input,s,name))))
cats <- read_all('category_counts.tsv')
libs <- unique(cats[,.(sample_id,variant,total_selected_reads)])
stopifnot(nrow(libs)==80,all(libs$total_selected_reads>0))
stopifnot(all(cats[,.(ok=sum(assigned_reads)==unique(total_selected_reads)),by=.(sample_id,variant)]$ok))
libs <- decorate(merge(libs,meta,by='sample_id'))
write_tab(libs,'library_sizes')
profiles <- merge(read_all('metaprofiles.tsv'),meta,by='sample_id')
stopifnot(nrow(profiles)==56000, !anyDuplicated(profiles[,.(sample_id,variant,kind,bin_index)]))
stopifnot(all(profiles[,.(ok=uniqueN(features_with_nonempty_bin)==1 && uniqueN(total_bin_nt)==1),by=.(kind,bin_index)]$ok))
profiles[,x:=(bin_start+bin_end)/2]
profiles[,group:=paste(construct,assay)]
profiles <- decorate(profiles)
write_tab(profiles,'metaprofiles_all_samples')
colors <- c('324 CLIP'='#008C76','324 SMI'='#64B8A6','328 CLIP'='#CE531C','328 SMI'='#DE9E6F')
titles <- c(splice_donor='Splice donor: expanded window',splice_acceptor='Splice acceptor: expanded window',
  transcript_end="Annotated transcript 3' end: expanded window",intron='Intron: standardized length',`3UTR`="Spliced 3' UTR: standardized length")
means <- profiles[,.(depth=mean(mean_depth_per_feature_per_million_reads),raw_depth=mean(mean_depth_per_feature)),
  by=.(variant,motif,mode,kind,bin_index,x,group,construct,assay)]
write_tab(means,'metaprofile_group_means')
for (k in names(titles)) {
  d <- profiles[kind==k]; a <- means[kind==k]
  scaled <- k %in% c('intron','3UTR')
  xlabel <- if(scaled) "Position along feature (5' to 3'; %)" else 'Position relative to anchor (nt; transcript orientation)'
  for (metric in c('normalized','raw')) {
    d[,value:=if(metric=='normalized') mean_depth_per_feature_per_million_reads else mean_depth_per_feature]
    a[,value:=if(metric=='normalized') depth else raw_depth]
    p <- ggplot(d,aes(x,value,color=group,linetype=assay))+
      geom_line(aes(group=sample_id),alpha=.3,linewidth=.35)+
      geom_line(data=a,aes(group=group),linewidth=.8)+
      facet_grid(mode~motif)+scale_color_manual(values=colors)+scale_linetype_manual(values=c(CLIP='solid',SMI='dashed'))+
      labs(title=titles[[k]],subtitle='Thin lines: biological replicates A and B; thick lines: their arithmetic mean',
        x=xlabel,y=if(metric=='normalized') 'Mean feature depth per million selected reads' else 'Mean feature depth (raw)',color=NULL,linetype=NULL,
        caption='Fixed annotation universe; zero-coverage features included. Axes shared across all ten variants.\nPeak-selected reads contribute their full aligned blocks; profiles are not clipped to peak intervals.')
    if(!scaled) p <- p+geom_vline(xintercept=0,linewidth=.35,color='grey45')
    save_plot(p,paste0(k,'_',metric))
  }
}
# Descriptive paired differences avoid a depth-scale pseudocount in metaprofiles.
paired <- dcast(profiles,variant+motif+mode+kind+bin_index+x+construct+biological_replicate~assay,
  value.var='mean_depth_per_feature_per_million_reads')
paired[,enrichment:=CLIP-SMI]
dd <- dcast(paired,variant+motif+mode+kind+bin_index+x+biological_replicate~construct,value.var='enrichment')
dd[,delta:=`324`-`328`]
ddmean <- dd[,.(delta=mean(delta)),by=.(variant,motif,mode,kind,bin_index,x)]
write_tab(paired,'paired_metaprofile_clip_minus_smi')
write_tab(dd,'paired_metaprofile_construct_differences')
for(k in names(titles)) {
  p <- ggplot(dd[kind==k],aes(x,delta))+geom_hline(yintercept=0,color='grey60')+
    geom_line(aes(group=biological_replicate),alpha=.35,color='#543B87')+
    geom_line(data=ddmean[kind==k],color='#543B87',linewidth=.8)+facet_grid(mode~motif)+
    labs(title=paste(titles[[k]],'- construct difference'),subtitle='(324 CLIP - SMI) - (328 CLIP - SMI); positive = greater additive enrichment in 324',
      x=if(k %in% c('intron','3UTR')) "Position along feature (%)" else 'Position relative to anchor (nt)',
      y='Difference in mean feature depth per million selected reads',caption='Thin lines: replicate labels A/B across constructs; thick line: mean. Descriptive, without significance testing.')
  save_plot(p,paste0(k,'_construct_difference'))
}
# One fixed, construct-independent peak universe across all variants.
counts <- read_all('peak_counts.tsv.gz')
stopifnot(!anyDuplicated(counts[,.(sample_id,variant,peak_id)]))
universe <- counts[variant=='all__genome_full',.(total=sum(overlapping_reads)),by=peak_id][total>0,peak_id]
counts <- merge(counts[peak_id %in% universe],libs[,.(sample_id,variant,total_selected_reads,construct,assay,biological_replicate)],by=c('sample_id','variant'))
stopifnot(nrow(counts)==length(universe)*80)
counts[,cpm:=overlapping_reads/total_selected_reads*1e6]
avg <- counts[,.(cpm=mean(cpm)),by=.(variant,peak_id,construct,assay)]
wide <- dcast(avg,variant+peak_id~construct+assay,value.var='cpm')
wide[,lfc324:=log2((`324_CLIP`+.5)/(`324_SMI`+.5))]
wide[,lfc328:=log2((`328_CLIP`+.5)/(`328_SMI`+.5))]
wide[,delta:=lfc324-lfc328]
wide[,delta_pc01:=log2((`324_CLIP`+.1)/(`324_SMI`+.1))-log2((`328_CLIP`+.1)/(`328_SMI`+.1))]
support <- counts[,.(supported=any(overlapping_reads>0)),by=.(variant,peak_id)]
ann <- fread(file.path(input,'reference_peaks','peak_annotations.tsv.gz'))
wide <- decorate(merge(merge(wide,support,by=c('variant','peak_id')),ann,by='peak_id'))
wide[,region:=factor(assigned_class,levels=c('5UTR','CDS','3UTR','intron','other'),labels=c("5' UTR",'CDS',"3' UTR",'Intron','Other'))]
stopifnot(!anyNA(wide$region),all(is.finite(wide$delta)))
write_tab(wide,'continuous_peak_preference')
region_colors <- c("5' UTR"='#4C78A8',CDS='#F58518',"3' UTR"='#54A24B',Intron='#B279A2',Other='#79706E')
for(supported_only in c(FALSE,TRUE)) {
  d <- if(supported_only) wide[supported==TRUE] else wide
  p <- ggplot(d,aes(region,delta,fill=region))+geom_hline(yintercept=0,color='#9E2A2B',linewidth=.4)+
    geom_violin(scale='width',linewidth=.25)+geom_boxplot(width=.13,fill='white',outlier.shape=NA,linewidth=.3)+
    facet_grid(mode~motif)+scale_fill_manual(values=region_colors,guide='none')+coord_flip()+
    labs(title='Continuous 324-versus-328 peak preference',subtitle=paste0(if(supported_only) 'Sensitivity: peaks with any counts within each variant' else paste('Fixed universe:',format(length(universe),big.mark=','),'peaks with any counts in all-read full genome')),
      x=NULL,y='Delta-delta log2 enrichment (positive favors 324)',caption='log2[(mean CLIP CPM + 0.5)/(mean SMI CPM + 0.5)] in 324 minus 328. CPM uses selected-library sizes.\nSame axes throughout. Fixed-universe plots retain zero-support peaks at zero; sensitivity plots exclude them.')
  save_plot(p,if(supported_only) 'continuous_preference_supported_peaks' else 'continuous_preference_all_variants',12,14)
}
summary <- wide[,.(n_peaks=.N,n_supported=sum(supported),median_delta=median(delta),q25=quantile(delta,.25),q75=quantile(delta,.75),
  fraction_positive=mean(delta>0),median_delta_pc01=median(delta_pc01)),by=.(variant,motif,mode,region)]
write_tab(summary,'continuous_region_summary')
write_tab(wide[supported==TRUE,.(n_peaks=.N,median_delta=median(delta)),by=.(variant,region)],'continuous_supported_region_summary')
# Genome/peak selection pairs must retain identical per-peak counts: only the denominator changes.
count_check <- dcast(counts,sample_id+peak_id~variant,value.var='overlapping_reads')
stopifnot(all(count_check$all__genome_full==count_check$all__peaks_full),
  all(count_check$uag__genome_full==count_check$uag__peaks_full),
  all(count_check$all__genome_equal==count_check$all__peaks_after_genome_equal),
  all(count_check$uag__genome_equal==count_check$uag__peaks_after_genome_equal))
profile_base <- means[mode=='Full genome',.(motif,kind,bin_index,group,baseline_depth=depth)]
profile_change <- merge(means,profile_base,by=c('motif','kind','bin_index','group'))
profile_change <- profile_change[,.(pearson=cor(depth,baseline_depth),
  mean_depth_ratio=mean(depth)/mean(baseline_depth),
  relative_rms_change=sqrt(mean((depth-baseline_depth)^2))/sqrt(mean(baseline_depth^2))),by=.(variant,motif,mode,kind,group)]
write_tab(profile_change,'metaprofile_changes_from_matching_full_genome')
base <- wide[variant=='all__genome_full',.(peak_id,baseline_delta=delta)]
comparison <- merge(wide,base,by='peak_id')
changes <- comparison[,.(spearman=cor(delta,baseline_delta,method='spearman'),median_absolute_change=median(abs(delta-baseline_delta)),
  fraction_sign_changed=mean(sign(delta)!=sign(baseline_delta)),n_unsupported=sum(!supported)),by=.(variant,motif,mode)]
write_tab(changes,'preference_change_from_all_genome_full')
p <- ggplot(summary,aes(mode,median_delta,color=region,group=region))+geom_hline(yintercept=0,color='grey65')+
  geom_line()+geom_point(size=2)+facet_wrap(~motif)+scale_color_manual(values=region_colors)+
  labs(title='How median continuous preference changes',x=NULL,y='Median delta-delta log2 enrichment',color=NULL,
    caption='Fixed peak universe; selected-read CPM; 0.5 CPM pseudocount. Lines connect analysis variants, not time points.')+
  theme(axis.text.x=element_text(angle=25,hjust=1))
save_plot(p,'continuous_preference_median_changes',13,5.5)
p <- ggplot(libs,aes(mode,total_selected_reads,color=paste(construct,assay),group=sample_id))+
  geom_line(alpha=.6)+geom_point()+facet_wrap(~motif)+scale_color_manual(values=colors)+scale_y_log10()+
  labs(title='Selected-library sizes across analysis variants',x=NULL,y='Selected reads (log scale)',color=NULL)+
  theme(axis.text.x=element_text(angle=25,hjust=1))
save_plot(p,'selected_library_sizes',13,5.5)
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
cat('Completed:',length(universe),'peaks;',nrow(profiles),'profile rows\n')
