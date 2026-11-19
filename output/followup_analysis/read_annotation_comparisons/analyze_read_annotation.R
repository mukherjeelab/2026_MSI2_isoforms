suppressPackageStartupMessages({library(data.table);library(ggplot2)})
out <- 'output/followup_analysis/read_annotation_comparisons'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
input <- 'output/followup_analysis/tables/read_annotation_summary.tsv'
hash <- tools::md5sum(input)
a <- fread(input)
stopifnot(nrow(a)==48,uniqueN(a$sample_id)==8,!anyDuplicated(a[,.(sample_id,category)]),
  all(abs(a[,sum(fractional_fraction),by=sample_id]$V1-1)<1e-8))
cats <- c('3UTR','5UTR','CDS','intron','intergenic','other')
labels <- c('3UTR'="3' UTR",'5UTR'="5' UTR",CDS='CDS',intron='Intron',intergenic='Intergenic',other='Other')
stopifnot(setequal(a$category,cats))
a[,isoform:=factor(construct,levels=c(324,328),labels=c('324','328'))]
fwrite(a,file.path(out,'source_library_fractions.tsv'),sep='\t')
smi <- a[assay=='SMI',.(SMI_mean=mean(fractional_fraction),SMI_sd=sd(fractional_fraction),
  SMI_min=min(fractional_fraction),SMI_max=max(fractional_fraction),n_SMI=.N),by=category]
stopifnot(all(smi$n_SMI==4),all(smi$SMI_mean>0),abs(sum(smi$SMI_mean)-1)<1e-8)
clip <- merge(a[assay=='CLIP'],smi,by='category')
clip[,fold_over_pooled_SMI:=fractional_fraction/SMI_mean]
summary <- clip[,.(mean_fraction=mean(fractional_fraction),sd_fraction=sd(fractional_fraction),
  min_fraction=min(fractional_fraction),max_fraction=max(fractional_fraction),n_CLIP=.N,
  SMI_mean=SMI_mean[1],SMI_sd=SMI_sd[1],n_SMI=n_SMI[1],
  fold_over_pooled_SMI=mean(fold_over_pooled_SMI),
  difference_from_SMI_pp=100*(mean(fractional_fraction)-SMI_mean[1])),by=.(category,isoform)]
stopifnot(all(summary$n_CLIP==2))

# Enumerate all label assignments of whole-library fractions; reads are not replicates.
permutation_test <- function(x,y) {
  n<-length(x); combined<-c(x,y)
  combinations<-combn(length(combined),n)
  diffs<-apply(combinations,2,function(ix) mean(combined[ix])-mean(combined[-ix]))
  observed<-mean(x)-mean(y)
  p<-mean(abs(diffs)>=abs(observed)-1e-14)
  list(observed=observed,p=p,null=diffs,assignments=ncol(combinations))
}
tests_smi <- list(); tests_iso <- list(); permutation_rows <- list()
for(cat in cats) {
  for(iso in c('324','328')) {
    x <- clip[category==cat & isoform==iso,fractional_fraction]
    y <- a[category==cat & assay=='SMI',fractional_fraction]
    t <- permutation_test(x,y)
    tests_smi[[paste(cat,iso)]] <- data.table(category=cat,isoform=iso,
      comparison=paste0(iso,' CLIP vs pooled SMI'),mean_CLIP=mean(x),mean_SMI=mean(y),
      difference_pp=100*t$observed,fold=mean(x)/mean(y),n_CLIP=length(x),n_SMI=length(y),
      p_two_sided=t$p,label_assignments=t$assignments)
    permutation_rows[[paste(cat,iso)]] <- data.table(family='CLIP_vs_SMI',category=cat,
      comparison=paste0(iso,' CLIP vs pooled SMI'),assignment=seq_along(t$null),difference=t$null)
  }
  x<-clip[category==cat & isoform=='324',fractional_fraction]
  y<-clip[category==cat & isoform=='328',fractional_fraction]
  t<-permutation_test(x,y)
  tests_iso[[cat]] <- data.table(category=cat,comparison='324 CLIP vs 328 CLIP',
    mean_324=mean(x),mean_328=mean(y),difference_324_minus_328_pp=100*t$observed,
    ratio_324_over_328=mean(x)/mean(y),n_324=length(x),n_328=length(y),
    p_two_sided=t$p,label_assignments=t$assignments)
  permutation_rows[[paste(cat,'isoforms')]] <- data.table(family='324_vs_328',category=cat,
    comparison='324 CLIP vs 328 CLIP',assignment=seq_along(t$null),difference=t$null)
}
ts <- rbindlist(tests_smi);ti <- rbindlist(tests_iso)
ts[,q_BH:=p.adjust(p_two_sided,'BH')];ti[,q_BH:=p.adjust(p_two_sided,'BH')]
fwrite(smi,file.path(out,'pooled_SMI_reference.tsv'),sep='\t')
fwrite(summary,file.path(out,'replicate_mean_fractions.tsv'),sep='\t')
fwrite(clip,file.path(out,'CLIP_plot_inputs.tsv'),sep='\t')
fwrite(ts,file.path(out,'CLIP_vs_pooled_SMI_statistics.tsv'),sep='\t')
fwrite(ti,file.path(out,'324_vs_328_statistics.tsv'),sep='\t')
fwrite(rbindlist(permutation_rows),file.path(out,'exact_permutation_distributions.tsv'),sep='\t')

colors <- c('324'='#4C78A8','328'='#E17C35')
theme_p <- theme_classic(base_size=12)+theme(legend.position='none',
  strip.background=element_rect(fill='#F2F2F2',color='grey70'),strip.text=element_text(face='bold'),
  plot.caption=element_text(hjust=0,size=10),panel.spacing=grid::unit(1,'lines'))
save_plot <- function(p,name) {
  ggsave(file.path(out,paste0(name,'.png')),p,width=12,height=8,dpi=250,bg='white')
  ggsave(file.path(out,paste0(name,'.pdf')),p,width=12,height=8,bg='white',useDingbats=FALSE)
}
summary[,panel:=factor(category,levels=cats,labels=labels[cats])]
clip[,panel:=factor(category,levels=cats,labels=labels[cats])]
clip[,plot_x:=as.numeric(isoform)+fifelse(biological_replicate=='A',-.055,.055)]
summary[,plot_x:=as.numeric(isoform)]
ts[,panel:=factor(category,levels=cats,labels=labels[cats])]
ti[,panel:=factor(category,levels=cats,labels=labels[cats])]
smi[,panel:=factor(category,levels=cats,labels=labels[cats])]
maxfold <- max(clip$fold_over_pooled_SMI)
ts[,`:=`(plot_x=match(isoform,c('324','328')),label=sprintf('q = %.3f',q_BH))]
tops <- clip[,.(top=max(fold_over_pooled_SMI)),by=.(category,isoform)]
ts<-merge(ts,tops,by=c('category','isoform'))
smi[,ref_label:=sprintf('SMI mean = %.1f%% (n = 4)',100*SMI_mean)]
p1 <- ggplot(summary,aes(plot_x,fold_over_pooled_SMI,fill=isoform))+
  geom_col(width=.58,alpha=.8)+geom_hline(yintercept=1,linetype=2,color='grey35')+
  geom_point(data=clip,aes(x=plot_x,y=fold_over_pooled_SMI),shape=21,fill='white',color='black',size=2.1)+
  geom_text(data=ts,aes(x=plot_x,y=pmax(top+.16,1.23),label=label),inherit.aes=FALSE,size=3.3)+
  geom_text(data=smi,aes(x=1.5,y=maxfold+.55,label=ref_label),inherit.aes=FALSE,size=3.3,color='grey30')+
  facet_wrap(~panel,ncol=3)+scale_fill_manual(values=colors)+
  scale_x_continuous(breaks=c(1,2),labels=c('324 CLIP','328 CLIP'),limits=c(.5,2.5))+
  scale_y_continuous(limits=c(0,maxfold+.8),expand=expansion(mult=c(0,.01)))+
  theme_p+labs(title='Read-category enrichment relative to the pooled SMI mean',
    subtitle='Bars: mean of two CLIP libraries / equally weighted mean of all four SMI libraries',
    x=NULL,y='Fold over pooled SMI',
    caption=paste('Dashed line: SMI reference = 1. White points: individual CLIP replicates divided by the same SMI mean.',
      'q: BH-adjusted exact two-sided library-label permutation P (12 tests; 2 CLIP vs 4 SMI).',
      'Tests use individual library fractions, including SMI variability; SMI uncertainty is not shown by the points.',sep='\n'))
save_plot(p1,'read_categories_CLIP_vs_pooled_SMI')

maxfrac <- max(clip$fractional_fraction)
ti[,label:=sprintf('324 - 328: %+.1f pp; q = %.3f',difference_324_minus_328_pp,q_BH)]
p2 <- ggplot(summary,aes(plot_x,mean_fraction,fill=isoform))+
  geom_col(width=.58,alpha=.8)+
  geom_point(data=clip,aes(x=plot_x,y=fractional_fraction),shape=21,fill='white',color='black',size=2.1)+
  geom_text(data=ti,aes(x=1.5,y=maxfrac+.045,label=label),inherit.aes=FALSE,size=3.2)+
  facet_wrap(~panel,ncol=3)+scale_fill_manual(values=colors)+
  scale_x_continuous(breaks=c(1,2),labels=c('324 CLIP','328 CLIP'),limits=c(.5,2.5))+
  scale_y_continuous(labels=scales::label_percent(),limits=c(0,maxfrac+.075),expand=expansion(mult=c(0,.01)))+
  theme_p+labs(title='CLIP read-category fractions: 324 versus 328',
    subtitle='Bars: equally weighted replicate means; no SMI normalization or SMI bars',
    x=NULL,y='Fractional annotation weight / eligible aligned reads',
    caption=paste('White points: individual biological replicates (n = 2 per isoform). pp = percentage-point difference.',
      'q: BH-adjusted exact two-sided library-label permutation P (6 category tests).',
      'Only six label assignments are possible for 2 vs 2 libraries; the minimum two-sided P is 0.333.',sep='\n'))
save_plot(p2,'read_categories_324_vs_328')
stopifnot(identical(hash,tools::md5sum(input)))
fwrite(data.table(path=input,md5=unname(hash)),file.path(out,'input_provenance.tsv'),sep='\t')
writeLines(capture.output(sessionInfo()),file.path(out,'session_info.txt'))
print(ts[,.(category,isoform,fold,difference_pp,p_two_sided,q_BH)])
print(ti[,.(category,mean_324,mean_328,difference_324_minus_328_pp,p_two_sided,q_BH)])
