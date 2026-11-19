suppressPackageStartupMessages({library(data.table);library(ggplot2)})
base <- 'output/followup_analysis/read_annotation_comparisons'
out <- file.path(base,'additional_figures')
source <- file.path(base,c('source_library_fractions.tsv','324_vs_328_statistics.tsv'))
hash <- tools::md5sum(source)
a <- fread(source[1])[assay=='CLIP']
original_tests <- fread(source[2])
colors <- c('324'='#4C78A8','328'='#E17C35')
maxfrac <- max(a$fractional_fraction)
for(version in c('intron_3UTR','combined_intergenic_other')) {
  if(version=='intron_3UTR') {
    z <- a[category %in% c('intron','3UTR')]
    tests <- original_tests[category %in% c('intron','3UTR')]
    cats <- c('3UTR','intron');labs <- c("3' UTR",'Intron')
    caption <- 'q values retain the original six-category BH correction; fractions are not renormalized.'
  } else {
    z <- copy(a)
    z[category %in% c('intergenic','other'),category:='intergenic_other']
    z <- z[,.(fractional_fraction=sum(fractional_fraction)),
      by=.(sample_id,construct,biological_replicate,category)]
    stopifnot(all(abs(z[,sum(fractional_fraction),by=sample_id]$V1-1)<1e-8))
    cats <- c('3UTR','5UTR','CDS','intron','intergenic_other')
    labs <- c("3' UTR","5' UTR",'CDS','Intron','Intergenic + Other')
    tests <- rbindlist(lapply(cats,function(cat) {
      x<-z[category==cat & construct==324,fractional_fraction]
      y<-z[category==cat & construct==328,fractional_fraction]
      v<-c(x,y);observed<-mean(x)-mean(y)
      perms<-apply(combn(4,2),2,function(ix) mean(v[ix])-mean(v[-ix]))
      data.table(category=cat,mean_324=mean(x),mean_328=mean(y),
        difference_324_minus_328_pp=100*observed,
        p_two_sided=mean(abs(perms)>=abs(observed)-1e-14))
    }))
    tests[,q_BH:=p.adjust(p_two_sided,'BH')]
    caption <- 'Intergenic + Other is summed within each library. Tests recomputed; BH correction across five categories.'
  }
  z[,isoform:=factor(construct,levels=c(324,328),labels=c('324','328'))]
  z[,panel:=factor(category,levels=cats,labels=labs)]
  z[,plot_x:=as.numeric(isoform)+fifelse(biological_replicate=='A',-.055,.055)]
  sm <- z[,.(mean_fraction=mean(fractional_fraction),n=.N),by=.(category,panel,isoform)]
  stopifnot(all(sm$n==2))
  sm[,plot_x:=as.numeric(isoform)]
  tests[,panel:=factor(category,levels=cats,labels=labs)]
  tests[,label:=sprintf('324 - 328: %+.1f pp; q = %.3f',difference_324_minus_328_pp,q_BH)]
  p <- ggplot(sm,aes(plot_x,mean_fraction,fill=isoform))+
    geom_col(width=.58,alpha=.8)+
    geom_point(data=z,aes(y=fractional_fraction),shape=21,fill='white',color='black',size=2.1)+
    geom_text(data=tests,aes(x=1.5,y=maxfrac+.045,label=label),inherit.aes=FALSE,size=3.2)+
    facet_wrap(~panel,ncol=if(version=='intron_3UTR') 2 else 3)+
    scale_fill_manual(values=colors)+
    scale_x_continuous(breaks=c(1,2),labels=c('324 CLIP','328 CLIP'),limits=c(.5,2.5))+
    scale_y_continuous(labels=scales::label_percent(),limits=c(0,maxfrac+.075),expand=expansion(mult=c(0,.01)))+
    theme_classic(base_size=12)+theme(legend.position='none',
      strip.background=element_rect(fill='#F2F2F2',color='grey70'),strip.text=element_text(face='bold'),
      plot.caption=element_text(hjust=0,size=9),panel.spacing=grid::unit(1,'lines'))+
    labs(title='CLIP read-category fractions: 324 versus 328',
      subtitle='Bars: replicate means; points: individual replicates. No SMI normalization.',
      x=NULL,y='Fractional annotation weight / eligible aligned reads',
      caption=paste(caption,'Exact two-sided library permutation tests; n = 2 per isoform. pp = percentage points.',sep='\n'))
  stem <- file.path(out,paste0('read_categories_324_vs_328_',version))
  for(ext in c('png','pdf')) ggsave(paste0(stem,'.',ext),p,
    width=if(version=='intron_3UTR') 10 else 12,
    height=if(version=='intron_3UTR') 5.5 else 8,dpi=250,bg='white')
  fwrite(z,paste0(stem,'_library_values.tsv'),sep='\t')
  fwrite(sm,paste0(stem,'_means.tsv'),sep='\t')
  fwrite(tests,paste0(stem,'_statistics.tsv'),sep='\t')
}
stopifnot(identical(hash,tools::md5sum(source)))
writeLines(c('# Additional direct CLIP comparisons','',
  'Both versions retain the original denominator (all eligible aligned reads), replicate averaging, colors and y-axis range. No SMI is used.',
  '', 'The intron/3UTR figure is a display subset: no renormalization and no change to the original six-category BH-adjusted q values.',
  '', 'The five-category figure sums intergenic and other fractions within each library before calculating means, effects and exact library-level permutation tests. All five categories are then BH-adjusted together. Other category fractions remain unchanged.',
  '', 'The exact two-sided 2-versus-2 library tests have minimum P = 1/3. The q values are exploratory; a lack of significance is not evidence of equivalence.',
  '', 'Reproduce from the project root: `Rscript output/followup_analysis/read_annotation_comparisons/additional_figures/plot_additional.R`.',
  '', 'Per-library values, means and statistics are saved beside each PNG/PDF. Original figures and source tables are unchanged.'),file.path(out,'README.md'))
fwrite(data.table(path=source,md5=unname(hash)),file.path(out,'input_provenance.tsv'),sep='\t')
