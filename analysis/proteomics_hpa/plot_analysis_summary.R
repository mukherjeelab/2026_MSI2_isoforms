#!/usr/bin/env Rscript
# Presentation figures from completed results only; no statistical analyses are rerun.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(readr);library(ggplot2)})
base<-'output/proteomics_hpa'
out<-file.path(base,'summary_figures')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
read<-function(n)read_tsv(file.path(base,n),show_col_types=FALSE)
save<-function(z,n)write_tsv(as.data.frame(z),file.path(out,paste0(n,'.tsv')),na='NA')
theme_set(theme_minimal(base_size=12)+theme(panel.grid.minor=element_blank(),panel.grid.major.y=element_blank(),
 plot.title=element_text(face='bold',size=18,colour='#202c38'),plot.subtitle=element_text(size=11,colour='#4b5866',margin=margin(b=14)),
 plot.caption=element_text(size=9,colour='#4b5866',hjust=0,margin=margin(t=12)),plot.margin=margin(16,22,14,14),
 axis.title=element_text(colour='#384657'),axis.text=element_text(colour='#384657'),strip.text=element_text(face='bold',size=11),
 legend.position='bottom',legend.title=element_blank()))
plots<-list()
export<-function(p,name,w,h) {
 plots[[name]]<<-p
 ggsave(file.path(out,paste0(name,'.png')),p,width=w,height=h,dpi=200,bg='white')
 ggsave(file.path(out,paste0(name,'.pdf')),p,width=w,height=h,device=grDevices::pdf,useDingbats=FALSE)
}
family<-function(d,group,label,question,pcol) {
 d %>% group_by(across(all_of(group))) %>% summarise(category_tests=sum(!is.na(.data[[pcol]])),
  min_FDR=min(.data[[pcol]],na.rm=TRUE),significant_005=sum(.data[[pcol]]<.05,na.rm=TRUE),.groups='drop') %>%
  unite('family',all_of(group),sep=' / ') %>% mutate(method=label,question=question)
}
g1<-read('gsea_all.tsv');g2<-read('sensitivity/gsea_nonzero_all.tsv')
o1<-read('ora_all.tsv');o2<-read('sensitivity/ora_annotated_all.tsv')
o3<-read('reactome_matched_ora/hpa_ora_all_exact_reactome_sets.tsv')
o4<-read('broad_exploratory_ora/broad_hpa_ora_all.tsv')
direct<-read('direct_localization_comparison/direct_comparison_all.tsv')
labels<-c('Ranked GSEA: all genes','Ranked GSEA: exclude exact zeros','ORA: original HPA-mapped sets','ORA: annotated-only background',
 'ORA: exact Reactome input sets','ORA: eight broad categories','Direct 324 versus 328 proportions')
fam<-bind_rows(family(g1,'analysis',labels[1],'Ranked association','padj'),family(g2,'analysis',labels[2],'Ranked association','padj'),
 family(o1,c('analysis','isoform'),labels[3],'Enrichment versus measured background','padj'),
 family(o2,c('analysis','isoform'),labels[4],'Enrichment versus measured background','padj'),
 family(o3,c('pair','resolution','isoform'),labels[5],'Enrichment versus measured background','BH_FDR'),
 family(o4,c('pair','isoform'),labels[6],'Enrichment versus measured background','BH_FDR'),
 family(direct,'pair',labels[7],'Between-set difference','BH_FDR'))
save(fam,'all_analysis_families')
summary<-fam %>% group_by(method,question) %>% summarise(families=n(),category_tests=sum(category_tests),minimum_FDR=min(min_FDR),
 maximum_family_minimum_FDR=max(min_FDR),significant_005=sum(significant_005),.groups='drop') %>%
 mutate(method=factor(method,levels=rev(labels)))
save(summary,'method_summary')
stopifnot(nrow(fam)==45,sum(fam$significant_005)==0,min(fam$min_FDR)>.1)
fam<-fam %>% mutate(method=factor(method,levels=rev(labels)))
colours<-c('Ranked association'='#487a98','Enrichment versus measured background'='#687c69','Between-set difference'='#92596e')
p<-ggplot(summary,aes(y=method))+
 annotate('rect',xmin=0,xmax=.1,ymin=-Inf,ymax=Inf,fill='#f1e3db',alpha=.6)+
 geom_vline(xintercept=c(.05,.1),linetype=c('dotted','dashed'),colour='#9e6551',linewidth=.5)+
 geom_segment(aes(x=minimum_FDR,xend=maximum_family_minimum_FDR,yend=method),colour='#d0d6da',linewidth=2)+
 geom_point(data=fam,aes(x=min_FDR,colour=question),size=2,alpha=.5)+
 geom_point(aes(x=minimum_FDR,colour=question),size=4)+
 geom_text(aes(x=minimum_FDR,label=sprintf('%.3f',minimum_FDR)),nudge_y=.27,size=3.5,fontface='bold',colour='#263441')+
 geom_text(aes(x=1.10,label=paste0(families,' families')),hjust=1,size=3.3,colour='#4b5866')+
 scale_colour_manual(values=colours)+scale_x_continuous(limits=c(0,1.13),breaks=c(0,.1,.25,.5,.75,1),labels=c('0','0.10','0.25','0.50','0.75','1.00'))+
 labs(title='No significant HPA localization result across the tested approaches',
 subtitle='Each small dot is the lowest FDR in one correction family. Large dots mark the minimum for each approach.',
 x='Lowest category FDR within each analysis family',y=NULL,
 caption='Shaded region: FDR < 0.10. Dotted line: 0.05; dashed line: 0.10. No result crosses either threshold.\nThese are related analyses of the same experiment, not independent replications. FDR is not an effect size.')+
 guides(colour=guide_legend(nrow=2,byrow=TRUE))
export(p,'01_analysis_overview',13,7.5)

terms<-c('Nuclear/nucleoplasmic','Nucleolus','Cytosolic','Mitochondrial','ER/Golgi','Vesicular/endolysosomal','Plasma membrane/junctional','Cytoskeletal')
pair_names<-c(cluster='Full clusters\n324: n=421 | 328: n=609',hits='EGFP-filtered clusters\n324: n=237 | 328: n=240',specific='Isoform-specific subsets\n324: n=172 | 328: n=300')
dd<-direct %>% mutate(localization=factor(localization,levels=rev(terms)),pair_label=factor(pair_names[pair],levels=pair_names))
save(dd %>% select(pair,localization,starts_with('full_'),starts_with('test_'),odds_ratio_324_vs_328,OR_CI95_lower,OR_CI95_upper,p_value,BH_FDR),'direct_comparison_plot_data')
p<-ggplot(dd,aes(y=localization))+
 geom_segment(aes(x=100*full_324_proportion,xend=100*full_328_proportion,yend=localization),colour='#c5cdd3',linewidth=.7)+
 geom_point(aes(x=100*full_324_proportion,colour='MSI2 324'),size=3)+
 geom_point(aes(x=100*full_328_proportion,colour='MSI2 328'),size=3)+
 facet_wrap(~pair_label,nrow=1)+scale_colour_manual(values=c('MSI2 324'='#b44a50','MSI2 328'='#287ba5'))+
 scale_x_continuous(limits=c(0,50),breaks=seq(0,50,10),labels=function(x)paste0(x,'%'))+
 labs(title='Direct comparison: no broad category differs significantly between hit sets',
 subtitle='Original-set membership proportions across all three matched pairs; eight two-sided tests per pair.',
 x='Share of original hit set annotated to category',y=NULL,
 caption='Categories overlap: a gene can contribute to several rows, so percentages do not sum to 100%.\nSpecific sets share EIF4A2. Full proportions are shown; inference uses nonshared genes (171 versus 299). Lowest direct-test FDR: 0.498.')
export(p,'02_direct_category_proportions',15,7)

forest<-dd %>% mutate(pair_label=as.character(pair_label),pair_label=ifelse(pair=='specific','Isoform-specific: nonshared genes\n324: n=171 | 328: n=299',pair_label),
 pair_label=factor(pair_label,levels=c(pair_names[1:2],'Isoform-specific: nonshared genes\n324: n=171 | 328: n=299')))
p<-ggplot(forest,aes(y=localization,x=odds_ratio_324_vs_328))+
 geom_vline(xintercept=1,linetype='dashed',colour='#6e7a85')+
 geom_segment(aes(x=OR_CI95_lower,xend=OR_CI95_upper,yend=localization),colour='#71899c',linewidth=.7)+
 geom_point(size=2.6,colour='#315873')+facet_wrap(~pair_label,nrow=1)+
 scale_x_log10(breaks=c(.25,.5,1,2,4),labels=c('0.25','0.5','1','2','4'))+
 labs(title='Direct comparisons remain compatible with equal localization odds',
 subtitle='Points: odds ratio for 324 relative to 328. Lines: exact, unadjusted 95% confidence intervals.',
 x='Odds ratio on a log scale (below 1: higher in 328; above 1: higher in 324)',y=NULL,
 caption='Every confidence interval includes OR = 1. BH correction is across eight categories per pair.\nSpecific-subset inference uses nonshared genes. Absence of significance does not establish equivalence; intervals show the remaining uncertainty.')
export(p,'03_direct_odds_ratios',15,7)

gc<-read('sensitivity/gsea_primary_vs_nonzero.tsv') %>% filter(grepl('^collapsed_',analysis),status_primary=='tested') %>%
 mutate(confidence=ifelse(analysis=='collapsed_primary','Enhanced + Supported + Approved','Enhanced + Supported'))
gl<-bind_rows(gc %>% transmute(localization,confidence,ranking='All 1,641 genes',NES=NES_primary,FDR=padj_primary),
 gc %>% transmute(localization,confidence,ranking='1,062 nonzero genes',NES=NES_nonzero,FDR=padj_nonzero))
order<-c('Nucleoplasm','Nucleoli','Nuclear membrane','Cytosol','Mitochondria','Endoplasmic reticulum','Golgi apparatus','Vesicular system',
 'Plasma membrane/cell junctions','Actin-filament associated','Microtubule associated','Intermediate filaments','Centrosome','Primary cilium','Sperm structures')
gl<-gl %>% mutate(localization=factor(localization,levels=rev(order)),confidence=factor(confidence,levels=c('Enhanced + Supported + Approved','Enhanced + Supported')))
save(gl,'ranked_sensitivity_plot_data')
p<-ggplot(gl,aes(NES,localization,colour=ranking))+geom_vline(xintercept=0,linetype='dashed',colour='#85919b')+
 geom_point(position=position_dodge(width=.45),size=2.7,na.rm=TRUE)+facet_wrap(~confidence,nrow=1)+
 scale_colour_manual(values=c('All 1,641 genes'='#758595','1,062 nonzero genes'='#178a8d'))+
 labs(title='Ranked enrichment remains nonsignificant under both sensitivity choices',
 subtitle='Collapsed HPA compartments; all genes versus removal of exact-zero scores, at two annotation-confidence levels.',
 x='Normalized enrichment score (negative: 328; positive: 324)',y=NULL,
 caption='No category shown has FDR < 0.10. NES direction alone is not significant enrichment.\nSperm structures has fewer than five retained high-confidence nonzero genes, so that sensitivity estimate is not tested.')
export(p,'04_ranked_sensitivity_comparison',13,8)

focus<-dd %>% filter(localization %in% c('Nucleolus','Mitochondrial')) %>%
 select(pair,pair_label,localization,full_324_in,full_324_total,full_328_in,full_328_total,p_value,BH_FDR) %>%
 pivot_longer(c(full_324_in,full_328_in),names_to='isoform',values_to='overlap') %>%
 mutate(isoform=ifelse(isoform=='full_324_in','MSI2 324','MSI2 328'),n=ifelse(isoform=='MSI2 324',full_324_total,full_328_total),proportion=100*overlap/n,
 label=paste0(overlap,'/',n,'\n',sprintf('%.1f%%',proportion)),localization=factor(localization,levels=c('Nucleolus','Mitochondrial')))
save(focus,'hypothesis_plot_data')
p<-ggplot(focus,aes(isoform,proportion,fill=isoform))+geom_col(width=.55)+
 geom_text(aes(label=label),vjust=-.3,size=3.2,lineheight=1)+facet_grid(localization~pair_label)+
 scale_fill_manual(values=c('MSI2 324'='#b44a50','MSI2 328'='#287ba5'))+scale_y_continuous(limits=c(0,21),breaks=c(0,5,10,15,20),labels=function(x)paste0(x,'%'))+
 labs(title='Nucleolar and mitochondrial membership: no detected isoform difference',
 subtitle='Original hit-set proportions. Neither hypothesis passes the direct comparison in any matched pair.',
 x=NULL,y='Share of original hit set',
 caption='EGFP-filtered mitochondrial example: 36/237 versus 37/240; odds ratio 0.983 (95% CI 0.578-1.670), Fisher p = 1.00, FDR = 1.00.\nBars are descriptive proportions, not replicate means. Exact odds-ratio intervals for the comparisons appear in Figure 3.')+
 theme(legend.position='none',panel.grid.major.x=element_blank())
export(p,'05_nucleolus_mitochondria_focus',13,8)

pdf(file.path(out,'HPA_localization_summary_all_figures.pdf'),width=15,height=8.5,useDingbats=FALSE)
for(p in plots)print(p)
dev.off()
stopifnot(all(dd$OR_CI95_lower<=1&dd$OR_CI95_upper>=1),all(direct$BH_FDR>.1))
save(tibble(check=c('45 correction families represented','All family minimum FDRs exceed 0.10','All 24 direct confidence intervals include 1','All 24 direct FDRs exceed 0.10'),passed=TRUE),'figure_data_checks')
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
print(summary %>% arrange(method) %>% select(method,families,minimum_FDR,significant_005))
