#!/usr/bin/env Rscript
# Run from repository root after run_analysis.R. Primary results are read-only.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(readr);library(ggplot2);library(fgsea)})
base <- 'output/proteomics_hpa'
out <- file.path(base,'sensitivity')
dir.create(file.path(out,'figures'),recursive=TRUE,showWarnings=FALSE)
primary_files <- list.files(base,full.names=TRUE)
primary_files <- primary_files[!dir.exists(primary_files)]
before <- tools::md5sum(primary_files)
x <- readRDS(file.path(base,'analysis_objects.rds'))
cg <- read_tsv(file.path(base,'isoform_classification_genes.tsv'),show_col_types=FALSE)
save_tsv <- function(z,n)write_tsv(as.data.frame(z),file.path(out,paste0(n,'.tsv')),na='NA')
join_text <- function(z)paste(sort(unique(z[!is.na(z)&nzchar(z)])),collapse=';')
# Exclude exact zeros AFTER the primary gene-level median collapse. No tolerance or p filter.
ranked <- x$ranked %>% filter(rank_stat!=0) %>% rename(primary_rank=rank) %>% mutate(rank=row_number())
stats <- x$stats[x$stats!=0]
stopifnot(nrow(ranked)==1062,length(stats)==1062,identical(names(stats),ranked$gene_id))
save_tsv(ranked,'ranked_nonzero_genes')
save_tsv(filter(x$ranked,rank_stat==0),'excluded_exact_zero_genes')
write_tsv(select(ranked,gene_id,rank_stat),file.path(out,'ranked_nonzero.rnk'),col_names=FALSE)
# Fixed primary-confidence-annotated universe for BOTH confidence collections.
bg <- x$ranked$gene_id[x$ranked$mapped_to_hpa & x$ranked$primary_hpa_localization]
stopifnot(length(bg)==1439)
save_tsv(filter(x$ranked,gene_id %in% bg),'ora_primary_annotated_background')
save_tsv(filter(cg,gene_id %in% bg),'ora_restricted_hit_classifications')
save_tsv(filter(cg,!gene_id %in% bg),'ora_excluded_unannotated_classifications')
hit_counts <- cg %>% mutate(in_annotated_universe=gene_id %in% bg) %>% group_by(classification) %>%
  summarise(primary_n=n(),annotated_universe_n=sum(in_annotated_universe),excluded_n=sum(!in_annotated_universe),.groups='drop')
save_tsv(hit_counts,'ora_hit_set_counts')
gg<-list();oo<-list();ll<-list()
for(key in names(x$collections)) {
  paths<-x$collections[[key]]
  sizes<-x$gsea %>% filter(analysis==key) %>% select(analysis,localization,hpa_size) %>%
    mutate(measured_size=vapply(localization,function(n)length(intersect(paths[[n]],names(stats))),integer(1)),
      status=case_when(measured_size<5~'fewer_than_5_measured_genes',measured_size>=.95*length(stats)~'at_least_95_percent_of_universe',TRUE~'tested'))
  eligible<-sizes$localization[sizes$status=='tested']
  set.seed(20260924)
  fg<-as_tibble(fgseaMultilevel(paths[eligible],stats,minSize=5,maxSize=length(stats)-1,eps=0,nPermSimple=10000,
    sampleSize=1001,BPPARAM=BiocParallel::SerialParam(),gseaParam=1)) %>% rename(localization=pathway,measured_size_fgsea=size) %>%
    mutate(padj=p.adjust(pval,'BH'),direction=ifelse(NES>0,'324','328'),analysis=key,
      leading_edge_gene_ids=vapply(leadingEdge,join_text,character(1)),
      leading_edge_symbols=vapply(leadingEdge,function(g)join_text(ranked$gene_symbol[match(g,ranked$gene_id)]),character(1)),
      leading_edge_proteins=vapply(leadingEdge,function(g)join_text(ranked$protein_labels[match(g,ranked$gene_id)]),character(1)))
  stopifnot(!anyNA(fg$pval),!anyNA(fg$NES))
  ll[[key]]<-fg %>% select(analysis,localization,NES,padj,direction,leadingEdge) %>% unnest_longer(leadingEdge,values_to='gene_id') %>% left_join(ranked,by='gene_id')
  gg[[key]]<-left_join(sizes,select(fg,-analysis,-leadingEdge),by='localization')
  save_tsv(gg[[key]],paste0('gsea_nonzero_',key))
  for(iso in c('324','328')) {
    hit<-intersect(cg$gene_id[cg$classification==paste0(iso,'-specific')],bg)
    rr<-bind_rows(lapply(sizes$localization,function(term) {
      loc<-intersect(paths[[term]],bg);genes<-intersect(hit,loc)
      a<-length(genes);b<-length(hit)-a;c<-length(loc)-a;dd<-length(bg)-a-b-c
      tested<-length(loc)>=5 && length(loc)<.95*length(bg) && length(hit)>0
      ft<-if(tested)fisher.test(matrix(c(a,b,c,dd),2,byrow=TRUE),alternative='greater') else NULL
      tibble(analysis=key,isoform=iso,localization=term,hit_in_location=a,total_hits=length(hit),background_in_location=length(loc),total_background=length(bg),
        odds_ratio=if(tested)unname(ft$estimate) else NA_real_,pval=if(tested)ft$p.value else NA_real_,status=if(tested)'tested' else 'insufficient_representation',
        contributing_gene_ids=join_text(genes),contributing_symbols=join_text(x$ranked$gene_symbol[match(genes,x$ranked$gene_id)]),
        contributing_proteins=join_text(x$ranked$protein_labels[match(genes,x$ranked$gene_id)]))
    })) %>% mutate(padj=p.adjust(pval,'BH'))
    oo[[paste(key,iso)]]<-rr
    save_tsv(rr,paste0('ora_annotated_',key,'_',iso))
  }
}
gsea<-bind_rows(gg);ora<-bind_rows(oo);leading<-bind_rows(ll)
save_tsv(gsea,'gsea_nonzero_all');save_tsv(ora,'ora_annotated_all')
save_tsv(leading,'leading_edges_nonzero_all');save_tsv(filter(leading,padj<.05),'leading_edges_nonzero_significant')
gc<-full_join(x$gsea %>% select(analysis,localization,hpa_size,measured_size,ES,NES,pval,padj,status),
  gsea %>% select(analysis,localization,measured_size,ES,NES,pval,padj,status),by=c('analysis','localization'),suffix=c('_primary','_nonzero')) %>%
  mutate(NES_change=NES_nonzero-NES_primary,direction_changed=sign(NES_primary)!=sign(NES_nonzero),
    FDR_conclusion_changed=(padj_primary<.05)!=(padj_nonzero<.05))
oc<-full_join(x$ora %>% select(analysis,isoform,localization,hit_in_location,total_hits,background_in_location,total_background,odds_ratio,pval,padj,status),
  ora %>% select(analysis,isoform,localization,hit_in_location,total_hits,background_in_location,total_background,odds_ratio,pval,padj,status),
  by=c('analysis','isoform','localization'),suffix=c('_primary','_annotated')) %>%
  mutate(FDR_conclusion_changed=(padj_primary<.05)!=(padj_annotated<.05))
save_tsv(gc,'gsea_primary_vs_nonzero');save_tsv(oc,'ora_primary_vs_annotated')
gs<-gsea %>% group_by(analysis) %>% summarise(tested=sum(status=='tested'),significant=sum(padj<.05,na.rm=TRUE),
  minimum_p=min(pval,na.rm=TRUE),minimum_FDR=min(padj,na.rm=TRUE),.groups='drop')
os<-ora %>% group_by(analysis,isoform) %>% summarise(total_background=first(total_background),total_hits=first(total_hits),tested=sum(status=='tested'),
  significant=sum(padj<.05,na.rm=TRUE),minimum_p=min(pval,na.rm=TRUE),minimum_FDR=min(padj,na.rm=TRUE),.groups='drop')
save_tsv(gs,'gsea_nonzero_summary');save_tsv(os,'ora_annotated_summary')
theme_set(theme_minimal(base_size=11)+theme(panel.grid.minor=element_blank(),plot.title=element_text(face='bold')))
save_plot<-function(p,n,w=10,h=7) {
  ggsave(file.path(out,'figures',paste0(n,'.png')),p,width=w,height=h,dpi=180,bg='white')
  ggsave(file.path(out,'figures',paste0(n,'.pdf')),p,width=w,height=h,device=grDevices::pdf,useDingbats=FALSE)
}
for(key in names(x$collections)) {
  z<-filter(gc,analysis==key,status_primary=='tested'|status_nonzero=='tested') %>% arrange(NES_primary) %>%
    mutate(localization=factor(localization,levels=localization))
  p<-ggplot(z,aes(y=localization))+geom_vline(xintercept=0,linetype=2,colour='grey65')+
    geom_segment(aes(x=NES_primary,xend=NES_nonzero,yend=localization),colour='grey65',na.rm=TRUE)+
    geom_point(aes(x=NES_primary,colour='All 1,641 genes'),size=2.5,na.rm=TRUE)+
    geom_point(aes(x=NES_nonzero,colour='1,062 nonzero genes'),size=2.5,na.rm=TRUE)+
    geom_text(aes(x=NES_nonzero,label=ifelse(!is.na(padj_nonzero)&padj_nonzero<.05,'*','')),nudge_y=.3,na.rm=TRUE)+
    scale_colour_manual(values=c('All 1,641 genes'='#707070','1,062 nonzero genes'='#11778a'),name=NULL)+
    labs(title=paste('Exact-zero exclusion:',gsub('_',' ',key)),subtitle='* Nonzero GSEA FDR < 0.05. Categories require at least 5 retained genes.',
      x='NES (negative: 328; positive: 324)',y=NULL)+theme(legend.position='bottom')
  save_plot(p,paste0('GSEA_comparison_',key),10,max(5,.25*nrow(z)+2))
}
for(level in c('fine','collapsed')) {
  z<-filter(oc,analysis==paste0(level,'_primary'),status_primary=='tested') %>%
    select(localization,isoform,padj_primary,padj_annotated) %>% pivot_longer(starts_with('padj'),names_to='universe',values_to='FDR') %>%
    mutate(column=paste(isoform,ifelse(universe=='padj_primary','all mapped','annotated only')),score=-log10(FDR))
  p<-ggplot(z,aes(column,localization,fill=score))+geom_tile(colour='white')+
    geom_text(aes(label=ifelse(FDR<.05,'*','')),size=5)+scale_fill_gradient(low='white',high='#11778a',name='-log10(FDR)')+
    labs(title=paste('ORA background sensitivity:',level),subtitle='* FDR < 0.05. Both hit sets are restricted to each tested universe.',x=NULL,y=NULL)+
    theme(axis.text.x=element_text(angle=20,hjust=1))
  save_plot(p,paste0('ORA_comparison_',level),9,max(5,.25*n_distinct(z$localization)+2))
}
# Independent numerical checks, without rerunning or mutating the primary analysis.
es_errors<-vapply(which(gsea$status=='tested'),function(i) {
  hit<-names(stats) %in% x$collections[[gsea$analysis[i]]][[gsea$localization[i]]]
  running<-c(0,cumsum(ifelse(hit,abs(stats)/sum(abs(stats[hit])),-1/sum(!hit))))
  es<-if(max(running)> -min(running))max(running) else min(running)
  abs(es-gsea$ES[i])
},numeric(1))
tested_ora<-filter(ora,status=='tested')
hyper<-with(tested_ora,phyper(hit_in_location-1,background_in_location,total_background-background_in_location,total_hits,lower.tail=FALSE))
checks<-tibble(check=c('Exactly 579 zero genes excluded after primary collapse','Nonzero scores unchanged and unique','Fixed 1439-gene ORA background',
 'All ORA hit sets restricted to background','Independent GSEA running sums','Fisher matches hypergeometric tail','Primary output checksums unchanged'),
 passed=c(sum(x$stats==0)==579 && length(stats)==1062,identical(stats,x$stats[names(stats)]) && !anyDuplicated(names(stats)),
 all(ora$total_background==1439),all(ora$hit_in_location<=ora$total_hits)&all(ora$total_hits<=1439),max(es_errors)<1e-10,
 max(abs(hyper-tested_ora$pval))<1e-10,identical(before,tools::md5sum(primary_files))))
stopifnot(all(checks$passed))
save_tsv(checks,'validation_checks')
save_tsv(tibble(file=names(before),md5=unname(before)),'primary_input_checksums')
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
print(gs);print(os);print(hit_counts)
print(filter(gsea,padj<.05) %>% select(analysis,localization,NES,pval,padj))
print(filter(ora,padj<.05) %>% select(analysis,isoform,localization,odds_ratio,pval,padj))
