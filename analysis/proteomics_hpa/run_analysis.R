#!/usr/bin/env Rscript
# Run from the repository root. Sources are frozen in data/hpa_25.1.
suppressPackageStartupMessages({library(readxl); library(dplyr); library(tidyr); library(readr); library(ggplot2); library(fgsea); library(xml2)})
out <- 'output/proteomics_hpa'
ref <- 'data/hpa_25.1'
dir.create(file.path(out,'figures'), recursive=TRUE, showWarnings=FALSE)
save_tsv <- function(x, name) write_tsv(as.data.frame(x), file.path(out,paste0(name,'.tsv')), na='NA')
read_zip <- function(name) {
  z <- file.path(ref,paste0(name,'.tsv.zip'))
  read_tsv(unz(z,unzip(z,list=TRUE)$Name[1]), show_col_types=FALSE, progress=FALSE)
}
join_text <- function(x) paste(sort(unique(x[!is.na(x) & nzchar(x)])),collapse=';')
tokens <- function(x) { if(is.na(x) || !nzchar(x)) return(character()); unique(trimws(unlist(strsplit(x,'[;,]')))) }
hpa <- read_zip('subcellular_location')
lookup <- read_zip('subcellular_locations')
atlas <- read_zip('proteinatlas')
stopifnot(nrow(lookup)==49, !anyDuplicated(atlas$Ensembl), !anyDuplicated(hpa$Gene))
raw <- read_xlsx('data/FredHutchProteomics.xlsx')
d <- raw %>% transmute(source_row=row_number()+1L, protein_label=label, source_symbol=`Gene Symbol`, source_ensembl=`Ensembl Gene ID`,
  fc_328_324=`f_MSI2_328_WT /MSI2_324_WT`, comparison_p=`p_MSI2_328_WT /MSI2_324_WT`,
  rank_stat=-log2(fc_328_324), testable=is.finite(rank_stat)&is.finite(comparison_p)&comparison_p>=0&comparison_p<=1)
# Unique exact Ensembl; resolve multi-ID candidates only with concordant symbol.
# Exact HPA symbol and UniProt accessions rescue absent/outdated/nonhuman IDs.
up <- atlas %>% select(Ensembl,Uniprot) %>% separate_rows(Uniprot,sep='[,;]') %>% mutate(Uniprot=trimws(Uniprot))
sy <- atlas %>% select(Ensembl,Gene)
mapping <- bind_rows(lapply(seq_len(nrow(d)),function(i) {
  r<-d[i,]; ids<-tokens(r$source_ensembl); ids<-sub('\\.[0-9]+$','',ids)
  human<-ids[grepl('^ENSG[0-9]{11}$',ids)]
  ec<-intersect(human,atlas$Ensembl)
  sc<-unique(sy$Ensembl[sy$Gene %in% tokens(r$source_symbol)])
  acc<-sub('^.*_','',r$protein_label)
  uc<-unique(up$Ensembl[!is.na(up$Uniprot)&up$Uniprot==acc])
  gene<-NA_character_; method<-'unmapped'
  if(length(ec)==1 && (length(sc)==0 || ec %in% sc)) {gene<-ec;method<-'exact_ensembl'}
  else if(length(ec)>1 && length(intersect(ec,sc))==1) {gene<-intersect(ec,sc);method<-'multi_id_resolved_by_symbol'}
  else if(length(ec)==0 && length(sc)==1) {gene<-sc;method<-'exact_symbol'}
  else if(length(ec)==0 && length(sc)==0 && length(uc)==1) {gene<-uc;method<-'exact_uniprot'}
  else if(length(ec)+length(sc)+length(uc)>0) method<-'ambiguous_or_conflicting'
  # Retain a uniquely identified human gene absent from HPA in the ranked universe.
  if(method=='unmapped' && length(human)==1) {gene<-human;method<-'source_human_ensembl_not_in_hpa'}
  tibble(gene_id=gene,mapping_method=method,input_id_count=length(ids),input_human_id_count=length(human),
    hpa_ensembl_candidates=join_text(ec),hpa_symbol_candidates=join_text(sc),hpa_uniprot_candidates=join_text(uc),
    mapped_to_hpa=!is.na(gene)&&gene %in% atlas$Ensembl)
}))
d<-bind_cols(d,mapping)
d$gene_symbol<-atlas$Gene[match(d$gene_id,atlas$Ensembl)]
d$gene_symbol<-coalesce(d$gene_symbol,d$source_symbol)
long <- hpa %>% select(gene_id=Gene,gene_symbol=`Gene name`,Enhanced,Supported,Approved,Uncertain) %>%
  pivot_longer(Enhanced:Uncertain,names_to='confidence',values_to='localization') %>%
  filter(!is.na(localization),localization!='') %>% separate_rows(localization,sep=';') %>% mutate(localization=trimws(localization)) %>% distinct()
stopifnot(all(long$localization %in% lookup$`Subcellular location`))
save_tsv(long,'hpa_gene_location_long')
d$any_hpa_localization<-d$gene_id %in% long$gene_id
d$primary_hpa_localization<-d$gene_id %in% long$gene_id[long$confidence!='Uncertain']
d$high_confidence_localization<-d$gene_id %in% long$gene_id[long$confidence %in% c('Enhanced','Supported')]
d$rank_included<-d$testable & !is.na(d$gene_id)
d$rank_exclusion_reason<-ifelse(!d$testable,'invalid_or_missing_comparison',ifelse(is.na(d$gene_id),d$mapping_method,'included'))
dup_counts<-d %>% filter(rank_included) %>% count(gene_id,name='rows_per_gene')
d<-left_join(d,dup_counts,by='gene_id')
save_tsv(d,'protein_identifier_mapping')
save_tsv(filter(d,!rank_included),'unmapped_or_excluded_proteins')
ranked<-d %>% filter(rank_included) %>% group_by(gene_id) %>% summarise(gene_symbol=first(gene_symbol),
  rank_stat=median(rank_stat),n_protein_rows=n(),protein_labels=join_text(protein_label),source_rows=join_text(as.character(source_row)),
  mapped_to_hpa=all(mapped_to_hpa),any_hpa_localization=any(any_hpa_localization),primary_hpa_localization=any(primary_hpa_localization),
  .groups='drop') %>% arrange(desc(rank_stat),gene_id) %>% mutate(rank=row_number())
# Avoid summarise name shadowing for diagnostic range.
ranges<-d %>% filter(rank_included) %>% group_by(gene_id) %>% summarise(min_row_log2FC=min(rank_stat),max_row_log2FC=max(rank_stat),.groups='drop')
ranked<-ranked %>% left_join(ranges,by='gene_id')
save_tsv(ranked,'ranked_324_vs_328_genes')
write_tsv(select(ranked,gene_id,rank_stat),file.path(out,'ranked_324_vs_328.rnk'),col_names=FALSE)
stats<-setNames(ranked$rank_stat,ranked$gene_id)
N<-sum(d$testable)
qc<-tibble(metric=c('input_protein_rows','tested_protein_rows','mapped_to_hpa_rows','rows_with_any_hpa_localization','rows_without_any_hpa_localization',
  'rows_with_primary_localization','rows_with_high_confidence_localization','input_multiple_id_rows','input_multiple_human_id_rows',
  'unresolved_ambiguous_rows','unmapped_rows','duplicate_genes_after_mapping','redundant_rows_collapsed','unique_ranked_genes','zero_ranked_genes','unique_hpa_mapped_background_genes'),
  count=c(nrow(d),N,sum(d$testable&d$mapped_to_hpa),sum(d$testable&d$any_hpa_localization),sum(d$testable&!d$any_hpa_localization),
  sum(d$testable&d$primary_hpa_localization),sum(d$testable&d$high_confidence_localization),sum(d$testable&d$input_id_count>1),sum(d$testable&d$input_human_id_count>1),
  sum(d$testable&d$mapping_method=='ambiguous_or_conflicting'),sum(d$testable&is.na(d$gene_id)),sum(dup_counts$rows_per_gene>1),sum(pmax(dup_counts$rows_per_gene-1,0)),
  nrow(ranked),sum(ranked$rank_stat==0),sum(ranked$mapped_to_hpa))) %>% mutate(percent_of_tested_rows=100*count/N)
save_tsv(qc,'mapping_QC')
# Parse the actual HPA hierarchy from its archived HTML; normalize term casing through lookup.
doc<-read_html(file.path(ref,'subcellular_overview.html'))
nodes<-xml_find_all(doc,"//span[@class='cell_organelle']")
hm<-tibble(term=trimws(sub('^- ','',xml_text(nodes))),path=xml_attr(xml_find_first(nodes,'.//a'),'href')) %>% distinct()
hm$localization<-lookup$`Subcellular location`[match(tolower(hm$term),tolower(lookup$`Subcellular location`))]
hm<-filter(hm,!is.na(localization))
hm$hpa_parent<-gsub('+',' ',sub('^.*/','',hm$path),fixed=TRUE)
parent_names<-c('nuclear membrane'='Nuclear membrane','nucleoli'='Nucleoli','nucleoplasm'='Nucleoplasm','actin filaments'='Actin-filament associated',
 'centrosome'='Centrosome','cytosol'='Cytosol','intermediate filaments'='Intermediate filaments','microtubules'='Microtubule associated',
 'mitochondria'='Mitochondria','endoplasmic reticulum'='Endoplasmic reticulum','golgi apparatus'='Golgi apparatus',
 'plasma membrane'='Plasma membrane/cell junctions','vesicles'='Vesicular system','primary cilium'='Primary cilium','sperm'='Sperm structures')
hm$collapsed_localization<-unname(parent_names[hm$hpa_parent])
hm$source_url<-paste0('https://www.proteinatlas.org',hm$path)
hm<-hm %>% select(localization,collapsed_localization,hpa_parent,source_url) %>% left_join(lookup,by=c('localization'='Subcellular location'))
stopifnot(nrow(hm)==49,!anyNA(hm$collapsed_localization),!anyDuplicated(hm$localization))
save_tsv(hm,'collapsed_location_mapping')
# Reproduce notebook clustering before HPA mapping, using original symbols and all workbook rows.
avg<-raw %>% filter(!is.na(`Gene Symbol`),`Gene Symbol`!='') %>% group_by(gene=`Gene Symbol`) %>% summarise(
  EGFP=mean(c(EGFP_Rep_B_n_n,EGFP_Rep_C_n_n),na.rm=TRUE),
  WT324=mean(c(FLAG_MSI2_324_WT_Rep_A_n_n,FLAG_MSI2_324_WT_Rep_B_n_n,FLAG_MSI2_324_WT_Rep_C_n_n),na.rm=TRUE),
  WT328=mean(c(FLAG_MSI2_328_WT_Rep_A_n_n,FLAG_MSI2_328_WT_Rep_B_n_n,FLAG_MSI2_328_WT_Rep_C_n_n),na.rm=TRUE),.groups='drop')
mat<-as.matrix(avg[,-1]);rownames(mat)<-avg$gene
mat<-mat[rowSums(is.na(mat))<ncol(mat),];z<-t(scale(t(mat)));z<-z[complete.cases(z),]
set.seed(123);km<-kmeans(z,centers=3)
cl<-tibble(source_symbol=names(km$cluster),cluster=as.integer(km$cluster))
save_tsv(left_join(cl,rename(avg,source_symbol=gene),by='source_symbol'),'reproduced_notebook_clusters')
save_tsv(bind_cols(tibble(cluster=1:3),as_tibble(km$centers)),'notebook_cluster_centers_scaled')
cls<-d %>% left_join(cl,by='source_symbol') %>% mutate(classification=case_when(
  cluster %in% c(2,3) & comparison_p>0 & comparison_p<0.5 & rank_stat>=1 ~ '324-specific',
  cluster %in% c(2,3) & comparison_p>0 & comparison_p<0.5 & rank_stat<= -1 ~ '328-specific',TRUE~'shared/non-specific'))
save_tsv(cls,'isoform_classification_proteins')
cg<-cls %>% filter(rank_included,mapped_to_hpa) %>% group_by(gene_id) %>% summarise(
  is324=any(classification=='324-specific'),is328=any(classification=='328-specific'),.groups='drop') %>% mutate(
    classification=case_when(is324&is328~'conflicting_duplicate_classes',is324~'324-specific',is328~'328-specific',TRUE~'shared/non-specific'))
save_tsv(cg,'isoform_classification_genes')
bg<-ranked$gene_id[ranked$mapped_to_hpa]
save_tsv(filter(ranked,mapped_to_hpa),'ora_measured_background')
save_tsv(count(cg,classification),'isoform_classification_counts')
collections<-list();gsea_all<-list();ora_all<-list();audit_all<-list();leading_all<-list()
for(conf in c('primary','high_confidence')) for(level in c('fine','collapsed')) {
  key<-paste(level,conf,sep='_')
  ann<-long %>% filter(confidence %in% if(conf=='primary') c('Enhanced','Supported','Approved') else c('Enhanced','Supported'))
  ann<-left_join(ann,hm,by='localization') %>% mutate(location=if(level=='fine') localization else collapsed_localization) %>% distinct(gene_id,location)
  paths<-split(ann$gene_id,ann$location); paths<-lapply(paths,sort)
  collections[[key]]<-paths
  save_tsv(ann,paste0('gene_sets_',key))
  writeLines(vapply(names(paths),function(n) paste(c(n,'HPA_25.1',paths[[n]]),collapse='\t'),character(1)),file.path(out,paste0('gene_sets_',key,'.gmt')))
  terms<-if(level=='fine') sort(lookup$`Subcellular location`) else sort(unique(hm$collapsed_localization))
  sizes<-tibble(localization=terms,hpa_size=vapply(terms,function(n) length(paths[[n]]),integer(1)),
    measured_size=vapply(terms,function(n) length(intersect(paths[[n]],names(stats))),integer(1))) %>% mutate(
    analysis=key,status=case_when(measured_size<5~'fewer_than_5_measured_genes',measured_size>=0.95*length(stats)~'at_least_95_percent_of_universe',TRUE~'tested'))
  audit_all[[key]]<-sizes
  eligible<-sizes$localization[sizes$status=='tested']
  set.seed(20260924)
  fg<-as_tibble(fgseaMultilevel(paths[eligible],stats,minSize=5,maxSize=length(stats)-1,eps=0,nPermSimple=10000,sampleSize=1001,BPPARAM=BiocParallel::SerialParam(),gseaParam=1))
  fg$padj<-p.adjust(fg$pval,'BH')
  fg<-fg %>% rename(localization=pathway,measured_size_fgsea=size) %>% mutate(direction=ifelse(NES>0,'324','328'),analysis=key,
    leading_edge_gene_ids=vapply(leadingEdge,join_text,character(1)),
    leading_edge_symbols=vapply(leadingEdge,function(g)join_text(ranked$gene_symbol[match(g,ranked$gene_id)]),character(1)),
    leading_edge_proteins=vapply(leadingEdge,function(g)join_text(ranked$protein_labels[match(g,ranked$gene_id)]),character(1)))
  stopifnot(!anyNA(fg$pval),!anyNA(fg$NES))
  le<-fg %>% select(analysis,localization,NES,padj,direction,leadingEdge) %>% unnest_longer(leadingEdge,values_to='gene_id') %>% left_join(ranked,by='gene_id')
  leading_all[[key]]<-le
  result<-left_join(sizes,select(fg,-analysis,-leadingEdge),by='localization')
  gsea_all[[key]]<-result
  save_tsv(result,paste0('gsea_',key))
  for(iso in c('324','328')) {
    hit<-cg$gene_id[cg$classification==paste0(iso,'-specific')];hit<-intersect(hit,bg)
    rr<-bind_rows(lapply(terms,function(term) {
      loc<-intersect(paths[[term]],bg);a<-length(intersect(hit,loc));b<-length(hit)-a;c<-length(loc)-a;dd<-length(bg)-a-b-c
      tested<-length(loc)>=5 && length(loc)<.95*length(bg) && length(hit)>0
      ft<-if(tested) fisher.test(matrix(c(a,b,c,dd),nrow=2,byrow=TRUE),alternative='greater') else NULL
      tibble(analysis=key,isoform=iso,localization=term,hit_in_location=a,total_hits=length(hit),background_in_location=length(loc),total_background=length(bg),
       odds_ratio=if(tested)unname(ft$estimate) else NA_real_,pval=if(tested)ft$p.value else NA_real_,status=if(tested)'tested' else 'insufficient_representation',
       contributing_gene_ids=join_text(intersect(hit,loc)),contributing_symbols=join_text(ranked$gene_symbol[match(intersect(hit,loc),ranked$gene_id)]),
       contributing_proteins=join_text(ranked$protein_labels[match(intersect(hit,loc),ranked$gene_id)]))
    })) %>% mutate(padj=p.adjust(pval,'BH'))
    ora_all[[paste(key,iso)]]<-rr
    save_tsv(rr,paste0('ora_',key,'_',iso))
  }
}
gsea<-bind_rows(gsea_all);ora<-bind_rows(ora_all);leading<-bind_rows(leading_all)
save_tsv(gsea,'gsea_all');save_tsv(ora,'ora_all');save_tsv(bind_rows(audit_all),'gene_set_eligibility')
save_tsv(leading,'leading_edges_all');save_tsv(filter(leading,padj<.05),'leading_edges_significant')
comparison<-gsea %>% mutate(confidence=ifelse(grepl('high_confidence',analysis),'high_confidence','primary'),level=ifelse(grepl('^fine',analysis),'fine','collapsed')) %>%
  select(level,localization,confidence,NES,padj,measured_size,status) %>% pivot_wider(names_from=confidence,values_from=c(NES,padj,measured_size,status)) %>%
  mutate(same_direction=sign(NES_primary)==sign(NES_high_confidence),significant_both=padj_primary<.05 & padj_high_confidence<.05)
save_tsv(comparison,'confidence_sensitivity_comparison')
theme_set(theme_minimal(base_size=12)+theme(panel.grid.minor=element_blank(),plot.title=element_text(face='bold'),plot.margin=margin(12,18,12,12)))
save_plot<-function(p,name,w=9,h=7) {
  ggsave(file.path(out,'figures',paste0(name,'.png')),p,width=w,height=h,dpi=180,bg='white')
  ggsave(file.path(out,'figures',paste0(name,'.pdf')),p,width=w,height=h,device=grDevices::pdf,bg='white',useDingbats=FALSE)
}
for(key in names(collections)) {
  x<-filter(gsea,analysis==key,status=='tested') %>% mutate(localization=factor(localization,levels=localization[order(NES)]),significant=padj<.05)
  p<-ggplot(x,aes(NES,localization))+geom_vline(xintercept=0,colour='grey65',linetype=2)+
    geom_point(aes(size=measured_size,colour=padj,shape=significant),stroke=1.1)+
    scale_color_viridis_c(option='C',direction=-1,limits=c(0,1),name='BH FDR')+
    scale_shape_manual(values=c('FALSE'=16,'TRUE'=8),name='FDR < 0.05')+scale_size_continuous(range=c(2.5,8),name='Measured genes')+
    labs(title=paste('HPA localization enrichment:',gsub('_',' ',key)),subtitle='All ranked genes; positive NES favors MSI2 324',x='Normalized enrichment score (negative: 328; positive: 324)',y=NULL)
  save_plot(p,paste0('NES_',key),10,max(5,0.28*nrow(x)+1.8))
}
# Standard running enrichment curves: top two FDR-ranked terms per direction at each primary level.
interest<-c('Nucleoplasm','Nuclear speckles','Nuclear bodies','Nucleoli','Nuclear membrane','Cytosol','Endoplasmic reticulum','Golgi apparatus',
 'Vesicles','Endosomes','Lysosomes','Plasma membrane','Cell Junctions','Actin filaments','Microtubules','Centrosome','Mitochondria','Intermediate filaments',
 'Actin-filament associated','Microtubule associated','Vesicular system','Plasma membrane/cell junctions')
chosen<-gsea %>% filter(grepl('_primary$',analysis),status=='tested',localization %in% interest) %>% group_by(analysis,direction) %>% arrange(pval,desc(abs(NES)),.by_group=TRUE) %>% slice_head(n=2) %>% ungroup()
save_tsv(chosen,'running_enrichment_plot_selection')
# Remove only superseded running plots produced in this analysis directory.
expected_running<-paste0('running_',chosen$analysis,'_',gsub('[^A-Za-z0-9]+','_',chosen$localization))
old_running<-list.files(file.path(out,'figures'),pattern='^running_.*\\.(png|pdf)$',full.names=TRUE)
unlink(old_running[!sub('\\.(png|pdf)$','',basename(old_running)) %in% expected_running])
for(i in seq_len(nrow(chosen))) {
  r<-chosen[i,]; p<-plotEnrichment(collections[[r$analysis]][[r$localization]],stats)+
    labs(title=paste(r$localization,'-',r$analysis),subtitle=sprintf('NES %.2f; BH FDR %.3g; direction %s%s',r$NES,r$padj,r$direction,ifelse(r$padj<.05,'',' (not significant)')),x='Rank: 324-associated (left) to 328-associated (right)',y='Running enrichment score')
  save_plot(p,paste0('running_',r$analysis,'_',gsub('[^A-Za-z0-9]+','_',r$localization)),9,5)
}
for(level in c('fine','collapsed')) {
  tile<-bind_rows(
    gsea %>% filter(analysis==paste0(level,'_primary')) %>% transmute(localization,display='GSEA primary',score=sign(NES)*-log10(padj),padj),
    gsea %>% filter(analysis==paste0(level,'_high_confidence')) %>% transmute(localization,display='GSEA high confidence',score=sign(NES)*-log10(padj),padj),
    ora %>% filter(analysis==paste0(level,'_primary')) %>% transmute(localization,display=paste0('ORA ',isoform),score=ifelse(isoform=='324',1,-1)*-log10(padj),padj))
  order<-gsea %>% filter(analysis==paste0(level,'_primary'),status=='tested') %>% arrange(NES) %>% pull(localization)
  tile<-tile %>% filter(localization %in% order) %>% mutate(localization=factor(localization,levels=order),display=factor(display,levels=c('GSEA primary','GSEA high confidence','ORA 324','ORA 328')))
  save_tsv(tile,paste0('analysis_tile_data_',level))
  lim<-max(abs(tile$score),na.rm=TRUE)
  p<-ggplot(tile,aes(display,localization,fill=score))+geom_tile(colour='white')+geom_text(aes(label=ifelse(!is.na(padj)&padj<.05,'*','')),size=5)+
    scale_fill_gradient2(low='#2864a0',mid='white',high='#bc3d39',limits=c(-lim,lim),na.value='grey85',name='Signed -log10(FDR)')+
    labs(title=paste('Localization evidence:',level),subtitle='Red: 324; blue: 328. * FDR < 0.05. Grey: not tested.',x=NULL,y=NULL)+theme(axis.text.x=element_text(angle=20,hjust=1))
  save_plot(p,paste0('analysis_tile_',level),9,max(5,.27*length(order)+2))
}
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
saveRDS(list(stats=stats,collections=collections,gsea=gsea,ora=ora,ranked=ranked,mapping=d),file.path(out,'analysis_objects.rds'))
stopifnot(!anyDuplicated(names(stats)),length(stats)==nrow(ranked),all(is.finite(stats)),all(leading$gene_id %in% names(stats)))
save_tsv(tibble(check=c('unique_gene_ranks','finite_rank_values','all_49_locations_mapped','all_GSEA_pvalues_finite','ORA_hits_subset_of_measured_HPA_background'),passed=TRUE),'validation_checks')
print(qc);print(filter(gsea,padj<.05) %>% select(analysis,localization,measured_size,NES,padj))
print(filter(ora,padj<.05) %>% select(analysis,isoform,localization,odds_ratio,padj))
