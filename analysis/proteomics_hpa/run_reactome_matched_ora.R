#!/usr/bin/env Rscript
# Exact notebook inputs reconstructed by audit_reactome_inputs.R; no hit reclassification.
suppressPackageStartupMessages({library(org.Hs.eg.db);library(dplyr);library(tidyr);library(readr);library(ggplot2)})
base <- 'output/proteomics_hpa'
out <- file.path(base,'reactome_matched_ora')
dir.create(file.path(out,'figures'),recursive=TRUE,showWarnings=FALSE)
save_tsv<-function(x,n)write_tsv(as.data.frame(x),file.path(out,paste0(n,'.tsv')),na='NA')
jt<-function(x)paste(sort(unique(x[!is.na(x)&nzchar(x)])),collapse=';')
inputs<-readRDS(file.path(base,'reactome_input_audit/original_reactome_inputs.rds'))
bg<-sort(unique(inputs$entrez_mapping$ENTREZID))
stopifnot(length(bg)==1638)
atlas<-read_tsv(unz('data/hpa_25.1/proteinatlas.tsv.zip','proteinatlas.tsv'),show_col_types=FALSE)
long<-read_tsv(file.path(base,'hpa_gene_location_long.tsv'),show_col_types=FALSE) %>% filter(confidence %in% c('Enhanced','Supported','Approved'))
hierarchy<-read_tsv(file.path(base,'collapsed_location_mapping.tsv'),show_col_types=FALSE)
# Bridge to the exact Entrez identifiers used by Reactome. Preserve unknowns in universe.
bridge<-AnnotationDbi::select(org.Hs.eg.db,keys=bg,keytype='ENTREZID',columns=c('ENSEMBL','SYMBOL')) %>% as_tibble() %>% distinct()
save_tsv(bridge,'orgdb_entrez_ensembl_bridge')
mapping<-bind_rows(lapply(bg,function(id) {
  b<-filter(bridge,ENTREZID==id)
  ec<-intersect(b$ENSEMBL,atlas$Ensembl)
  syms<-unique(inputs$entrez_mapping$SYMBOL[inputs$entrez_mapping$ENTREZID==id])
  sc<-unique(atlas$Ensembl[atlas$Gene %in% c(syms,b$SYMBOL)])
  gene<-NA_character_;method<-'unmapped_to_hpa'
  if(length(ec)==1 && (length(sc)==0||ec %in% sc)) {gene<-ec;method<-'orgdb_ensembl'}
  else if(length(ec)>1&&length(intersect(ec,sc))==1) {gene<-intersect(ec,sc);method<-'ensembl_resolved_by_symbol'}
  else if(length(ec)==0&&length(sc)==1) {gene<-sc;method<-'exact_symbol'}
  else if(length(ec)+length(sc)>0) method<-'ambiguous_or_conflicting'
  tibble(ENTREZID=id,source_symbols=jt(syms),orgdb_symbols=jt(b$SYMBOL),hpa_ensembl_candidates=jt(ec),hpa_symbol_candidates=jt(sc),
         hpa_gene_id=gene,hpa_symbol=atlas$Gene[match(gene,atlas$Ensembl)],mapping_method=method)
}))
# Do not assign one HPA gene's annotations to multiple distinct Entrez genes.
ambiguous_reverse<-mapping %>% filter(!is.na(hpa_gene_id)) %>% count(hpa_gene_id) %>% filter(n>1) %>% pull(hpa_gene_id)
mapping<-mapping %>% mutate(ambiguous_reverse_mapping=hpa_gene_id %in% ambiguous_reverse,
  annotation_usable=!is.na(hpa_gene_id)&!ambiguous_reverse_mapping,
  has_primary_localization=annotation_usable & hpa_gene_id %in% long$gene_id)
save_tsv(mapping,'measured_universe_hpa_mapping')
save_tsv(filter(mapping,!annotation_usable|!has_primary_localization),'unannotated_or_unmapped_universe_genes')
ann<-mapping %>% filter(annotation_usable) %>% select(ENTREZID,hpa_gene_id) %>%
  inner_join(long,by=c('hpa_gene_id'='gene_id'),relationship='many-to-many') %>%
  left_join(hierarchy,by='localization') %>% distinct(ENTREZID,hpa_gene_id,hpa_symbol=gene_symbol,localization,collapsed_localization)
save_tsv(ann,'measured_gene_localization_membership')
sets<-lapply(inputs$sets,function(symbols)sort(unique(inputs$entrez_mapping$ENTREZID[inputs$entrez_mapping$SYMBOL %in% symbols])))
set_members<-bind_rows(lapply(names(sets),function(n)tibble(set_name=n,ENTREZID=sets[[n]]))) %>% left_join(mapping,by='ENTREZID')
save_tsv(set_members,'exact_mapped_hit_memberships')
save_tsv(bind_rows(lapply(names(inputs$sets),function(n)tibble(set_name=n,SYMBOL=inputs$sets[[n]]))) %>%
  left_join(inputs$entrez_mapping,by='SYMBOL',relationship='many-to-many'),'original_symbol_hit_memberships')
save_tsv(tibble(SYMBOL=inputs$universe_genes) %>% left_join(inputs$entrez_mapping,by='SYMBOL'),'original_measured_universe')
qc<-set_members %>% group_by(set_name) %>% summarise(mapped_hit_genes=n(),hpa_gene_mapped=sum(annotation_usable),primary_localized=sum(has_primary_localization),.groups='drop')
save_tsv(qc,'hit_mapping_QC')
save_tsv(tibble(metric=c('original_universe_symbols','entrez_mapped_universe','unique_HPA_gene_mapped','primary_localized','no_primary_localization'),
  n=c(length(inputs$universe_genes),length(bg),sum(mapping$annotation_usable),sum(mapping$has_primary_localization),sum(!mapping$has_primary_localization))),'universe_QC')
pair_names<-c(cluster='Full clusters',specific='Isoform-specific subsets',hits='EGFP-filtered cluster subsets')
overlap<-bind_rows(lapply(names(pair_names),function(pair) {
  z<-intersect(sets[[paste0('genes_324_',pair)]],sets[[paste0('genes_328_',pair)]])
  tibble(pair=pair,n_shared_between_hit_sets=length(z),shared_entrez=jt(z),shared_symbols=jt(mapping$source_symbols[match(z,mapping$ENTREZID)]))
}))
save_tsv(overlap,'between_set_overlap_preserved')
results<-list()
for(pair in names(pair_names)) for(level in c('fine','collapsed')) for(iso in c('324','328')) {
  set_name<-paste0('genes_',iso,'_',pair);hit<-sets[[set_name]]
  term_col<-if(level=='fine')'localization' else 'collapsed_localization'
  terms<-sort(unique(hierarchy[[term_col]]))
  full_ann<-long %>% left_join(hierarchy,by='localization')
  rr<-bind_rows(lapply(terms,function(term) {
    loc<-unique(ann$ENTREZID[ann[[term_col]]==term]);genes<-intersect(hit,loc)
    a<-length(genes);b<-length(hit)-a;c<-length(loc)-a;dd<-length(bg)-a-b-c
    tested<-length(loc)>=5&&length(loc)<.95*length(bg)
    ft<-if(tested)fisher.test(matrix(c(a,b,c,dd),2,byrow=TRUE),alternative='greater') else NULL
    tibble(pair=pair,pair_description=unname(pair_names[pair]),set_name=set_name,isoform=iso,resolution=level,localization=term,
      hpa_gene_set_size=n_distinct(full_ann$gene_id[full_ann[[term_col]]==term]),
      overlap_count=a,hit_set_size=length(hit),background_overlap=length(loc),background_size=length(bg),
      GeneRatio=paste0(a,'/',length(hit)),gene_ratio=a/length(hit),BgRatio=paste0(length(loc),'/',length(bg)),background_ratio=length(loc)/length(bg),
      odds_ratio=if(tested)unname(ft$estimate) else NA_real_,p_value=if(tested)ft$p.value else NA_real_,
      status=case_when(length(loc)<5~'fewer_than_5_measured_genes',length(loc)>=.95*length(bg)~'at_least_95_percent_of_background',TRUE~'tested'),
      contributing_entrez=jt(genes),contributing_genes=jt(mapping$source_symbols[match(genes,mapping$ENTREZID)]),
      contributing_hpa_genes=jt(mapping$hpa_gene_id[match(genes,mapping$ENTREZID)]))
  })) %>% mutate(BH_FDR=p.adjust(p_value,'BH'))
  results[[paste(pair,level,iso,sep='_')]]<-rr
}
all<-bind_rows(results) %>% group_by(pair,isoform) %>% mutate(BH_FDR_both_resolutions=p.adjust(p_value,'BH')) %>% ungroup()
save_tsv(all,'hpa_ora_all_exact_reactome_sets')
for(n in names(results)) {
  rr<-results[[n]]
  save_tsv(filter(all,pair==rr$pair[1],resolution==rr$resolution[1],isoform==rr$isoform[1]),paste0('ora_',n))
}
summary<-all %>% group_by(pair,resolution,isoform) %>% summarise(hit_set_size=first(hit_set_size),background_size=first(background_size),
  tested=sum(status=='tested'),significant_005=sum(BH_FDR<.05,na.rm=TRUE),significant_010=sum(BH_FDR<.1,na.rm=TRUE),
  minimum_p=min(p_value,na.rm=TRUE),minimum_FDR=min(BH_FDR,na.rm=TRUE),.groups='drop')
save_tsv(summary,'analysis_summary');save_tsv(filter(all,BH_FDR<.05),'significant_localizations')
theme_set(theme_minimal(base_size=11)+theme(panel.grid.minor=element_blank(),plot.title=element_text(face='bold')))
for(pair_name in names(pair_names)) for(level in c('fine','collapsed')) {
  z<-filter(all,pair==pair_name,resolution==level,status=='tested') %>% mutate(isoform=paste('MSI2',isoform))
  p<-ggplot(z,aes(gene_ratio,localization))+geom_point(aes(size=overlap_count,colour=BH_FDR,shape=BH_FDR<.05),stroke=1)+
    facet_wrap(~isoform,nrow=1)+scale_colour_viridis_c(option='C',direction=-1,limits=c(0,1),name='BH FDR')+
    scale_shape_manual(values=c('FALSE'=16,'TRUE'=8),name='FDR < 0.05')+scale_size_continuous(range=c(2,7),name='Overlap')+
    labs(title=paste(unname(pair_names[pair_name]),'-',level,'HPA'),subtitle='Original Reactome input sets; 1,638 measured Entrez genes in background',x='GeneRatio (overlap / original mapped hit set)',y=NULL)
  name<-paste0('ORA_',pair_name,'_',level);height<-max(5,.27*n_distinct(z$localization)+2)
  ggsave(file.path(out,'figures',paste0(name,'.png')),p,width=12,height=height,dpi=180,bg='white')
  ggsave(file.path(out,'figures',paste0(name,'.pdf')),p,width=12,height=height,device=grDevices::pdf,useDingbats=FALSE)
}
# Validate exact mapped memberships, counts, hypergeometric equivalence and BH scopes.
t<-filter(all,status=='tested')
hp<-with(t,phyper(overlap_count-1,background_overlap,background_size-background_overlap,hit_set_size,lower.tail=FALSE))
checks<-tibble(check=c('All six hit sets are exact original notebook Entrez mappings','All hits are inside the original Entrez universe',
 'Every result uses the original 1638-gene background','Fisher p-values match hypergeometric upper tails','Per-set per-resolution BH verified','Multi-localization membership retained'),
 passed=c(all(vapply(names(sets),function(n)setequal(sets[[n]],inputs$entrez_mapping$ENTREZID[inputs$entrez_mapping$SYMBOL %in% inputs$sets[[n]]]),logical(1))),
 all(unlist(sets) %in% bg),all(all$background_size==length(bg)),max(abs(hp-t$p_value))<1e-10,
 all(vapply(split(all,interaction(all$pair,all$resolution,all$isoform,drop=TRUE)),function(z)all(abs(p.adjust(z$p_value,'BH')-z$BH_FDR)<1e-12,na.rm=TRUE),logical(1))),
 any(table(ann$ENTREZID)>1)))
stopifnot(all(checks$passed));save_tsv(checks,'validation_checks')
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
save_tsv(tibble(file=c('code/Proetomics.Rmd',file.path(base,'reactome_input_audit/original_reactome_inputs.rds'),'data/hpa_25.1/proteinatlas.tsv.zip'),
 md5=unname(tools::md5sum(c('code/Proetomics.Rmd',file.path(base,'reactome_input_audit/original_reactome_inputs.rds'),'data/hpa_25.1/proteinatlas.tsv.zip')))),'input_checksums')
print(summary);print(qc);print(filter(all,BH_FDR<.1) %>% select(pair,resolution,isoform,localization,GeneRatio,BgRatio,odds_ratio,p_value,BH_FDR))
