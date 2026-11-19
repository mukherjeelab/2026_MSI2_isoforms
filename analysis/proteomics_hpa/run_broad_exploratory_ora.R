#!/usr/bin/env Rscript
# Fixed eight-category exploratory follow-up. Never reselect Reactome hit genes.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(readr);library(ggplot2)})
base<-'output/proteomics_hpa'
prior<-file.path(base,'reactome_matched_ora')
out<-file.path(base,'broad_exploratory_ora')
dir.create(file.path(out,'figures'),recursive=TRUE,showWarnings=FALSE)
read<-function(path)read_tsv(path,col_types=cols(.default=col_character()),show_col_types=FALSE)
save<-function(z,n)write_tsv(as.data.frame(z),file.path(out,paste0(n,'.tsv')),na='NA')
jt<-function(z)paste(sort(unique(z[!is.na(z)&nzchar(z)])),collapse=';')
input_paths<-file.path(prior,c('exact_mapped_hit_memberships.tsv','measured_universe_hpa_mapping.tsv','measured_gene_localization_membership.tsv'))
hash_before<-tools::md5sum(input_paths)
hits<-read(input_paths[1]);universe<-read(input_paths[2]);ann<-read(input_paths[3])
bg<-sort(unique(universe$ENTREZID));stopifnot(length(bg)==1638)
hierarchy<-read(file.path(base,'collapsed_location_mapping.tsv'))
groups<-c('Nucleoplasm'='Nuclear/nucleoplasmic','Nuclear membrane'='Nuclear/nucleoplasmic','Nucleoli'='Nucleolus',
 'Cytosol'='Cytosolic','Mitochondria'='Mitochondrial','Endoplasmic reticulum'='ER/Golgi','Golgi apparatus'='ER/Golgi',
 'Vesicular system'='Vesicular/endolysosomal','Plasma membrane/cell junctions'='Plasma membrane/junctional',
 'Actin-filament associated'='Cytoskeletal','Microtubule associated'='Cytoskeletal','Intermediate filaments'='Cytoskeletal',
 'Centrosome'='Cytoskeletal','Primary cilium'='Cytoskeletal')
terms<-c('Nuclear/nucleoplasmic','Nucleolus','Cytosolic','Mitochondrial','ER/Golgi','Vesicular/endolysosomal','Plasma membrane/junctional','Cytoskeletal')
mapping<-hierarchy %>% mutate(broad_category=unname(groups[collapsed_localization]),
 included=!is.na(broad_category),decision=ifelse(included,'Included in fixed eight-category scheme','Sperm-specific term outside the requested broad scheme'))
stopifnot(nrow(mapping)==49,setequal(na.omit(mapping$broad_category),terms))
# Save the mapping before conducting any enrichment tests.
save(mapping,'fine_to_broad_mapping')
membership<-ann %>% left_join(select(mapping,localization,broad_category),by='localization') %>% filter(!is.na(broad_category)) %>%
 group_by(ENTREZID,hpa_gene_id,hpa_symbol,broad_category) %>% summarise(contributing_fine_terms=jt(localization),.groups='drop')
save(membership,'measured_gene_broad_membership')
save(hits,'unchanged_hit_memberships');save(universe,'unchanged_measured_background')
hpa_long<-read(file.path(base,'hpa_gene_location_long.tsv')) %>% filter(confidence %in% c('Enhanced','Supported','Approved')) %>%
 left_join(select(mapping,localization,broad_category),by='localization') %>% filter(!is.na(broad_category)) %>% distinct(gene_id,broad_category)
save(hpa_long,'hpa_broad_gene_sets')
results<-list();contributors<-list()
for(set in unique(hits$set_name)) {
 hit<-sort(unique(hits$ENTREZID[hits$set_name==set]));iso<-strsplit(set,'_')[[1]][2];pair<-sub('^genes_[0-9]+_','',set)
 stopifnot(all(hit %in% bg))
 rr<-bind_rows(lapply(terms,function(term) {
  loc<-unique(membership$ENTREZID[membership$broad_category==term]);genes<-intersect(hit,loc)
  a<-length(genes);b<-length(hit)-a;c<-length(loc)-a;dd<-length(bg)-a-b-c
  stopifnot(length(loc)>=5,length(loc)<.95*length(bg))
  ft<-fisher.test(matrix(c(a,b,c,dd),2,byrow=TRUE),alternative='greater')
  contributors[[paste(set,term)]]<<-membership %>% filter(ENTREZID %in% genes,broad_category==term) %>%
    left_join(select(universe,ENTREZID,source_symbols),by='ENTREZID') %>% mutate(set_name=set,pair=pair,isoform=iso)
  tibble(pair=pair,set_name=set,isoform=iso,localization=term,hpa_gene_set_size=sum(hpa_long$broad_category==term),
   overlap_count=a,hit_set_size=length(hit),background_overlap=length(loc),background_size=length(bg),
   GeneRatio=paste0(a,'/',length(hit)),gene_ratio=a/length(hit),BgRatio=paste0(length(loc),'/',length(bg)),background_ratio=length(loc)/length(bg),
   odds_ratio=unname(ft$estimate),p_value=ft$p.value,contributing_genes=jt(universe$source_symbols[match(genes,universe$ENTREZID)]),
   contributing_entrez=jt(genes),contributing_hpa_genes=jt(universe$hpa_gene_id[match(genes,universe$ENTREZID)]))
 })) %>% mutate(BH_FDR=p.adjust(p_value,'BH'))
 stopifnot(nrow(rr)==8)
 results[[set]]<-rr;save(rr,paste0('ora_',set))
}
all_results<-bind_rows(results)
save(all_results,'broad_hpa_ora_all');save(bind_rows(contributors),'contributing_genes_long')
save(filter(all_results,BH_FDR<.05),'significant_FDR005');save(filter(all_results,BH_FDR<.1),'significant_FDR010')
summary<-all_results %>% group_by(pair,isoform) %>% summarise(hit_set_size=first(hit_set_size),background_size=first(background_size),tests=n(),
 significant_005=sum(BH_FDR<.05),significant_010=sum(BH_FDR<.1),minimum_p=min(p_value),minimum_FDR=min(BH_FDR),.groups='drop')
save(summary,'analysis_summary')
save(all_results %>% filter(localization %in% c('Nucleolus','Mitochondrial')),'nucleolus_mitochondria_hypotheses')
theme_set(theme_minimal(base_size=12)+theme(panel.grid.minor=element_blank(),plot.title=element_text(face='bold')))
pair_labels<-c(cluster='Full clusters',hits='EGFP-filtered cluster subsets',specific='Isoform-specific subsets')
for(pair_name in names(pair_labels)) {
 z<-filter(all_results,pair==pair_name) %>% mutate(localization=factor(localization,levels=rev(terms)),isoform=paste('MSI2',isoform))
 p<-ggplot(z,aes(gene_ratio,localization))+geom_point(aes(size=overlap_count,colour=BH_FDR,shape=BH_FDR<.05),stroke=1)+
  facet_wrap(~isoform,nrow=1)+scale_colour_viridis_c(option='C',direction=-1,limits=c(0,1),name='BH FDR')+
  scale_shape_manual(values=c('FALSE'=16,'TRUE'=8),name='FDR < 0.05')+scale_size_continuous(range=c(3,8),name='Overlap')+
  labs(title=paste('Exploratory broad HPA ORA:',pair_labels[pair_name]),subtitle='Eight tests per original hit set; 1,638 measured genes in background',x='GeneRatio (overlap / original mapped hit-set size)',y=NULL)
 ggsave(file.path(out,'figures',paste0('broad_ORA_',pair_name,'.png')),p,width=12,height=5.5,dpi=180,bg='white')
 ggsave(file.path(out,'figures',paste0('broad_ORA_',pair_name,'.pdf')),p,width=12,height=5.5,device=grDevices::pdf,useDingbats=FALSE)
}
hp<-with(all_results,phyper(overlap_count-1,background_overlap,background_size-background_overlap,hit_set_size,lower.tail=FALSE))
prior_results<-read_tsv(file.path(prior,'hpa_ora_all_exact_reactome_sets.tsv'),show_col_types=FALSE) %>% filter(resolution=='collapsed')
unchanged<-all_results %>% filter(localization %in% c('Nucleolus','Mitochondrial')) %>% mutate(old_location=ifelse(localization=='Nucleolus','Nucleoli','Mitochondria')) %>%
 left_join(select(prior_results,set_name,old_location=localization,old_overlap=overlap_count,old_background=background_overlap,old_p=p_value),by=c('set_name','old_location'))
checks<-tibble(check=c('Original input files unchanged','48 tests: eight categories for each of six original sets','Fisher matches hypergeometric tail',
 'BH applied independently across eight tests per original set','Nucleolus and mitochondrial counts and nominal p-values unchanged','Multilabel membership preserved'),
 passed=c(identical(hash_before,tools::md5sum(input_paths)),nrow(all_results)==48&&all(table(all_results$set_name)==8),max(abs(hp-all_results$p_value))<1e-10,
 all(vapply(split(all_results,all_results$set_name),function(z)all(abs(z$BH_FDR-p.adjust(z$p_value,'BH'))<1e-12),logical(1))),
 all(unchanged$overlap_count==unchanged$old_overlap)&all(unchanged$background_overlap==unchanged$old_background)&all(abs(unchanged$p_value-unchanged$old_p)<1e-12),
 any(table(membership$ENTREZID)>1)))
stopifnot(all(checks$passed));save(checks,'validation_checks')
save(tibble(path=names(hash_before),md5=unname(hash_before)),'input_checksums')
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
print(summary);print(all_results %>% arrange(BH_FDR,p_value) %>% select(pair,isoform,localization,GeneRatio,BgRatio,odds_ratio,p_value,BH_FDR) %>% head(12))
