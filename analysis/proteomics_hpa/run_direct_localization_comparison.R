#!/usr/bin/env Rscript
# Compare localization membership between the existing 324 and 328 sets.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(readr);library(ggplot2)})
base<-'output/proteomics_hpa/broad_exploratory_ora'
out<-'output/proteomics_hpa/direct_localization_comparison'
dir.create(file.path(out,'figures'),recursive=TRUE,showWarnings=FALSE)
read<-function(n)read_tsv(file.path(base,n),col_types=cols(.default=col_character()),show_col_types=FALSE)
save<-function(z,n)write_tsv(as.data.frame(z),file.path(out,paste0(n,'.tsv')),na='NA')
jt<-function(z)paste(sort(unique(z[!is.na(z)&nzchar(z)])),collapse=';')
input_paths<-file.path(base,c('unchanged_hit_memberships.tsv','measured_gene_broad_membership.tsv','fine_to_broad_mapping.tsv'))
hash_before<-tools::md5sum(input_paths)
hits<-read('unchanged_hit_memberships.tsv');membership<-read('measured_gene_broad_membership.tsv');mapping<-read('fine_to_broad_mapping.tsv')
terms<-c('Nuclear/nucleoplasmic','Nucleolus','Cytosolic','Mitochondrial','ER/Golgi','Vesicular/endolysosomal','Plasma membrane/junctional','Cytoskeletal')
stopifnot(setequal(na.omit(mapping$broad_category),terms))
pairs<-c('cluster','hits','specific');rr<-list();members<-list();contrib<-list()
for(pair in pairs) {
 h324<-unique(hits$ENTREZID[hits$set_name==paste0('genes_324_',pair)])
 h328<-unique(hits$ENTREZID[hits$set_name==paste0('genes_328_',pair)])
 shared<-intersect(h324,h328);u324<-setdiff(h324,shared);u328<-setdiff(h328,shared)
 members[[pair]]<-tibble(pair=pair,ENTREZID=union(h324,h328)) %>% mutate(in_324=ENTREZID %in% h324,in_328=ENTREZID %in% h328,
  comparison_group=case_when(in_324&in_328~'shared_fixed',in_324~'324_only',TRUE~'328_only')) %>%
  left_join(hits %>% select(ENTREZID,source_symbols,hpa_gene_id) %>% distinct(),by='ENTREZID')
 rr[[pair]]<-bind_rows(lapply(terms,function(term) {
  loc<-unique(membership$ENTREZID[membership$broad_category==term])
  a<-length(intersect(u324,loc));b<-length(u324)-a;c<-length(intersect(u328,loc));d<-length(u328)-c
  full_a<-length(intersect(h324,loc));full_c<-length(intersect(h328,loc))
  ft<-fisher.test(matrix(c(a,b,c,d),2,byrow=TRUE),alternative='two.sided',conf.level=.95)
  symbols<-function(ids)jt(hits$source_symbols[match(ids,hits$ENTREZID)])
  contrib[[paste(pair,term)]]<<-members[[pair]] %>% filter(ENTREZID %in% loc) %>% mutate(localization=term)
  tibble(pair=pair,localization=term,full_324_in=full_a,full_324_out=length(h324)-full_a,full_324_total=length(h324),
   full_328_in=full_c,full_328_out=length(h328)-full_c,full_328_total=length(h328),
   full_324_proportion=full_a/length(h324),full_328_proportion=full_c/length(h328),
   full_difference_percentage_points=100*(full_a/length(h324)-full_c/length(h328)),
   shared_total=length(shared),shared_in_location=length(intersect(shared,loc)),shared_genes=symbols(shared),
   test_324_in=a,test_324_out=b,test_324_total=length(u324),test_328_in=c,test_328_out=d,test_328_total=length(u328),
   test_324_proportion=a/length(u324),test_328_proportion=c/length(u328),
   test_difference_percentage_points=100*(a/length(u324)-c/length(u328)),
   odds_ratio_324_vs_328=unname(ft$estimate),OR_CI95_lower=ft$conf.int[1],OR_CI95_upper=ft$conf.int[2],
   p_value=ft$p.value,test_method=ifelse(length(shared)>0,'Two-sided Fisher conditional on shared genes fixed; test nonshared genes','Two-sided Fisher; disjoint original sets'),
   genes_324_in_location=symbols(intersect(h324,loc)),genes_328_in_location=symbols(intersect(h328,loc)),
   test_genes_324_in_location=symbols(intersect(u324,loc)),test_genes_328_in_location=symbols(intersect(u328,loc)))
 })) %>% mutate(BH_FDR=p.adjust(p_value,'BH'))
 save(rr[[pair]],paste0('direct_',pair))
}
results<-bind_rows(rr) %>% mutate(BH_FDR_all_24=p.adjust(p_value,'BH'))
save(results,'direct_comparison_all');save(bind_rows(members),'original_membership_and_shared_gene_audit');save(bind_rows(contrib),'contributing_genes_long')
save(filter(results,localization %in% c('Nucleolus','Mitochondrial')),'nucleolus_mitochondria_comparisons')
summary<-results %>% group_by(pair) %>% summarise(full_324_total=first(full_324_total),full_328_total=first(full_328_total),shared=first(shared_total),
 test_324_total=first(test_324_total),test_328_total=first(test_328_total),tests=n(),significant_005=sum(BH_FDR<.05),significant_010=sum(BH_FDR<.1),
 minimum_p=min(p_value),minimum_FDR=min(BH_FDR),.groups='drop')
save(summary,'analysis_summary');save(filter(results,BH_FDR<.05),'significant_FDR005')
labels<-c(cluster='Full clusters',hits='EGFP-filtered clusters',specific='Isoform-specific subsets')
theme_set(theme_minimal(base_size=12)+theme(panel.grid.minor=element_blank(),plot.title=element_text(face='bold')))
for(pair_name in pairs) {
 z<-filter(results,pair==pair_name) %>% mutate(localization=factor(localization,levels=rev(terms)))
 p<-ggplot(z,aes(y=localization))+geom_segment(aes(x=100*full_324_proportion,xend=100*full_328_proportion,yend=localization),colour='grey65')+
  geom_point(aes(x=100*full_324_proportion,colour='324'),size=3)+geom_point(aes(x=100*full_328_proportion,colour='328'),size=3)+
  scale_colour_manual(values=c('324'='#b14449','328'='#2673a3'),name='Original hit set')+
  labs(title=paste('Direct localization comparison:',labels[pair_name]),subtitle=paste0('Original-set proportions; two-sided tests across 8 categories. Lowest FDR = ',sprintf('%.3f',min(z$BH_FDR))),
   x='Genes annotated to category (% of original hit set)',y=NULL)+theme(legend.position='bottom')
 for(ext in c('png','pdf')) {
  path<-file.path(out,'figures',paste0('proportions_',pair_name,'.',ext))
  if(ext=='png')ggsave(path,p,width=10,height=5.5,dpi=180,bg='white') else ggsave(path,p,width=10,height=5.5,device=grDevices::pdf,useDingbats=FALSE)
 }
}
# Independent two-sided conditional hypergeometric p-values, matching Fisher's definition.
manual<-vapply(seq_len(nrow(results)),function(i) {
 r<-results[i,];M<-r$test_324_in+r$test_328_in;N<-r$test_324_out+r$test_328_out;k<-r$test_324_total
 possible<-seq.int(max(0,k-N),min(k,M));pr<-dhyper(possible,M,N,k);observed<-dhyper(r$test_324_in,M,N,k)
 sum(pr[pr<=observed*(1+1e-7)])
},numeric(1))
checks<-tibble(check=c('Frozen input files unchanged','24 comparisons and eight per pair','No shared gene counted as independent observations in both test rows',
 'Independent two-sided hypergeometric p-values','BH independently applied per pair','All CI directions consistent with odds ratios'),
 passed=c(identical(hash_before,tools::md5sum(input_paths)),nrow(results)==24&&all(table(results$pair)==8),
 all(results$full_324_total-results$test_324_total==results$shared_total)&all(results$full_328_total-results$test_328_total==results$shared_total),
 max(abs(manual-results$p_value))<1e-10,all(vapply(split(results,results$pair),function(z)all(abs(p.adjust(z$p_value,'BH')-z$BH_FDR)<1e-12),logical(1))),
 all(results$OR_CI95_lower<=results$odds_ratio_324_vs_328&results$OR_CI95_upper>=results$odds_ratio_324_vs_328)))
stopifnot(all(checks$passed));save(checks,'validation_checks');save(tibble(path=names(hash_before),md5=unname(hash_before)),'input_checksums')
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
print(summary);print(results %>% arrange(BH_FDR,p_value) %>% select(pair,localization,test_324_in,test_328_in,odds_ratio_324_vs_328,p_value,BH_FDR) %>% head(12))
