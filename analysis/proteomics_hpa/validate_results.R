suppressPackageStartupMessages({library(dplyr);library(readr)})
o<-'output/proteomics_hpa'
read<-function(n)read_tsv(file.path(o,paste0(n,'.tsv')),show_col_types=FALSE)
x<-readRDS(file.path(o,'analysis_objects.rds'))
m<-x$mapping;r<-x$ranked;g<-x$gsea;ora<-x$ora
checks<-list()
record<-function(name,passed,detail) {stopifnot(passed);checks[[name]]<<-tibble(check=name,passed=passed,detail=detail)}
record('No comparison significance filter',sum(m$testable)==1677 && sum(m$rank_included)==1654,'1677 tested rows; only 23 unresolved identifier rows excluded')
record('Rank direction and duplicate median',all(vapply(seq_len(nrow(r)),function(i)abs(r$rank_stat[i]-median(-log2(m$fc_328_324[m$gene_id %in% r$gene_id[i]&m$testable])))<1e-12,logical(1))), 'Independent median of negative log2 source 328/324 ratios')
record('Gene rank uniqueness',!anyDuplicated(r$gene_id)&&nrow(r)==1641,'1641 unique genes')
errs<-numeric()
for(i in which(g$status=='tested')) {
  ids<-x$collections[[g$analysis[i]]][[g$localization[i]]]
  hit<-names(x$stats) %in% ids
  inc<-ifelse(hit,abs(x$stats)/sum(abs(x$stats[hit])),-1/sum(!hit))
  run<-c(0,cumsum(inc));es<-if(max(run)> -min(run))max(run) else min(run)
  errs<-c(errs,abs(es-g$ES[i]))
}
record('Independent running-sum enrichment scores',max(errs)<1e-10,paste('Maximum absolute ES difference',format(max(errs),digits=3)))
idx<-which(ora$status=='tested')
hp<-with(ora[idx,],phyper(hit_in_location-1,background_in_location,total_background-background_in_location,total_hits,lower.tail=FALSE))
record('ORA equals hypergeometric upper-tail test',max(abs(hp-ora$pval[idx]))<1e-10,'Fisher greater-tail p-values independently checked with phyper')
bh_ok<-function(df,cols) all(vapply(split(df,interaction(df[,cols],drop=TRUE)),function(z)all(abs(p.adjust(z$pval,'BH')-z$padj)<1e-12,na.rm=TRUE),logical(1)))
record('BH correction families',bh_ok(g,'analysis')&&bh_ok(ora,c('analysis','isoform')),'Separate GSEA collections; separate ORA collections and isoform sets')
record('ORA background is measured HPA-mapped genes',all(ora$total_background==1641)&&all(x$ranked$mapped_to_hpa),'No unmeasured HPA genes in Fisher background')
record('No significant results',sum(g$padj<.05,na.rm=TRUE)==0&&sum(ora$padj<.05,na.rm=TRUE)==0,'All four GSEA and all eight ORA families checked')
pdfs<-list.files(file.path(o,'figures'),pattern='\\.pdf$',full.names=TRUE)
pngs<-list.files(file.path(o,'figures'),pattern='\\.png$',full.names=TRUE)
record('Paired PDF and PNG figures',length(pdfs)==14&&length(pngs)==14&&setequal(sub('\\.pdf$','',basename(pdfs)),sub('\\.png$','',basename(pngs))),'14 figure pairs')
write_tsv(bind_rows(checks),file.path(o,'validation_checks.tsv'))
print(bind_rows(checks))
