#!/usr/bin/env Rscript
# Reconstruct the notebook inputs by executing its original definitions verbatim.
suppressPackageStartupMessages({library(readxl);library(readr);library(dplyr);library(tidyr);library(tibble);library(ggplot2);library(pheatmap);library(clusterProfiler);library(org.Hs.eg.db)})
out <- 'output/proteomics_hpa/reactome_input_audit'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
src <- readLines('code/Proetomics.Rmd')
fragments <- list()
chunk <- function(name) {
  start <- which(startsWith(src,paste0('```{r ',name)))[1]
  stopifnot(!is.na(start))
  end <- which(seq_along(src)>start & src=='```')[1]
  src[seq.int(start+1,end-1)]
}
run <- function(lines,name) {
  fragments[[name]] <<- lines
  for(expr in parse(text=lines)) eval(expr,envir=.GlobalEnv)
}
data_raw <- read_xlsx('data/FredHutchProteomics.xlsx')
run(chunk('clean-data'),'clean-data')
run(chunk('heatmap-functions'),'heatmap-functions')
run(chunk('heatmaps-and-clusters'),'heatmaps-and-clusters')
start <- which(startsWith(src,'filter_cluster_hits <- function'))[1]
end <- which(src=='length(genes_324v328_hits)')[1]
lines <- src[start:end]
run(lines[!startsWith(lines,'```')],'original-hit-definitions')
start <- which(startsWith(src,'map_symbols_to_entrez <- function'))[1]
end <- which(startsWith(src,'make_ora <- function'))[1]-1
run(src[start:end],'original-entrez-mapping-function')
universe_genes <- prot_df %>% pull(gene) %>% unique()
sets <- list(genes_324_cluster=genes_324_cluster,genes_328_cluster=genes_328_cluster,
             genes_324_specific=genes_324_specific,genes_328_specific=genes_328_specific,
             genes_324_hits=genes_324_hits,genes_328_hits=genes_328_hits)
mapped <- map_symbols_to_entrez(universe_genes)
universe_audit <- left_join(tibble(SYMBOL=universe_genes),mapped,by='SYMBOL')
members <- bind_rows(lapply(names(sets),function(n)tibble(set_name=n,SYMBOL=sets[[n]])))
mapped_members <- left_join(members,mapped,by='SYMBOL',relationship='many-to-many')
write_tsv(universe_audit,file.path(out,'original_universe_symbols_entrez.tsv'))
write_tsv(mapped_members,file.path(out,'original_set_memberships.tsv'))
summary <- mapped_members %>% group_by(set_name) %>% summarise(original_symbols=n_distinct(SYMBOL),
  mapped_entrez=n_distinct(ENTREZID,na.rm=TRUE),unmapped_symbols=n_distinct(SYMBOL[is.na(ENTREZID)]),.groups='drop')
write_tsv(summary,file.path(out,'input_set_summary.tsv'))
saveRDS(list(sets=sets,universe_genes=universe_genes,entrez_mapping=mapped),file.path(out,'original_reactome_inputs.rds'))
writeLines(unlist(lapply(names(fragments),function(n)c(paste0('# Notebook fragment: ',n),fragments[[n]],''))),file.path(out,'executed_notebook_definitions.R'))
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
stopifnot(all(unlist(sets) %in% universe_genes))
print(summary)
cat('Original universe:',length(universe_genes),'symbols;',n_distinct(mapped$ENTREZID),'Entrez IDs\n')
