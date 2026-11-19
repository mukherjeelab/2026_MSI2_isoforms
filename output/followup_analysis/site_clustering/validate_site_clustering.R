suppressPackageStartupMessages(library(data.table))
out <- 'output/followup_analysis/site_clustering'
checks <- list()
for (kind_name in c('intron','3UTR')) {
  d <- readRDS(file.path(out,'cache',paste0(kind_name,'_site_distances.rds')))
  groups <- split(seq_len(nrow(d$site_map)),d$site_map$feature_id)
  stopifnot(all(lengths(groups)>=2),all(d$site_map$n_sites==rep(lengths(groups),lengths(groups))),
            ncol(d$randomized_positions)==500,
            all(d$randomized_positions>=12.5),
            all(d$randomized_positions<=d$site_map$length_nt-12.5),
            all(d$randomized_positions%%1==.5))
  # Independent pairwise-distance implementation, not the adjacent-gap algorithm.
  for (it in c(0L,1L,250L,500L)) {
    pos <- if(it==0L) d$observed_positions else d$randomized_positions[,it]
    expected <- if(it==0L) d$observed_nearest else d$nearest_neighbor_null[,it]
    expected_gaps <- if(it==0L) d$observed_gaps else d$consecutive_gap_null[,it]
    recomputed_gaps <- unlist(lapply(groups,function(ix) {
      pairwise <- abs(outer(pos[ix],pos[ix],'-')); diag(pairwise) <- Inf
      stopifnot(identical(as.numeric(apply(pairwise,1,min)),as.numeric(expected[ix])))
      diff(sort(pos[ix]))
    }),use.names=FALSE)
    stopifnot(identical(recomputed_gaps,as.numeric(expected_gaps)))
  }
  checks[[kind_name]] <- data.table(kind=kind_name,features=length(groups),
    summits=nrow(d$site_map),iterations_checked='observed;1;250;500',
    check='Independent pairwise NN and sorted-gap reconstruction; all 500 position bounds and exact feature counts',
    result='PASS')
}
source <- fread(file.path(out,'input_provenance.tsv'))
stopifnot(identical(unname(tools::md5sum(source$path)),source$md5))
fwrite(rbindlist(checks),file.path(out,'validation.tsv'),sep='\t')
cat('PASS: independent distance reconstruction, null bounds, feature counts, and source hashes.\n')
