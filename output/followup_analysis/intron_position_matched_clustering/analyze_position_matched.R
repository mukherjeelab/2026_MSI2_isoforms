suppressPackageStartupMessages(library(data.table))
out <- 'output/followup_analysis/intron_position_matched_clustering'
dir.create(file.path(out,'cache'),recursive=TRUE,showWarnings=FALSE)
position_path <- 'output/followup_analysis/tables/summit_normalized_positions.tsv.gz'
included_path <- 'output/followup_analysis/site_clustering/MSI2_site_clustering_included_summits.tsv.gz'
uniform_path <- 'output/followup_analysis/site_clustering/MSI2_site_clustering_all_statistics.tsv'
paths <- c(position_path,included_path,uniform_path)
hashes <- tools::md5sum(paths)
B <- 500L
pool <- fread(position_path)[kind=='intron' & length>500]
pool[,position_nt:=offset+.5]
stopifnot(all(abs(pool$fraction-pool$position_nt/pool$length)<1e-12))
pool[,bin:=pmin(99L,floor(100*fraction))]
bin_counts <- tabulate(pool$bin+1L,nbins=100)
prob <- bin_counts/sum(bin_counts)
fit <- data.table(bin=0:99,position_percent=(0:99)+.5,summits=bin_counts,probability=prob)
fwrite(fit,file.path(out,'empirical_intron_position_distribution.tsv'),sep='\t')
fwrite(pool,file.path(out,'empirical_position_source_summits.tsv.gz'),sep='\t')
zfull <- fread(included_path)[kind=='intron']
setorder(zfull,feature_id,observed_position_nt,peak_id)
stopifnot(nrow(zfull)==14598,uniqueN(zfull$feature_id)==3799,
          all(zfull$length_nt>500),all(zfull$n_sites>=2),
          all(zfull$peak_id %in% pool$peak_id))
zfull[,total_intron_sites:=n_sites]
zfull[,analysis_sites:=.N,by=feature_id]
zone <- zfull[observed_position_nt/length_nt<.2]
zone[,analysis_sites:=.N,by=feature_id]
zone_counts <- unique(zone[,.(feature_id,length_nt,total_intron_sites,analysis_sites)])
fwrite(zone_counts,file.path(out,'first20_zone_feature_counts_before_minimum.tsv'),sep='\t')
zone <- zone[analysis_sites>=2]
zone[,stratum:=fifelse(analysis_sites==2,'2 sites',fifelse(analysis_sites<=4,'3-4 sites','>=5 sites'))]
fwrite(zfull,file.path(out,'full_intron_analysis_summits.tsv.gz'),sep='\t')
fwrite(zone,file.path(out,'first20_analysis_summits.tsv.gz'),sep='\t')

nn_distances <- function(pos,group) {
  o <- order(group,pos); delta <- diff(pos[o]); delta[diff(group[o])!=0] <- Inf
  ans <- numeric(length(pos)); ans[o] <- pmin(c(Inf,delta),c(delta,Inf)); ans
}
metrics <- function(x) c(median_distance_nt=median(x),within_50_nt=mean(x<=50),
  within_100_nt=mean(x<=100),within_250_nt=mean(x<=250),within_500_nt=mean(x<=500))

# Draw from the empirical piecewise-constant relative-position density, then
# map into the recipient intron's base-centered half-integer coordinate grid.
# Rejection conditions this density on the original 12.5 nt end margins.
draw_position_matched <- function(lengths) {
  result <- numeric(length(lengths)); pending <- seq_along(lengths)
  while(length(pending)) {
    bin <- sample.int(100,length(pending),replace=TRUE,prob=prob)-1L
    u <- (bin+runif(length(pending)))/100
    p <- floor(u*lengths[pending])+.5
    ok <- p>=12.5 & p<=lengths[pending]-12.5
    result[pending[ok]] <- p[ok]
    pending <- pending[!ok]
  }
  result
}
observed_rows <- list(); iteration_rows <- list(); profile_rows <- list(); diagnostics <- list()
models <- c('position_matched','first20_uniform')
seeds <- c(position_matched=26092301L,first20_uniform=26092302L)
for(model_name in models) {
  z <- if(model_name=='position_matched') copy(zfull) else copy(zone)
  setorder(z,feature_id,observed_position_nt,peak_id)
  group <- match(z$feature_id,unique(z$feature_id)); n <- nrow(z)
  observed <- nn_distances(z$observed_position_nt,group)
  strata <- c('All','2 sites','3-4 sites','>=5 sites')
  masks <- lapply(strata,function(s) if(s=='All') rep(TRUE,n) else z$stratum==s)
  counts <- rbindlist(lapply(seq_along(strata),function(j) data.table(stratum=strata[j],
    n_features=uniqueN(z$feature_id[masks[[j]]]),n_summits=sum(masks[[j]]))))
  stopifnot(all(counts$n_features>0),all(is.finite(observed)))
  observed_rows[[model_name]] <- merge(rbindlist(lapply(seq_along(strata),function(j) {
    v <- metrics(observed[masks[[j]]])
    data.table(model=model_name,stratum=strata[j],metric=names(v),observed=unname(v))
  })),counts,by='stratum')
  posmat <- matrix(NA_real_,n,B); nnmat <- posmat; profiles <- matrix(0,100,B)
  rows <- vector('list',B); set.seed(seeds[[model_name]])
  for(it in 1:B) {
    if(model_name=='position_matched') {
      pos <- draw_position_matched(z$length_nt)
    } else {
      # Allowed half-integer sites start at 12.5 and end strictly before .2*L.
      choices <- ceiling(.2*z$length_nt-12.5)
      pos <- 12.5+floor(runif(n)*choices)
      stopifnot(all(pos/z$length_nt<.2))
    }
    stopifnot(all(pos>=12.5),all(pos<=z$length_nt-12.5))
    dd <- nn_distances(pos,group)
    posmat[,it] <- pos; nnmat[,it] <- dd
    profiles[,it] <- tabulate(pmin(99L,floor(100*pos/z$length_nt))+1L,100)/n
    rows[[it]] <- rbindlist(lapply(seq_along(strata),function(j) {
      v <- metrics(dd[masks[[j]]])
      data.table(model=model_name,stratum=strata[j],iteration=it,metric=names(v),value=unname(v))
    }))
    if(it %% 100==0) cat(model_name,it,'of',B,'complete\n')
  }
  iteration_rows[[model_name]] <- rbindlist(rows)
  profile_rows[[model_name]] <- data.table(model=model_name,bin=0:99,
    position_percent=(0:99)+.5,global_empirical_probability=prob,
    analysis_observed=tabulate(pmin(99L,floor(100*z$observed_position_nt/z$length_nt))+1L,100)/n,
    null_mean=rowMeans(profiles),null_low=apply(profiles,1,quantile,.025),
    null_high=apply(profiles,1,quantile,.975))
  diagnostics[[model_name]] <- data.table(model=model_name,n_features=uniqueN(group),n_summits=n,
    outside_null_margin=sum(z$observed_position_nt<12.5 | z$observed_position_nt>z$length_nt-12.5),
    observed_first20_fraction=mean(z$observed_position_nt/z$length_nt<.2),
    null_mean_first20_fraction=mean(colSums(profiles[1:20,,drop=FALSE])))
  saveRDS(list(model=model_name,seed=seeds[[model_name]],iterations=B,
    site_map=z[,.(feature_id,peak_id,strand,length_nt,total_intron_sites,analysis_sites,stratum)],
    observed_positions=z$observed_position_nt,observed_nearest=observed,
    randomized_positions=posmat,nearest_neighbor_null=nnmat,relative_position_null=profiles,
    empirical_bin_probabilities=prob,final_rng_state=.Random.seed),
    file.path(out,'cache',paste0(model_name,'_distances.rds')))
}
iterations <- rbindlist(iteration_rows); observed <- rbindlist(observed_rows)
fwrite(iterations,file.path(out,'intron_position_matched_null_iterations.tsv.gz'),sep='\t')
fwrite(rbindlist(profile_rows),file.path(out,'positional_null_diagnostics.tsv'),sep='\t')
fwrite(rbindlist(diagnostics),file.path(out,'cohort_diagnostics.tsv'),sep='\t')
joined <- merge(iterations,observed,by=c('model','stratum','metric'),allow.cartesian=TRUE)
stats <- joined[,.(observed=observed[1],null_mean=mean(value),null_median=median(value),
  null_low_95=quantile(value,.025),null_high_95=quantile(value,.975),
  p_lower=(1+sum(value<=observed[1]))/(B+1),p_upper=(1+sum(value>=observed[1]))/(B+1),
  n_features=n_features[1],n_summits=n_summits[1]),by=.(model,stratum,metric)]
stats[,null_reference:=fifelse(metric=='median_distance_nt',null_median,null_mean)]
stats[,`:=`(observed_null_ratio=observed/null_reference,
  p_clustering=fifelse(metric=='median_distance_nt',p_lower,p_upper),
  p_two_sided=pmin(1,2*pmin(p_lower,p_upper)),iterations=B,seed=unname(seeds[model]))]
stats[,q_BH_clustering:=p.adjust(p_clustering,'BH')]
fwrite(stats,file.path(out,'intron_position_matched_all_statistics.tsv'),sep='\t')
fwrite(stats[stratum=='All'],file.path(out,'intron_position_matched_summary.tsv'),sep='\t')
fwrite(stats[stratum!='All'],file.path(out,'intron_position_matched_by_sitecount.tsv'),sep='\t')
old <- fread(uniform_path)[kind=='intron' & distance_type=='nearest_neighbor']
old[,model:='original_uniform']
common <- intersect(names(stats),names(old))
comparison <- rbindlist(list(old[,..common],stats[,..common]),use.names=TRUE)
fwrite(comparison,file.path(out,'uniform_vs_position_matched_vs_first20.tsv'),sep='\t')
stopifnot(identical(hashes,tools::md5sum(paths)))
fwrite(data.table(path=paths,md5=unname(hashes)),file.path(out,'input_provenance.tsv'),sep='\t')
writeLines(capture.output(sessionInfo()),file.path(out,'session_info.txt'))
print(comparison[stratum=='All' & metric %in% c('median_distance_nt','within_250_nt'),
  .(model,metric,n_features,n_summits,observed,null_reference,null_low_95,null_high_95,
    observed_null_ratio,p_clustering)])
