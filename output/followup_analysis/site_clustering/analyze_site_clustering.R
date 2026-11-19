suppressPackageStartupMessages(library(data.table))
out <- 'output/followup_analysis/site_clustering'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out,'cache'),showWarnings=FALSE)
input <- 'output/followup_analysis/tables/summit_normalized_positions.tsv.gz'
feature_input <- 'data/results_round2/features.tsv.gz'
source_script <- 'output/followup_analysis/scripts/site_spatial.R'
paths <- c(input,feature_input,source_script)
hashes <- tools::md5sum(paths)
seed <- 26092201L; B <- 500L; set.seed(seed)
hits <- fread(input)[kind %in% c('intron','3UTR')]
features <- fread(feature_input)
stopifnot(!anyDuplicated(features$feature_id),!anyDuplicated(hits[,.(kind,peak_id)]))
hits <- merge(hits,features[,.(feature_id,chrom,strand,length_nt,segments_0based_half_open,
                                gene_ids,gene_names,transcript_ids)],by='feature_id',sort=FALSE)
stopifnot(all(hits$length==hits$length_nt),all(hits$offset>=0),all(hits$offset<hits$length))
hits[,eligible := (kind=='intron' & length>500) | (kind=='3UTR' & length>=100)]
eligibility <- hits[,.(assigned_features=uniqueN(feature_id),assigned_summits=.N,
  length_eligible_features=uniqueN(feature_id[eligible]),length_eligible_summits=sum(eligible)),by=kind]
hits <- hits[eligible==TRUE]
hits[,n_sites:=.N,by=.(kind,feature_id)]
stopifnot(all(hits$n_sites==hits$sites_in_feature))
included <- hits[n_sites>=2]
included[,stratum:=fifelse(n_sites==2,'2 sites',fifelse(n_sites<=4,'3-4 sites','>=5 sites'))]
setorder(included,kind,feature_id,offset,peak_id)
included[,observed_position_nt:=offset+.5]
fwrite(included,file.path(out,'MSI2_site_clustering_included_summits.tsv.gz'),sep='\t')
feature_table <- unique(included[,.(kind,feature_id,chrom,strand,length_nt,n_sites,stratum,
                          segments_0based_half_open,gene_ids,gene_names,transcript_ids)])
fwrite(feature_table,file.path(out,'MSI2_site_clustering_features.tsv.gz'),sep='\t')
eligibility <- merge(eligibility,included[,.(multi_site_features=uniqueN(feature_id),
                       contributing_summits=.N),by=kind],by='kind')
fwrite(eligibility,file.path(out,'MSI2_site_clustering_eligibility.tsv'),sep='\t')

# Positions are transcript-oriented; nearest distances are returned in input row order.
distances <- function(pos,group) {
  ord <- order(group,pos)
  delta <- diff(pos[ord]); same <- diff(group[ord])==0
  delta[!same] <- Inf
  nearest_sorted <- pmin(c(Inf,delta),c(delta,Inf))
  nearest <- numeric(length(pos)); nearest[ord] <- nearest_sorted
  list(nearest=nearest,gaps=delta[same],left=ord[-length(ord)][same],right=ord[-1][same])
}
# Small nonuniform example confirms NN and adjacent gaps are not interchangeable.
check <- distances(c(0,10,100,1000,1025),c(1,1,1,2,2))
stopifnot(identical(check$nearest,c(10,10,90,25,25)),identical(check$gaps,c(10,90,25)))
metrics <- function(x) c(median_distance_nt=median(x),within_50_nt=mean(x<=50),
  within_100_nt=mean(x<=100),within_250_nt=mean(x<=250),within_500_nt=mean(x<=500))
all_iterations <- list(); all_observed <- list(); all_distances <- list(); all_gaps <- list()
for (kind_name in c('intron','3UTR')) {
  z <- included[kind==kind_name]
  group <- match(z$feature_id,unique(z$feature_id)); n <- nrow(z)
  obs <- distances(z$observed_position_nt,group)
  gap_stratum <- z$stratum[obs$left]
  strata <- c('All','2 sites','3-4 sites','>=5 sites')
  masks <- lapply(strata,function(s) if(s=='All') rep(TRUE,n) else z$stratum==s)
  gap_masks <- lapply(strata,function(s) if(s=='All') rep(TRUE,length(obs$gaps)) else gap_stratum==s)
  count_table <- rbindlist(lapply(seq_along(strata),function(j) data.table(
    stratum=strata[j],n_features=uniqueN(z$feature_id[masks[[j]]]),n_summits=sum(masks[[j]]),
    n_consecutive_gaps=sum(gap_masks[[j]]))))
  stopifnot(all(count_table$n_features>0),all(is.finite(obs$nearest)),
            length(obs$gaps)==n-uniqueN(group))
  observed_rows <- rbindlist(lapply(seq_along(strata),function(j) rbindlist(lapply(
    c('nearest_neighbor','consecutive_gap'),function(metric_type) {
      v <- if(metric_type=='nearest_neighbor') metrics(obs$nearest[masks[[j]]]) else metrics(obs$gaps[gap_masks[[j]]])
      data.table(kind=kind_name,stratum=strata[j],distance_type=metric_type,metric=names(v),observed=unname(v))
    }))))
  all_observed[[kind_name]] <- merge(observed_rows,count_table,by='stratum',sort=FALSE)
  all_distances[[kind_name]] <- cbind(z[,.(kind,feature_id,peak_id,strand,length_nt,n_sites,stratum,
                         observed_position_nt)],nearest_neighbor_nt=obs$nearest)
  all_gaps[[kind_name]] <- data.table(kind=kind_name,feature_id=z$feature_id[obs$left],
    stratum=gap_stratum,left_peak_id=z$peak_id[obs$left],right_peak_id=z$peak_id[obs$right],
    left_position_nt=z$observed_position_nt[obs$left],right_position_nt=z$observed_position_nt[obs$right],
    consecutive_gap_nt=obs$gaps)
  random_positions <- matrix(NA_real_,n,B)
  nearest_null <- matrix(NA_real_,n,B)
  gaps_null <- matrix(NA_real_,length(obs$gaps),B)
  iteration_rows <- vector('list',B)
  for (it in seq_len(B)) {
    # Exact same independent-grid distribution as the existing spatial adaptation.
    pos <- 12.5 + floor(runif(n)*(floor(z$length-25)+1))
    d <- distances(pos,group)
    stopifnot(all(pos>=12.5),all(pos<=z$length-12.5),all(is.finite(d$nearest)))
    random_positions[,it] <- pos; nearest_null[,it] <- d$nearest; gaps_null[,it] <- d$gaps
    iteration_rows[[it]] <- rbindlist(lapply(seq_along(strata),function(j) {
      v1 <- metrics(d$nearest[masks[[j]]]); v2 <- metrics(d$gaps[gap_masks[[j]]])
      data.table(kind=kind_name,stratum=strata[j],iteration=it,
        distance_type=rep(c('nearest_neighbor','consecutive_gap'),each=length(v1)),
        metric=rep(names(v1),2),value=c(unname(v1),unname(v2)))
    }))
    if(it %% 100==0) cat(kind_name,it,'of',B,'randomizations complete\n')
  }
  all_iterations[[kind_name]] <- rbindlist(iteration_rows)
  saveRDS(list(kind=kind_name,seed=seed,iterations=B,site_map=z[,.(feature_id,peak_id,strand,length_nt,n_sites,stratum)],
    gap_feature_id=z$feature_id[obs$left],gap_stratum=gap_stratum,
    observed_positions=z$observed_position_nt,observed_nearest=obs$nearest,observed_gaps=obs$gaps,
    randomized_positions=random_positions,nearest_neighbor_null=nearest_null,
    consecutive_gap_null=gaps_null),file.path(out,'cache',paste0(kind_name,'_site_distances.rds')))
}
iterations <- rbindlist(all_iterations)
observed <- rbindlist(all_observed)
fwrite(iterations,file.path(out,'MSI2_site_clustering_null_iterations.tsv.gz'),sep='\t')
fwrite(rbindlist(all_distances),file.path(out,'MSI2_observed_nearest_neighbor_distances.tsv.gz'),sep='\t')
fwrite(rbindlist(all_gaps),file.path(out,'MSI2_observed_consecutive_distances.tsv.gz'),sep='\t')
keys <- c('kind','stratum','distance_type','metric')
joined <- merge(iterations,observed,by=keys,allow.cartesian=TRUE)
stats <- joined[,.(observed=observed[1],null_mean=mean(value),null_median=median(value),
  null_low_95=quantile(value,.025),null_high_95=quantile(value,.975),
  p_lower=(1+sum(value<=observed[1]))/(B+1),p_upper=(1+sum(value>=observed[1]))/(B+1),
  n_features=n_features[1],n_summits=n_summits[1],n_consecutive_gaps=n_consecutive_gaps[1]),by=keys]
stats[,null_reference:=fifelse(metric=='median_distance_nt',null_median,null_mean)]
stats[,`:=`(observed_null_ratio=observed/null_reference,
  p_clustering=fifelse(metric=='median_distance_nt',p_lower,p_upper),
  p_two_sided=pmin(1,2*pmin(p_lower,p_upper)),iterations=B,seed=seed)]
stats[,q_BH_clustering:=p.adjust(p_clustering,'BH'),by=distance_type]
fwrite(stats[distance_type=='nearest_neighbor' & stratum=='All'],
       file.path(out,'MSI2_site_clustering_summary.tsv'),sep='\t')
fwrite(stats[distance_type=='nearest_neighbor' & stratum!='All'],
       file.path(out,'MSI2_site_clustering_by_sitecount.tsv'),sep='\t')
fwrite(stats[distance_type=='nearest_neighbor' & metric!='median_distance_nt'],
       file.path(out,'MSI2_site_clustering_thresholds.tsv'),sep='\t')
fwrite(stats[distance_type=='consecutive_gap'],file.path(out,'MSI2_consecutive_site_distance_summary.tsv'),sep='\t')
fwrite(stats,file.path(out,'MSI2_site_clustering_all_statistics.tsv'),sep='\t')
margin <- included[,.(summits=.N,observed_sites_outside_null_margin=sum(
  observed_position_nt<12.5 | observed_position_nt>length-12.5)),by=kind]
fwrite(margin,file.path(out,'MSI2_site_clustering_margin_diagnostics.tsv'),sep='\t')
stopifnot(identical(hashes,tools::md5sum(paths)))
fwrite(data.table(path=paths,md5=unname(hashes)),file.path(out,'input_provenance.tsv'),sep='\t')
saveRDS(.Random.seed,file.path(out,'cache','final_R_rng_state.rds'))
writeLines(c(capture.output(sessionInfo()),paste('seed',seed)),file.path(out,'session_info.txt'))
print(stats[distance_type=='nearest_neighbor' & stratum=='All',
  .(kind,metric,n_features,n_summits,observed,null_reference,null_low_95,null_high_95,
    observed_null_ratio,p_clustering)])
