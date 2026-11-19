suppressPackageStartupMessages(library(data.table))
out <- 'output/followup_analysis/intron_position_matched_clustering'
checks <- list(); expected_rows <- list()
for(model in c('position_matched','first20_uniform')) {
  d <- readRDS(file.path(out,'cache',paste0(model,'_distances.rds')))
  z <- d$site_map; groups <- split(seq_len(nrow(z)),z$feature_id)
  stopifnot(all(lengths(groups)>=2),all(z$analysis_sites==rep(lengths(groups),lengths(groups))),
    ncol(d$randomized_positions)==500,all(d$randomized_positions>=12.5),
    all(d$randomized_positions<=z$length_nt-12.5),all(d$randomized_positions%%1==.5))
  if(model=='first20_uniform') stopifnot(all(d$randomized_positions/z$length_nt<.2),
                                         all(d$observed_positions/z$length_nt<.2))
  for(it in c(0L,1L,250L,500L)) {
    pos <- if(it==0) d$observed_positions else d$randomized_positions[,it]
    nn <- if(it==0) d$observed_nearest else d$nearest_neighbor_null[,it]
    for(ix in groups) {
      m <- abs(outer(pos[ix],pos[ix],'-'));diag(m)<-Inf
      stopifnot(identical(as.numeric(apply(m,1,min)),as.numeric(nn[ix])))
    }
  }
  # Analytic binned expectation accounts for discretization and margin rejection.
  pp <- d$empirical_bin_probabilities; cumulative <- c(0,cumsum(pp))
  cdf <- function(u) {
    t <- pmax(0,pmin(1,u))*100; b <- pmin(99L,floor(t))
    cumulative[b+1L]+(t-b)*pp[b+1L]
  }
  lengths <- z[,.(n=.N),by=length_nt]
  expected_counts <- variance_counts <- numeric(100)
  for(j in seq_len(nrow(lengths))) {
    L <- lengths$length_nt[j]; n <- lengths$n[j]
    lo <- pmax(12,ceiling((0:99)*L/100-.5))
    hi <- pmin(L-13,ceiling((1:100)*L/100-.5)-1)
    if(model=='position_matched') {
      probs <- pmax(0,cdf((hi+1)/L)-cdf(lo/L))/(cdf((L-12)/L)-cdf(12/L))
    } else {
      max_k <- ceiling(.2*L-.5)-1
      hi <- pmin(hi,max_k)
      probs <- pmax(0,hi-lo+1)/(max_k-12+1)
    }
    stopifnot(abs(sum(probs)-1)<1e-10)
    expected_counts <- expected_counts+n*probs
    variance_counts <- variance_counts+n*probs*(1-probs)
  }
  expected <- expected_counts/nrow(z)
  mc_se <- sqrt(variance_counts/500)/nrow(z)
  actual <- rowMeans(d$relative_position_null)
  stopifnot(all(abs(actual-expected)<=6*mc_se+1e-12))
  expected_rows[[model]] <- data.table(model=model,bin=0:99,position_percent=(0:99)+.5,
    empirical_global=pp,analytic_expected_probability=expected,
    simulated_mean_probability=actual,monte_carlo_se=mc_se,
    standardized_error=ifelse(mc_se>0,(actual-expected)/mc_se,0))
  checks[[model]] <- data.table(model=model,features=length(groups),summits=nrow(z),
    max_position_bin_MC_z=max(abs((actual[mc_se>0]-expected[mc_se>0])/mc_se[mc_se>0])),
    checks='Pairwise NN on observed/1/250/500; all random-position bounds; exact feature counts; analytic positional expectation within 6 MC SE',result='PASS')
}
fwrite(rbindlist(expected_rows),file.path(out,'position_matching_validation.tsv'),sep='\t')
fwrite(rbindlist(checks),file.path(out,'validation.tsv'),sep='\t')
original <- readRDS('output/followup_analysis/site_clustering/cache/intron_site_distances.rds')
new <- readRDS(file.path(out,'cache','position_matched_distances.rds'))
stopifnot(identical(original$site_map$feature_id,new$site_map$feature_id),
  identical(original$site_map$peak_id,new$site_map$peak_id),
  identical(original$observed_positions,new$observed_positions),
  identical(original$observed_nearest,new$observed_nearest))
manifest <- fread(file.path(out,'input_provenance.tsv'))
stopifnot(identical(unname(tools::md5sum(manifest$path)),manifest$md5))
print(rbindlist(checks))
cat('PASS: original full-intron cohort, assignments, observed distances and source files unchanged.\n')
