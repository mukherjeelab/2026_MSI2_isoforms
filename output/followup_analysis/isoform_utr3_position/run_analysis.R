suppressPackageStartupMessages({library(data.table); library(ggplot2); library(GenomicRanges); library(GenomicAlignments); library(jsonlite)})
out <- 'output/followup_analysis/isoform_utr3_position'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
setDTthreads(2)
seed <- 260922L
savefig <- function(p, name, w=10, h=6.8) {
  for (ext in c('png','pdf')) ggsave(file.path(out,paste0(name,'.',ext)),p,width=w,height=h,dpi=240,bg='white')
}
sty <- theme_classic(base_size=12) + theme(legend.position='bottom',legend.box='vertical',
  plot.caption=element_text(hjust=0,size=9),plot.title=element_text(size=16),
  strip.background=element_rect(fill='#F2F2F2',color='grey70'),plot.margin=margin(12,16,12,12))
profile <- function(bin,w) {
  v <- rowsum(w,bin+1L,reorder=FALSE); ans <- numeric(100)
  ans[as.integer(rownames(v))] <- v[,1]
  ans/sum(w)
}
refsites <- fread('output/followup_analysis/tables/summit_normalized_positions.tsv.gz')[kind=='3UTR' & length>=100]
refids <- unique(refsites$feature_id)
refsites[,n_feature:=.N,by=feature_id]
rw <- 1/refsites$n_feature
set.seed(seed)
rn <- matrix(0,100,500)
for (it in 1:500) {
  rp <- 12.5+floor(runif(nrow(refsites))*(refsites$length-25+1))
  rn[,it] <- profile(pmin(99L,floor(100*rp/refsites$length)),rw)
}
ref <- data.table(position=(0:99)+.5,observed=profile(floor(100*refsites$fraction),rw),
  background_median=apply(rn,1,median),null_p10=apply(rn,1,quantile,.1),null_p25=apply(rn,1,quantile,.25),
  null_p75=apply(rn,1,quantile,.75),null_p90=apply(rn,1,quantile,.9))
stopifnot(abs(sum(ref$observed)-1)<1e-12,max(abs(colSums(rn)-1))<1e-12)
saveRDS(list(seed=seed,iterations=500,null=rn),file.path(out,'endogenous_equal_feature_nulls.rds'))
fwrite(refsites,file.path(out,'endogenous_assigned_summits.tsv.gz'),sep='\t')
cloud <- data.table(position=rep(ref$position,500),draw=rep(1:500,each=100),value=as.vector(rn))
null_layers <- function(p, cl) list(
  geom_ribbon(data=p,aes(position,ymin=null_p10,ymax=null_p90,fill='10th-90th percentiles')),
  geom_ribbon(data=p,aes(position,ymin=null_p25,ymax=null_p75,fill='25th-75th percentiles')),
  geom_line(data=cl,aes(position,value,group=draw),color='#46B8D5',alpha=.035,linewidth=.18),
  geom_line(data=p,aes(position,background_median),color='black',linewidth=.5),
  scale_fill_manual(name='Null bands',values=c('10th-90th percentiles'='#D6EBF5','25th-75th percentiles'='#7CAACD')))
reference_plot <- ggplot() + null_layers(ref,cloud) +
  geom_line(data=ref,aes(position,observed),color='#D62728',linewidth=.85) + sty +
  labs(title='Endogenous MSI2: equal 3UTR weight',
    subtitle=sprintf('%s summits in %s 3UTR features >=100 nt',format(nrow(refsites),big.mark=','),format(length(refids),big.mark=',')),x='3UTR position (%)',y='Weighted fraction per 1% bin',
    caption='Red: observed summits. Cyan: 500 randomized profiles. Black: null median.\nEach represented 3UTR contributes total weight one; original assignments retained; 500 new equal-feature null draws.')
savefig(reference_plot,'endogenous_utr3_equal_weight',9,6.6)
fwrite(ref,file.path(out,'endogenous_plot_inputs.tsv'),sep='\t')

features <- fread('data/results_round2/features.tsv.gz')[kind=='3UTR' & length_nt>=100]
segs <- lapply(features$segments_0based_half_open,fromJSON)
segments <- rbindlist(lapply(seq_along(segs),function(i) {
  z <- segs[[i]]
  stopifnot(sum(z[,2]-z[,1])==features$length_nt[i],all(z[,2]>z[,1]),!is.unsorted(z[,1]))
  data.table(fi=i,chrom=features$chrom[i],strand=features$strand[i],start0=z[,1],end0=z[,2],
    preceding=c(0,head(cumsum(z[,2]-z[,1]),-1)))
}))
igr <- GRanges(segments$chrom,IRanges(segments$start0+1,segments$end0),segments$strand)
# Validate segmented-coordinate reconstruction against the saved reference assignments.
anno <- fread('data/results/peak_analysis/figures/all_nonzero_peaks/transcript_region_bias/peak_transcript_region_annotations.tsv')
ra <- anno[match(refsites$peak_id,peak_id)]
rpg <- GRanges(ra$summit_chromosome,IRanges(ra$summit_start_0based+1,width=1),ra$summit_strand)
rh <- suppressWarnings(findOverlaps(rpg,igr))
check <- data.table(ri=queryHits(rh),si=subjectHits(rh))
check[,fi:=segments$fi[si]]; check[,length_nt:=features$length_nt[fi]]
setorder(check,ri,-length_nt,fi); check <- check[!duplicated(ri)]
check[,offset:=segments$preceding[si]+ra$summit_start_0based[ri]-segments$start0[si]]
check[features$strand[fi]=='-',offset:=length_nt-1-offset]
stopifnot(nrow(check)==nrow(refsites),all(features$feature_id[check$fi]==refsites$feature_id),
  max(abs((check$offset+.5)/check$length_nt-refsites$fraction))<1e-12)
fwrite(segments,file.path(out,'feature_segments.tsv.gz'),sep='\t')
meta <- fread('data/results_round2/sample_metadata.tsv')[assay=='CLIP'][order(construct,biological_replicate)]
assigned <- list(); qc <- list()
for (s in seq_len(nrow(meta))) {
  sid <- meta$sample_id[s]; message('Reconstructing aligned read centers: ',sid)
  d <- fread(file.path('data/results_round2',sid,'read_annotations.tsv.gz'))
  q <- fromJSON(file.path('data/results_round2',sid,'QC.json'))
  stopifnot(nrow(d)==q$census$all$genome,all(d$sample_id==sid),!anyDuplicated(d$read_name))
  blocks <- cigarRangesAlongReferenceSpace(d$cigar,pos=d$position_1based,ops=c('M','=','X'),reduce.ranges=TRUE)
  stopifnot(all(lengths(blocks)>0))
  br <- unlist(blocks)
  b <- data.table(ri=rep(seq_len(nrow(d)),lengths(blocks)),start1=start(br),width=width(br))
  b[,cum:=cumsum(width),by=ri]
  totals <- b[,.(aligned_length=sum(width)),by=ri]
  # Select the central aligned nucleotide; for even lengths take the 5'-central nucleotide.
  # Insertions, deletions, skipped reference intervals, and soft clips do not count as aligned bases.
  target <- ifelse(d$strand=='+',ceiling(totals$aligned_length/2),floor(totals$aligned_length/2)+1)
  centers <- b[cum>=target[ri] & cum-width<target[ri]]
  setorder(centers,ri)
  stopifnot(identical(centers$ri,seq_len(nrow(d))))
  centers[,center0:=start1+target[ri]-(cum-width)-2]
  pg <- GRanges(d$chrom,IRanges(centers$center0+1,width=1),d$strand)
  hh <- suppressWarnings(findOverlaps(pg,igr,ignore.strand=FALSE))
  a <- data.table(ri=queryHits(hh),si=subjectHits(hh))
  a[,fi:=segments$fi[si]]
  a[,length_nt:=features$length_nt[fi]]
  setorder(a,ri,-length_nt,fi)
  a <- a[!duplicated(ri)]
  a[,`:=`(feature_id=features$feature_id[fi],read_name=d$read_name[ri],chrom=d$chrom[ri],
    strand=d$strand[ri],center0=centers$center0[ri],aligned_length=totals$aligned_length[ri])]
  a[,offset0:=segments$preceding[si]+center0-segments$start0[si]]
  a[strand=='-',offset0:=length_nt-1-offset0]
  a[,`:=`(fraction=(offset0+.5)/length_nt,sample_id=sid,construct=as.character(meta$construct[s]),replicate=meta$biological_replicate[s])]
  stopifnot(all(a$fraction>0 & a$fraction<1),!anyDuplicated(a$ri))
  assigned[[sid]] <- a
  qc[[sid]] <- data.table(sample_id=sid,construct=meta$construct[s],replicate=meta$biological_replicate[s],
    eligible_alignments=nrow(d),assigned_read_centers=nrow(a),supported_features=uniqueN(a$feature_id),
    centers_within_12_5nt_of_end=sum(pmin(a$offset0+.5,a$length_nt-a$offset0-.5)<12.5))
}
reads <- rbindlist(assigned)
fwrite(reads,file.path(out,'assigned_read_centers.tsv.gz'),sep='\t')
fwrite(rbindlist(qc),file.path(out,'read_assignment_qc.tsv'),sep='\t')
fwrite(features,file.path(out,'candidate_utr3.tsv.gz'),sep='\t')
common <- Reduce(intersect,lapply(assigned,function(x)unique(x$feature_id)))
cohorts <- list(all_assigned=unique(reads$feature_id),reference_utr3=refids,common_four_libraries=common)
fwrite(rbindlist(lapply(names(cohorts),function(k)data.table(cohort=k,feature_id=cohorts[[k]]))),file.path(out,'cohort_memberships.tsv.gz'),sep='\t')
profiles <- list(); clouds <- list(); summaries <- list(); cache <- list()
set.seed(seed)
for (cohort in names(cohorts)) {
  message('500 randomizations per library: ',cohort)
  cache[[cohort]] <- list()
  for (sid in meta$sample_id) {
    a <- copy(reads[sample_id==sid & feature_id %in% cohorts[[cohort]]])
    stopifnot(nrow(a)>0)
    a[,n_feature:=.N,by=feature_id]; w <- 1/a$n_feature
    obs <- profile(floor(100*a$fraction),w)
    nu <- matrix(0,100,500)
    for(it in 1:500) {
      rp <- 12.5+floor(runif(nrow(a))*(a$length_nt-25+1))
      nu[,it] <- profile(pmin(99L,floor(100*rp/a$length_nt)),w)
    }
    stopifnot(abs(sum(obs)-1)<1e-10,max(abs(colSums(nu)-1))<1e-10)
    cache[[cohort]][[sid]] <- list(observed=obs,null=nu,reads=nrow(a),features=uniqueN(a$feature_id))
  }
  for (construct in c('324','328')) {
    ids <- meta$sample_id[as.character(meta$construct)==construct]
    # Equal weight for biological replicates, not pooled alignment counts.
    ob <- rowMeans(sapply(ids,function(s)cache[[cohort]][[s]]$observed))
    nm <- (cache[[cohort]][[ids[1]]]$null+cache[[cohort]][[ids[2]]]$null)/2
    pp <- data.table(cohort=cohort,construct=construct,position=(0:99)+.5,observed=ob,
      background_median=apply(nm,1,median),null_p10=apply(nm,1,quantile,.1),null_p25=apply(nm,1,quantile,.25),
      null_p75=apply(nm,1,quantile,.75),null_p90=apply(nm,1,quantile,.9))
    for (i in 1:2) {
      repvalues <- cache[[cohort]][[ids[i]]]$observed
      set(pp,j=paste0('replicate_',meta$biological_replicate[match(ids[i],meta$sample_id)]),value=repvalues)
    }
    profiles[[paste(cohort,construct)]] <- pp
    clouds[[paste(cohort,construct)]] <- data.table(cohort=cohort,construct=construct,position=rep(pp$position,500),draw=rep(1:500,each=100),value=as.vector(nm))
    for (label in c(ids,'replicate_mean')) {
      o <- if(label=='replicate_mean')ob else cache[[cohort]][[label]]$observed
      nn <- if(label=='replicate_mean')nm else cache[[cohort]][[label]]$null
      for (zone in c('first10','first20','last10')) {
        ix <- switch(zone,first10=1:10,first20=1:20,last10=91:100)
        nv <- colSums(nn[ix,,drop=FALSE]); ov <- sum(o[ix])
        summaries[[length(summaries)+1]] <- data.table(cohort=cohort,construct=construct,sample_id=label,zone=zone,
          observed_fraction=ov,null_mean=mean(nv),null_low=quantile(nv,.025),null_high=quantile(nv,.975),
          observed_null_ratio=ov/mean(nv),features=if(label=='replicate_mean')NA_integer_ else cache[[cohort]][[label]]$features,
          reads=if(label=='replicate_mean')NA_integer_ else cache[[cohort]][[label]]$reads)
      }
    }
  }
}
pr <- rbindlist(profiles); cl <- rbindlist(clouds); su <- rbindlist(summaries)
saveRDS(list(seed=seed,iterations=500,cohorts=cache),file.path(out,'read_position_nulls.rds'))
fwrite(pr,file.path(out,'read_position_profiles.tsv'),sep='\t')
fwrite(cl,file.path(out,'replicate_mean_null_curves.tsv.gz'),sep='\t')
fwrite(su,file.path(out,'positional_zone_summary.tsv'),sep='\t')
for (cohort in names(cohorts)) {
  sel <- pr$cohort==cohort; p <- pr[sel]
  sel_cloud <- cl$cohort==cohort; cc <- cl[sel_cloud]
  rr <- melt(p,id.vars=c('construct','position'),measure.vars=c('replicate_A','replicate_B'),variable.name='replicate',value.name='value')
  rr[,replicate:=sub('replicate_','Replicate ',replicate)]
  ttl <- switch(cohort,all_assigned='All assigned 3UTR features >=100 nt',reference_utr3='Subset assigned to endogenous-reference 3UTR features >=100 nt',
    common_four_libraries=paste(format(length(common),big.mark=','),'3UTR features supported in all four CLIP libraries'))
  plot <- ggplot() + null_layers(p,cc) +
    geom_line(data=rr,aes(position,value,linetype=replicate),color='#D62728',alpha=.48,linewidth=.45) +
    geom_line(data=p,aes(position,observed,linetype='Replicate mean'),color='#D62728',linewidth=.9) +
    facet_wrap(~construct,nrow=1,labeller=as_labeller(c('324'='324 CLIP','328'='328 CLIP'))) +
    scale_linetype_manual(name=NULL,values=c('Replicate mean'='solid','Replicate A'='dashed','Replicate B'='dotted')) + sty +
    labs(title='Isoform CLIP read positions: equal 3UTR weight',subtitle=ttl,x='3UTR position (%)',y='Weighted fraction per 1% bin',
      caption='Red: read centers, thin replicate curves and thick average. Cyan: 500 null averages; black: null median.\nEach 3UTR contributes total weight one within a library; libraries averaged equally. No SMI normalization.\nRead centers reconstructed from BAM alignment exports; no smoothing. Null uses the prior 12.5-nt end-margin convention.')
  savefig(plot,paste0('isoform_read_positions_',cohort),12,6.8)
}
print(su[cohort=='all_assigned' & zone %in% c('first10','first20')])
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
inputs <- c('data/results_round2/features.tsv.gz','data/results_round2/sample_metadata.tsv',
  'output/followup_analysis/tables/summit_normalized_positions.tsv.gz',
  'data/results/peak_analysis/figures/all_nonzero_peaks/transcript_region_bias/peak_transcript_region_annotations.tsv',
  file.path('data/results_round2',meta$sample_id,'read_annotations.tsv.gz'))
fwrite(data.table(path=inputs,md5=unname(tools::md5sum(inputs))),file.path(out,'input_manifest.tsv'),sep='\t')
