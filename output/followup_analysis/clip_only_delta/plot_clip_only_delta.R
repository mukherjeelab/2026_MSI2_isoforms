suppressPackageStartupMessages({library(data.table);library(ggplot2);library(preprocessCore)})
out <- 'output/followup_analysis/clip_only_delta'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
base <- 'data/results/peak_analysis'
paths <- c(counts=file.path(base,'peak_count_matrix.tsv'),metadata='data/results/sample_metadata.tsv',
  library_sizes=file.path(base,'bam_library_sizes.tsv'),qc=file.path(base,'figures/QC/VST_counts_for_QC.tsv'),
  annotations=file.path(base,'figures/all_nonzero_peaks/transcript_region_bias/peak_transcript_region_annotations.tsv'))
before <- tools::md5sum(paths)
counts <- fread(paths['counts']); meta <- fread(paths['metadata'])
libs <- fread(paths['library_sizes']); annotations <- fread(paths['annotations']);qc <- fread(paths['qc'])
all_samples <- setdiff(names(counts),'peak_id')
stopifnot(length(all_samples)==8,!anyDuplicated(counts$peak_id),!anyDuplicated(annotations$peak_id))
total <- rowSums(as.matrix(counts[,..all_samples]))
ids <- list(peaks_16046=counts$peak_id[total>0],peaks_4385=counts$peak_id[total>=10])
stopifnot(length(ids[[1]])==16046,length(ids[[2]])==4385,setequal(ids[[2]],qc$peak_id))
clip <- meta[assay=='CLIP'][order(construct,biological_replicate)]
stopifnot(nrow(clip)==4,all(table(clip$construct)==2),setequal(clip$construct,c(324,328)))
clip <- merge(clip,libs,by.x='sample_id',by.y='sample',sort=FALSE)
setorder(clip,construct,biological_replicate)
clip_samples <- clip$sample_id
stopifnot(all(clip_samples %in% all_samples),all(clip$primary_mapped_alignments>0))
fwrite(clip[,.(sample_id,construct,assay,biological_replicate,primary_mapped_alignments)],
       file.path(out,'CLIP_library_metadata.tsv'),sep='\t')
levels_region <- c("5' UTR",'CDS',"3' UTR",'Intron','Other')
colors <- c("5' UTR"='#4C78A8',CDS='#F58518',"3' UTR"='#54A24B',Intron='#B279A2',Other='#79706E')
all_values <- list(); sample_stats <- list(); diagnostics <- list()
for(dataset in names(ids)) {
  rows <- match(ids[[dataset]],counts$peak_id)
  raw <- as.matrix(counts[rows,..clip_samples]);storage.mode(raw)<-'double'
  rownames(raw) <- ids[[dataset]]
  cpm <- sweep(raw,2,clip$primary_mapped_alignments,'/')*1e6
  logcpm <- log2(cpm+.5)
  qn <- preprocessCore::normalize.quantiles(logcpm,copy=TRUE)
  dimnames(qn) <- dimnames(logcpm)
  stopifnot(all(is.finite(qn)),identical(dim(qn),dim(logcpm)))
  # Independent rank-target check on a tie-free toy matrix.
  toy <- matrix(c(4,1,9,3,2,8,6,5,12),3,3)
  target <- rowMeans(apply(toy,2,sort)); expected <- apply(toy,2,function(x) target[rank(x)])
  stopifnot(max(abs(preprocessCore::normalize.quantiles(toy)-expected))<1e-12)
  # Ties in sparse columns are expected; check that tied input values remain tied after QN.
  for(j in 1:4) {
    tied <- data.table(x=logcpm[,j],y=qn[,j])[,.(spread=max(y)-min(y)),by=x]
    stopifnot(all(tied$spread<1e-10))
  }
  for(nm in c('raw_CLIP_counts','CLIP_CPM','log2_CLIP_CPM','QN_log2_CLIP_CPM')) {
    mat <- switch(nm,raw_CLIP_counts=raw,CLIP_CPM=cpm,log2_CLIP_CPM=logcpm,QN_log2_CLIP_CPM=qn)
    fwrite(cbind(data.table(peak_id=rownames(mat)),as.data.table(mat)),
           file.path(out,paste0(dataset,'_',nm,'.tsv.gz')),sep='\t')
  }
  annotation <- annotations[match(ids[[dataset]],peak_id),.(peak_id,transcript_region,transcript_region_detail,
    summit_chromosome,summit_start_0based,summit_strand,representative_gene_id)]
  stopifnot(identical(annotation$peak_id,ids[[dataset]]),all(annotation$transcript_region %in% levels_region))
  for(method in c('CPM_only','quantile_normalized')) {
    m <- if(method=='CPM_only') logcpm else qn
    v <- copy(annotation)
    v[,`:=`(dataset=dataset,method=method,
      mean_324_log2_CLIP=rowMeans(m[,clip$construct==324,drop=FALSE]),
      mean_328_log2_CLIP=rowMeans(m[,clip$construct==328,drop=FALSE]),
      total_CLIP_counts=rowSums(raw),zero_in_all_four_CLIP=rowSums(raw)==0)]
    v[,delta:=mean_324_log2_CLIP-mean_328_log2_CLIP]
    # Independently reconstruct each sample-pair contrast and average across replicates.
    pair_delta <- sapply(c('A','B'),function(r) {
      m[,clip$construct==324 & clip$biological_replicate==r] -
        m[,clip$construct==328 & clip$biological_replicate==r]
    })
    stopifnot(max(abs(v$delta-rowMeans(pair_delta)))<1e-12)
    fwrite(v,file.path(out,paste0(dataset,'_',method,'_peak_deltas.tsv.gz')),sep='\t')
    all_values[[paste(dataset,method)]] <- v
    diagnostics[[paste(dataset,method)]] <- v[,.(dataset=dataset[1],method=method[1],n_peaks=.N,
      n_zero_all_CLIP=sum(zero_in_all_four_CLIP),
      n_zero_all_CLIP_nonzero_delta=sum(zero_in_all_four_CLIP & abs(delta)>1e-12),
      overall_median=median(delta),mean_delta=mean(delta),positive=sum(delta>1e-12),
      negative=sum(delta< -1e-12),zero=sum(abs(delta)<=1e-12))]
    sample_stats[[paste(dataset,method)]] <- rbindlist(lapply(1:4,function(j) data.table(
      dataset=dataset,method=method,sample_id=clip$sample_id[j],construct=clip$construct[j],
      mean_log2=mean(m[,j]),median_log2=median(m[,j]),q25=quantile(m[,j],.25),q75=quantile(m[,j],.75))))
  }
}
values <- rbindlist(all_values)
summary <- values[,.(n_peaks=.N,median_delta=median(delta),q25_delta=quantile(delta,.25),
  q75_delta=quantile(delta,.75),mean_delta=mean(delta),n_zero_all_CLIP=sum(zero_in_all_four_CLIP)),
  by=.(dataset,method,transcript_region)]
fwrite(summary,file.path(out,'CLIP_delta_region_summary.tsv'),sep='\t')
fwrite(rbindlist(diagnostics),file.path(out,'CLIP_delta_diagnostics.tsv'),sep='\t')
fwrite(rbindlist(sample_stats),file.path(out,'normalization_sample_summary.tsv'),sep='\t')
axis_limit <- ceiling(max(abs(values$delta)))+2
for(dataset_name in names(ids)) for(method_name in c('CPM_only','quantile_normalized')) {
  v <- values[dataset==dataset_name & method==method_name]
  ss <- summary[dataset==dataset_name & method==method_name]
  region_labels <- setNames(sprintf('%s  (n = %s; median = %.3f)',ss$transcript_region,
    format(ss$n_peaks,big.mark=',',trim=TRUE),ss$median_delta),ss$transcript_region)
  p <- ggplot(v,aes(x=factor(transcript_region,levels=levels_region),y=delta,fill=transcript_region))+
    geom_hline(yintercept=0,color='#9E2A2B',linewidth=.55)+
    geom_violin(scale='width',trim=FALSE,color='grey45',linewidth=.35,alpha=.75)+
    geom_boxplot(width=.13,outlier.shape=NA,fill='white',color='grey20',linewidth=.4)+
    stat_summary(fun=median,geom='point',shape=21,size=2.2,fill='black',color='white',stroke=.4)+
    scale_fill_manual(values=colors,guide='none')+
    scale_x_discrete(labels=region_labels,limits=levels_region)+
    coord_flip(ylim=c(-axis_limit,axis_limit))+
    theme_classic(base_size=11)+theme(plot.title=element_text(face='bold',size=14),
      plot.subtitle=element_text(color='grey30'),plot.caption=element_text(color='grey35',hjust=0),
      axis.title=element_text(face='bold'),panel.grid.major=element_line(color='grey90',linewidth=.3))+
    labs(title='CLIP-only 324WT-versus-328WT delta by transcript region',
      subtitle=paste0(format(nrow(v),big.mark=','),' peaks | ',
        if(method_name=='CPM_only') 'Library-depth CPM; no quantile normalization' else 'Quantile-normalized log2 CPM; four CLIP libraries only'),
      x=NULL,y='Delta (324 - 328; positive favors 324WT; negative favors 328WT)',
      caption=paste('Delta = mean log2(CLIP CPM + 0.5) in 324 - mean log2(CLIP CPM + 0.5) in 328;',
        if(method_name=='CPM_only') 'no SMI subtraction. Two replicates per construct.' else 'quantile normalization precedes replicate averaging. No SMI subtraction.',
        'Boxplots: median and IQR. Labels: peak count and median. Shared x-axis range across all four plots.',sep='\n'))
  stem <- file.path(out,paste0('CLIP_delta_',dataset_name,'_',method_name))
  ggsave(paste0(stem,'.png'),p,width=12,height=7.3,dpi=250,bg='white')
  ggsave(paste0(stem,'.pdf'),p,width=12,height=7.3,bg='white',useDingbats=FALSE)
}
stopifnot(identical(before,tools::md5sum(paths)))
fwrite(data.table(path=unname(paths),md5=unname(before)),file.path(out,'input_provenance.tsv'),sep='\t')
writeLines(capture.output(sessionInfo()),file.path(out,'session_info.txt'))
print(summary)
cat('Four plots saved; original sources unchanged.\n')
