suppressPackageStartupMessages({library(data.table);library(Matrix)})
p='output/followup_analysis/';x=readRDS(paste0(p,'cache/resampling_arrays.rds'));orig=fread('data/results/peak_analysis/peak_count_matrix.tsv');orig=orig[match(x$ids,peak_id)];obs=x$counts[,,1];save=list()
for(s in 1:8){r2=fread(paste0('data/results_round2/',x$meta$sample_id[s],'/peak_counts.tsv.gz'));print(names(r2));break}
qc=data.table(sample_id=x$meta$sample_id,original_fixed_peak_counts=colSums(as.matrix(orig[,-1])),reconstructed_counts=colSums(obs),peaks_nonzero=colSums(obs>0),peak_count_entries_equal_original=colMeans(obs==as.matrix(orig[,-1])))
fwrite(qc,paste0(p,'tables/reconstructed_count_QC.tsv'),sep='\t')
# Support denominator includes all eligible introns, even zero-support ones.
f=fread(paste0(p,'tables/eligible_introns.tsv.gz'));sp=fread(paste0(p,'tables/intron_nonzero_support.tsv.gz'));q=merge(sp,f[,.(feature_id,quartile,length_nt)],by='feature_id');sum=q[,.(nonzero_introns=uniqueN(feature_id)),by=.(sample_id,quartile)];den=f[,.(eligible_introns=.N),by=quartile];sum=merge(sum,den,by='quartile');sum[,fraction_supported:=nonzero_introns/eligible_introns];fwrite(sum,paste0(p,'tables/intron_support_summary.tsv'),sep='\t')
fwrite(data.table(check=c('fixed_peak_N','reconstructed_any_nonzero_N','sampled_CLIP_depth','sampled_SMI_depth','iterations','seed'),value=c(length(x$ids),sum(rowSums(obs)>0),min(x$meta[assay=='CLIP',target_depth]),min(x$meta[assay=='SMI',target_depth]),50,260921)),paste0(p,'tables/analysis_checks.tsv'),sep='\t')
