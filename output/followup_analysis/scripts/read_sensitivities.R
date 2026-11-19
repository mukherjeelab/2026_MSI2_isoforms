suppressPackageStartupMessages({library(data.table);library(GenomicRanges);library(GenomicAlignments);library(Matrix);library(jsonlite)})
options(warn=1); set.seed(260921)
out='output/followup_analysis'; dir.create(paste0(out,'/cache'),showWarnings=FALSE)
meta=fread('data/results_round2/sample_metadata.tsv'); depth=c(148054,257457,173852,320014,93659,571306,103812,512565);stopifnot(nrow(meta)==8)
meta[,original_depth:=depth];meta[,target_depth:=min(original_depth),by=assay];meta[,seed:=260921L];meta[,iterations:=50L];fwrite(meta,paste0(out,'/tables/downsampling_depths.tsv'),sep='\t')
f=fread('data/results_round2/features.tsv.gz'); intr=f[kind=='intron' & length_nt>=100]; seg=lapply(intr$segments_0based_half_open,fromJSON);intr[,start:=vapply(seg,function(x)x[1,1],numeric(1))];intr[,end:=vapply(seg,function(x)x[1,2],numeric(1))];intr[,quartile:=pmin(4L,ceiling(frank(length_nt,ties.method='average')/.N*4))]
ig=GRanges(intr$chrom,IRanges(intr$start+1,intr$end),intr$strand)
fwrite(intr[,.(feature_id,chrom,strand,start,end,length_nt,quartile)],paste0(out,'/tables/eligible_introns.tsv.gz'),sep='\t')
N=nrow(intr); nq=tabulate(intr$quartile,4); longN=sum(intr$length_nt>=2000)
profilemeta=rbindlist(list(data.table(profile='metaintron',bin=0:99,position=(0:99+.5),denominator=N),rbindlist(lapply(1:4,function(q)data.table(profile=paste0('metaintron_Q',q),bin=0:99,position=(0:99+.5),denominator=nq[q]))),data.table(profile='donor_into_intron',bin=0:99,position=(0:99+.5)*10,denominator=longN),data.table(profile='acceptor_into_intron',bin=0:99,position=(0:99+.5)*10,denominator=longN)))
fwrite(profilemeta,paste0(out,'/tables/spatial_denominators.tsv'),sep='\t')
cnt=fread('data/results/peak_analysis/peak_count_matrix.tsv'); ids=cnt[rowSums(as.matrix(cnt[,-1]))>0,peak_id];stopifnot(length(ids)==16046)
anno=fread('data/results/peak_analysis/figures/all_nonzero_peaks/transcript_region_bias/peak_transcript_region_annotations.tsv'); regions=anno$transcript_region[match(ids,anno$peak_id)]
counts=array(0,c(16046,8,51)); spatial=array(0,c(700,8,51)); annotations=list();supports=list()
for(s in 1:8){
 sid=meta$sample_id[s];message('Reading ',sid);d=fread(paste0('data/results_round2/',sid,'/read_annotations.tsv.gz')); n=nrow(d);stopifnot(n==depth[s]);
 # Each exported row is one eligible BAM alignment. Classification used CIGAR M/= /X blocks upstream.
 cats=strsplit(d$compatible_classes,';',fixed=TRUE); cats=lapply(seq_along(cats),function(i){z=cats[[i]];z=z[nzchar(z)];if(!length(z))z='other';if('other'%in%z && !nzchar(d$gene_ids[i]))z[z=='other']='intergenic';unique(z)})
 catsdt=data.table(read_index=rep(seq_len(n),lengths(cats)),category=unlist(cats));catsdt[,weight:=1/.N,by=read_index]
 a=catsdt[,.(nonexclusive_reads=.N,fractional_reads=sum(weight)),by=category];a[,`:=`(sample_id=sid,total_reads=n,nonexclusive_fraction=nonexclusive_reads/n,fractional_fraction=fractional_reads/n)];annotations[[s]]=a
 fwrite(catsdt,paste0(out,'/cache/',sid,'_fractional_read_memberships.tsv.gz'),sep='\t')
 pids=strsplit(d$overlapping_peak_ids,';',fixed=TRUE); pi=match(unlist(pids),ids); ri=rep(seq_len(n),lengths(pids)); keep=!is.na(pi); pm=sparseMatrix(i=ri[keep],j=pi[keep],x=1,dims=c(n,length(ids)))
 cache=paste0(out,'/cache/',sid,'_spatial_matrix.rds')
 if(file.exists(cache)){sm=readRDS(cache)}else{
 bl=cigarRangesAlongReferenceSpace(d$cigar,pos=d$position_1based,ops=c('M','=','X'),reduce.ranges=TRUE);rid=rep(seq_len(n),lengths(bl)); bg=GRanges(d$chrom[rid],unlist(bl),d$strand[rid]); h=findOverlaps(bg,ig,ignore.strand=FALSE);bi=queryHits(h);fi=subjectHits(h);a=pmax(start(bg)[bi]-1,intr$start[fi]);b=pmin(end(bg)[bi],intr$end[fi]);L=intr$length_nt[fi];minus=intr$strand[fi]=='-';aa=ifelse(minus,intr$end[fi]-b,a-intr$start[fi]);bb=ifelse(minus,intr$end[fi]-a,b-intr$start[fi]);r=rid[bi]
 k0=pmin(99L,ceiling((aa+1)*100/L)-1L);k1=pmin(99L,ceiling(bb*100/L)-1L);nn=k1-k0+1L;ix=rep(seq_along(nn),nn);kk=rep(k0,nn)+sequence(nn)-1L;le=floor(kk*L[ix]/100);re=floor((kk+1)*L[ix]/100);v=(pmin(bb[ix],re)-pmax(aa[ix],le))/(re-le)
 ii=c(r[ix],r[ix]);jj=c(kk+1L,intr$quartile[fi[ix]]*100L+kk+1L);vv=c(v,v)
 support=data.table(feature_id=intr$feature_id[fi],read_index=r);supports[[s]]=unique(support)[,.(reads=.N),by=feature_id][,sample_id:=sid]
 for(side in 1:2){ok=L>=2000;ta=if(side==1)aa[ok] else L[ok]-bb[ok];tb=if(side==1)bb[ok] else L[ok]-aa[ok];rr=r[ok];ok2=ta<1000 & tb>0;ta=pmax(0,ta[ok2]);tb=pmin(1000,tb[ok2]);rr=rr[ok2];lo=floor(ta/10);hi=ceiling(tb/10)-1L;nn=hi-lo+1L;ix=rep(seq_along(nn),nn);kk=rep(lo,nn)+sequence(nn)-1L;vv=c(vv,(pmin(tb[ix],(kk+1)*10)-pmax(ta[ix],kk*10))/10);ii=c(ii,rr[ix]);jj=c(jj,500L+(side-1L)*100L+kk+1L)}
 sm=sparseMatrix(i=ii,j=jj,x=vv,dims=c(n,700));saveRDS(sm,cache);rm(bg,bl,h,ii,jj,vv);gc()
 }
 for(it in 1:51){sel=if(it==1)seq_len(n) else sample.int(n,meta$target_depth[s],replace=FALSE);counts[,s,it]=Matrix::colSums(pm[sel,,drop=FALSE]);spatial[,s,it]=Matrix::colSums(sm[sel,,drop=FALSE])/profilemeta$denominator/length(sel)*1e6}
 fwrite(data.table(peak_id=ids,counts[,,1][,s]),paste0(out,'/cache/',sid,'_reconstructed_peak_counts.tsv'),sep='\t');rm(d,pm,sm);gc()
}
saveRDS(list(counts=counts,spatial=spatial,meta=meta,ids=ids,regions=regions,profilemeta=profilemeta),paste0(out,'/cache/resampling_arrays.rds'))
a=merge(rbindlist(annotations),meta[,.(sample_id,construct,assay,biological_replicate)],by='sample_id');fwrite(a,paste0(out,'/tables/read_annotation_summary.tsv'),sep='\t');stopifnot(all(abs(a[,sum(fractional_reads),by=sample_id]$V1-depth)<1e-6))
if(length(supports))fwrite(rbindlist(supports),paste0(out,'/tables/intron_nonzero_support.tsv.gz'),sep='\t')
curves=list();scores=list();reg=list()
for(it in 1:51){dep=if(it==1)depth else meta$target_depth;cpm=sweep(counts[,,it],2,dep,'/')*1e6; e324=log2((rowMeans(cpm[,c(5,7)])+.5)/(rowMeans(cpm[,c(6,8)])+.5));e328=log2((rowMeans(cpm[,c(1,3)])+.5)/(rowMeans(cpm[,c(2,4)])+.5));dd=e324-e328;scores[[it]]=dd;reg[[it]]=data.table(region=regions,dd=dd)[,.(n=.N,median=median(dd),positive_fraction=mean(dd>0)),by=region][,iteration:=it-1L];curves[[it]]=rbindlist(lapply(1:8,function(s)cbind(profilemeta,sample_id=meta$sample_id[s],construct=meta$construct[s],assay=meta$assay[s],replicate=meta$biological_replicate[s],iteration=it-1L,value=spatial[,s,it])))}
fwrite(rbindlist(curves),paste0(out,'/tables/downsampling_spatial_iterations.tsv.gz'),sep='\t');fwrite(rbindlist(reg),paste0(out,'/tables/downsampling_region_iterations.tsv'),sep='\t')
scores=do.call(cbind,scores);fwrite(data.table(peak_id=ids,region=regions,full_depth_dd=scores[,1],resampled_median=apply(scores[,-1],1,median),low=apply(scores[,-1],1,quantile,.025),high=apply(scores[,-1],1,quantile,.975),fraction_same_direction=rowMeans(sign(scores[,-1])==sign(scores[,1])),nonzero_full=rowSums(counts[,,1])>0),paste0(out,'/tables/downsampling_peak_summary.tsv'),sep='\t')
