suppressPackageStartupMessages({library(data.table);library(ggplot2)})
p='output/followup_analysis/raw_count_quantile/'
d=fread(paste0(p,'continuous_delta_delta_peaks.tsv'))
ref=d[transcript_region=='Intron',delta_delta]
res=rbindlist(lapply(c('transcript_region','transcript_region_detail'),function(cl)rbindlist(lapply(unique(d[[cl]]),function(z){x=d[get(cl)==z,delta_delta]; if(length(x)<20 || z=='Intron')return(NULL); w=wilcox.test(x,ref,exact=FALSE); k=suppressWarnings(ks.test(x,ref));data.table(grouping=cl,region=z,n=length(x),n_intron=length(ref),median=median(x),intron_median=median(ref),median_difference=median(x)-median(ref),wilcoxon_U=unname(w$statistic),p=w$p.value,cliffs_delta=2*unname(w$statistic)/(length(x)*length(ref))-1,ks_D=unname(k$statistic),ks_p=k$p.value)}))))
res[,q_BH:=p.adjust(p,'BH')]
fwrite(res,'output/followup_analysis/tables/region_vs_intron.tsv',sep='\t')
savep=function(x,n,w=9,h=6){for(e in c('png','pdf'))ggsave(paste0('output/followup_analysis/figures/',n,'.',e),x,width=w,height=h,dpi=220)}
savep(ggplot(res,aes(cliffs_delta,reorder(region,cliffs_delta),color=grouping))+geom_vline(xintercept=0,lty=2)+geom_point(size=3)+facet_wrap(~grouping,scales='free_y')+theme_bw()+labs(x="Cliff's delta relative to intron",y=NULL,title='Raw-count quantile normalization: regional effects',subtitle='Positive: more 324-favoring than intron; negative: more 328-favoring'),'region_intron_effects',11,6)
savep(ggplot(d,aes(delta_delta,color=transcript_region))+geom_density()+theme_bw()+labs(x='324 − 328 delta-delta (QN log2 raw counts)',y='Probability density'),'all_region_density')
savep(ggplot(d,aes(delta_delta,color=transcript_region))+stat_ecdf()+theme_bw()+labs(x='324 − 328 delta-delta (QN log2 raw counts)',y='Cumulative fraction'),'all_region_ecdf')
