suppressPackageStartupMessages({library(data.table);library(ggplot2)})
p='output/followup_analysis/'
savep=function(x,n,w=10,h=7){for(e in c('png','pdf'))ggsave(paste0(p,'figures/',n,'.',e),x,width=w,height=h,dpi=220)}
a=fread(paste0(p,'tables/read_annotation_summary.tsv'));a[,group:=paste(construct,assay)];a[,library:=paste(group,biological_replicate)]
savep(ggplot(a,aes(library,fractional_fraction,fill=category))+geom_col()+theme_bw()+theme(axis.text.x=element_text(angle=45,hjust=1))+labs(x=NULL,y='Fractional weight / eligible aligned reads',title='Aligned-read annotation: one total weight per read'),'read_fractional_annotation')
savep(ggplot(a,aes(group,nonexclusive_fraction,color=biological_replicate))+geom_point(size=3)+facet_wrap(~category,scales='free_y')+theme_bw()+labs(x=NULL,y='Reads overlapping category / eligible reads',title='Nonexclusive annotation: biological replicates'),'read_nonexclusive_annotation')
x=fread(paste0(p,'tables/downsampling_spatial_iterations.tsv.gz'))
full=x[iteration==0];sm=x[iteration>0,.(median=median(value),low=quantile(value,.025),high=quantile(value,.975)),by=.(profile,bin,position,sample_id,construct,assay,replicate)]
fwrite(sm,paste0(p,'tables/downsampling_spatial_summary.tsv'),sep='\t')
savep(ggplot(sm[profile%in%c('donor_into_intron','acceptor_into_intron')],aes(position,median,color=factor(construct),fill=factor(construct),linetype=replicate))+geom_ribbon(aes(ymin=low,ymax=high,group=interaction(construct,replicate)),alpha=.15,color=NA)+geom_line()+facet_grid(assay~profile)+theme_bw()+labs(x='Distance into intron (nt)',y='Mean aligned-base depth per intron per million sampled reads',title='50 exact read subsamples: intron ends',color='Construct',fill='Construct'),'intron_inward_resampling')
savep(ggplot(full[grepl('metaintron',profile)],aes(position,value,color=factor(construct),linetype=replicate))+geom_line()+facet_grid(assay~profile,scales='free_y')+theme_bw()+labs(x='Intron position (%)',y='Mean depth per intron per million eligible reads',title='Equal intron weighting; full-depth biological replicates'),'metaintron_read_quartiles',14,7)
# Additive differences avoid arbitrary ratios in sparse/zero profile bins.
g=x[,.(value=mean(value)),by=.(profile,bin,position,construct,assay,iteration)]
w=dcast(g,profile+bin+position+iteration~construct+assay,value.var='value');w[,dd:=(`324_CLIP`-`324_SMI`)-(`328_CLIP`-`328_SMI`)]
d=w[iteration>0,.(median=median(dd),low=quantile(dd,.025),high=quantile(dd,.975),fraction_positive=mean(dd>0)),by=.(profile,bin,position)]
fwrite(w,paste0(p,'tables/spatial_components_and_delta_iterations.tsv.gz'),sep='\t');fwrite(d,paste0(p,'tables/spatial_delta_summary.tsv'),sep='\t')
savep(ggplot(d[profile%in%c('metaintron','donor_into_intron','acceptor_into_intron')],aes(position,median))+geom_hline(yintercept=0,lty=2)+geom_ribbon(aes(ymin=low,ymax=high),fill='#5f9da8',alpha=.3)+geom_line()+facet_wrap(~profile,scales='free')+theme_bw()+labs(x='Relative position (%) or inward distance (nt)',y='(324 CLIP−SMI) − (328 CLIP−SMI), per million',title='Spatial construct contrast: median and 95% resampling interval'),'spatial_construct_sensitivity',12,5)
r=fread(paste0(p,'tables/downsampling_region_iterations.tsv'));z=r[iteration>0,.(median=median(median),low=quantile(median,.025),high=quantile(median,.975),fraction_negative=mean(median<0)),by=region];fwrite(z,paste0(p,'tables/downsampling_region_summary.tsv'),sep='\t')
savep(ggplot(r[iteration>0],aes(region,median))+geom_boxplot()+geom_hline(yintercept=0,lty=2)+theme_bw()+labs(x=NULL,y='Median peak delta-delta within each subsample',title='50 equal-depth subsamples; fixed16046 peak universe'),'downsampling_regions')
# Predefined zones; paired within-intron biological interpretation remains descriptive with two replicates.
full[,zone:=fifelse(bin<10,'first10',fifelse(bin>=90,'last10','middle80'))]
z=full[grepl('metaintron',profile),.(mean_density=mean(value)),by=.(profile,sample_id,construct,assay,replicate,zone)]
fwrite(z,paste0(p,'tables/metaintron_read_zone_summary.tsv'),sep='\t')
