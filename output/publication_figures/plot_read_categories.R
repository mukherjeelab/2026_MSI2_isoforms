suppressPackageStartupMessages({library(data.table); library(ggplot2)})
# Run from the project root. Fractions retain the original all-read denominator.
input <- 'output/followup_analysis/read_annotation_comparisons/source_library_fractions.tsv'
out <- 'output/publication_figures'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
d <- fread(input)[assay == 'CLIP' & category %in% c('3UTR', 'intron')]
d[, isoform := factor(construct, levels=c(324,328))]
d[, region := factor(category, levels=c('3UTR','intron'), labels=c("3' UTR", 'Intron'))]
d[, x := as.numeric(region) + fifelse(construct == 324, -.19, .19)]
d[, point_x := x + fifelse(biological_replicate == 'A', -.045, .045)]
m <- d[, .(fraction=mean(fractional_fraction), n=.N), by=.(region,isoform,x)]
stopifnot(nrow(d)==8, all(m$n==2), all(is.finite(d$fractional_fraction)))
p <- ggplot(m, aes(x, fraction, fill=isoform)) +
  geom_col(width=.33, color='black', linewidth=.45) +
  geom_point(data=d, aes(x=point_x,y=fractional_fraction), shape=21,
             fill='white', color='black', size=1.65, stroke=.4) +
  scale_fill_manual(values=c('324'='#D55E00','328'='#0072B2'),
                    labels=c('MSI2-324','MSI2-328')) +
  scale_x_continuous(breaks=c(1,2), labels=c("3' UTR",'Intron'), limits=c(.5,2.5)) +
  scale_y_continuous(limits=c(0,.5), breaks=seq(0,.5,.1),
                     minor_breaks=seq(0,.5,.05), labels=scales::label_percent(accuracy=1),
                     expand=expansion(mult=0)) +
  labs(x=NULL, y='Fraction of aligned CLIP reads', fill=NULL) +
  theme_minimal(base_size=10, base_family='Helvetica') +
  theme(panel.grid.major=element_line(color='#E5E5E5',linewidth=.45),
        panel.grid.minor.y=element_line(color='#EEEEEE',linewidth=.3),
        panel.grid.minor.x=element_blank(),
        axis.text=element_text(color='black'),
        axis.title.y=element_text(margin=margin(r=6)),
        legend.position='top', legend.direction='vertical',
        legend.justification='left', legend.margin=margin(0,0,3,0),
        legend.key.size=grid::unit(3.4,'mm'), legend.spacing.y=grid::unit(0,'mm'),
        plot.margin=margin(5,7,5,5))
stem <- file.path(out,'read_categories_324_vs_328')
ggsave(paste0(stem,'.png'),p,width=3.35,height=3.5,dpi=600,bg='white')
ggsave(paste0(stem,'.pdf'),p,width=3.35,height=3.5,device='pdf',useDingbats=FALSE,bg='white')
fwrite(d,paste0(stem,'_source_data.tsv'),sep='\t')
