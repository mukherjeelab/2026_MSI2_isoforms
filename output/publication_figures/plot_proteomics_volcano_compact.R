suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ggrepel)})
out <- 'output/publication_figures'
label_genes <- c('SF3A2','SF3A3','WDR33','ALYREF','EIF3G','EIF3D',
                 'PABPC4','RPL13','RPS13','RPL36A')
d <- fread('output/proteomics_hpa/isoform_classification_proteins.tsv')
d <- d[(cluster %in% c(2,3) | source_symbol %in% label_genes) & !is.na(source_symbol) & source_symbol!='' &
         is.finite(fc_328_324) & fc_328_324>0 & is.finite(comparison_p) & comparison_p>0]
d[, `:=`(gene=source_symbol,x=-log2(fc_328_324),y=-log10(comparison_p))]
d[, direction := fifelse(x>=1 & comparison_p<.5,'324',
                         fifelse(x<= -1 & comparison_p<.5,'328','Other'))]
d[, selected := gene %in% label_genes & !(gene=='ALYREF' & direction=='328')]
lab <- d[selected==TRUE]
stopifnot(setequal(lab$gene,label_genes),nrow(lab)==10)
# Nudges encourage labels to use open space above the crowded lower points.
p <- ggplot(d,aes(x,y)) +
  geom_point(color='grey75',alpha=.65,size=1.6) +
  geom_vline(xintercept=c(-1,1),linetype='dashed',color='grey45',linewidth=.4) +
  geom_hline(yintercept=-log10(.5),linetype='dashed',color='grey45',linewidth=.4) +
  geom_point(data=lab,aes(color=direction),size=2.1) +
  geom_text_repel(data=lab,aes(label=gene),size=2.35,family='Helvetica',
    seed=328, max.overlaps=Inf,max.iter=100000,max.time=10,
    box.padding=.3,point.padding=.15,min.segment.length=0,
    segment.color='grey55',segment.size=.25,force=2,force_pull=.15,
    nudge_y=ifelse(lab$y<2,1,0),ylim=c(-.15,6.25),
    show.legend=FALSE) +
  annotate('text',x=-7,y=6.55,label='MSI2-328 enriched',size=3,family='Helvetica') +
  annotate('text',x=8,y=6.55,label='MSI2-324 enriched',size=3,family='Helvetica') +
  scale_color_manual(values=c('324'='#D55E00','328'='#0072B2','Other'='black')) +
  scale_x_continuous(breaks=seq(-10,10,5),expand=expansion(mult=.045)) +
  scale_y_continuous(breaks=c(0,2,4,6),limits=c(-.35,6.75),expand=expansion(mult=0)) +
  labs(x=expression(log[2]*' fold change: MSI2-324 / MSI2-328'),y=expression(-log[10](p))) +
  theme_minimal(base_size=10,base_family='Helvetica') +
  theme(legend.position='none',panel.grid.minor=element_blank(),
        panel.grid.major=element_line(color='#E8E8E8',linewidth=.4),
        axis.text=element_text(color='black'),plot.margin=margin(5,7,5,5))
stem <- file.path(out,'proteomics_volcano_compact')
ggsave(paste0(stem,'.pdf'),p,width=5.2,height=4.5,device='pdf',useDingbats=FALSE,bg='white')
# Render the PDF for identical label positions in the PNG.
fontconf <- tempfile(fileext='.conf')
writeLines(paste0('<?xml version="1.0"?><!DOCTYPE fontconfig SYSTEM "fonts.dtd"><fontconfig>',
  '<dir>/System/Library/Fonts</dir><dir>/Library/Fonts</dir><cachedir>',
  tempdir(),'</cachedir></fontconfig>'),fontconf)
status <- system2('pdftoppm',c('-r','600','-singlefile','-png',
  shQuote(paste0(stem,'.pdf')),shQuote(stem)),
  env=paste0('FONTCONFIG_FILE=',shQuote(fontconf)))
stopifnot(status==0)
fwrite(d,paste0(stem,'_source_data.tsv'),sep='\t')
fwrite(lab[,.(gene,x,y,direction)],paste0(stem,'_labels.tsv'),sep='\t')
cat(nrow(d),'points;',nrow(lab),'labels\n')
