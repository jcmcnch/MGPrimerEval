library("ggplot2")
taxCompDF <- read.delim(snakemake@input[[1]], header=TRUE)
lenTaxa <- length(unique(taxCompDF$Taxon))
p <- ggplot(taxCompDF, aes(fill=Taxon, y=Fractional_abundance, x=Hits_or_misses)) + 
   geom_bar(stat="identity", color="black") + xlab(NULL) + ylab("Fractional abundance (for taxa >= 0.01)") + 
   scale_y_continuous(limits = c(0,1), expand = c(0, 0)) + labs(title=snakemake@params[[1]]) +
   scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "dodgerblue3", "lightskyblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey"))
p <- p + guides(fill=guide_legend(nrow=lenTaxa, byrow=TRUE)) 
p <- p + theme_bw() + theme(
                        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), 
                        axis.title = element_text(size = rel(1.5)), plot.title = element_text(size=rel(2),hjust = 0.5),
                        axis.text = element_text(size = rel(1.3)), strip.text = element_text(size = rel(1.3)),
                        legend.text = element_text(size = rel(0.8)), legend.title = element_text(size = rel(1.5),hjust = 0.5),
                        legend.box.background = element_rect(colour = "black", fill=NA, size=0.5) )
ggsave(file=snakemake@output[[1]], plot=p, width=10, height=10)
