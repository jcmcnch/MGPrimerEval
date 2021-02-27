library(ggplot2)

all.summary <- read.delim(snakemake@input[[1]], header=FALSE)
euk.all <- all.summary[with(all.summary, V3=="EUK" & V5!="6-mismatch" & V4!="27F" & V4!="338R" & V4!="341F" & V4!="806RB" & V4!="785R"), ]
euk.all$V4=factor(euk.all$V4, levels=c("515Y","926R","V4F","V4R","V4RB","926wF","1392R","1389F","1510R"))

euk.all$V10[euk.all$V4=="515Y"] <- "V4-V5"
euk.all$V10[euk.all$V4=="926R"] <- "V4-V5"

euk.all$V10[euk.all$V4=="V4F"] <- "V4"
euk.all$V10[euk.all$V4=="V4R"] <- "V4"
euk.all$V10[euk.all$V4=="V4RB"] <- "V4"

euk.all$V10[euk.all$V4=="926wF"] <- "V6-V8"
euk.all$V10[euk.all$V4=="1392R"] <- "V6-V8"

euk.all$V10[euk.all$V4=="1389F"] <- "V9"
euk.all$V10[euk.all$V4=="1510R"] <- "V9"

strDatasetName = snakemake@params[[1]]
ylabelEUK = paste("Coverage of ", strDatasetName, " SSU rRNA (Eukaryotes)",sep="") 
ylabelBACT = paste("Coverage of ", strDatasetName, " SSU rRNA (Bacteria)",sep="")
ylabelCYANO = paste("Coverage of ", strDatasetName, " SSU rRNA (Cyano + plastid)",sep="")
ylabelARCH = paste("Coverage of ", strDatasetName, " SSU rRNA (Archaea)",sep="")

#Euks (18S)
p <- ggplot(na.omit(euk.all), aes(x=V4, y=V9, fill=V5)) +
geom_boxplot() + labs(x="Primer", y=ylabelEUK, fill="Mismatches") +
facet_grid(. ~ V10, scales="free")
p<-p+scale_y_continuous(breaks=seq(0,1,0.1))
p<-p+ theme_bw() + theme(plot.title = element_text(hjust = 0.5), panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
                        panel.grid.major.y = element_line(linetype="dashed", colour = "grey", size=0.5),
                        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), 
                        axis.title = element_text(size = rel(1.5)), #axis.line = element_line(colour = "black"),
                        axis.text = element_text(size = rel(1.3)), strip.text = element_text(size = rel(1.3)),
                        legend.text = element_text(size = rel(1.3)), legend.title = element_text(size = rel(1.3),hjust = 0.5),
                        legend.box.background = element_rect(colour = "black", fill=NA, size=0.5) )
ggsave(file=snakemake@output[[1]], plot=p, width=11, height=8)

#Bacteria (non-cyano)

bact.all <- all.summary[with(all.summary, V3=="BACT-NON-CYANO" & V6>=10 & V5!="6-mismatch" & V4!="V4F" & V4!="1389F" & V4!="1510R" & V4!="V4RB" & V4!="V4R"), ]
bact.all$V4=factor(bact.all$V4, levels=c("27F","338R","341F","785R","515Y","806RB","926R","V4F","V4R","V4RB","926wF","1392R"))

bact.all$V10[bact.all$V4=="515Y"] <- "V4 / V4-V5"
bact.all$V10[bact.all$V4=="806RB"] <- "V4 / V4-V5"
bact.all$V10[bact.all$V4=="926R"] <- "V4 / V4-V5"

bact.all$V10[bact.all$V4=="27F"] <- "V1-V2"
bact.all$V10[bact.all$V4=="338R"] <- "V1-V2"

bact.all$V10[bact.all$V4=="341F"] <- "V3-V4"
bact.all$V10[bact.all$V4=="785R"] <- "V3-V4"

bact.all$V10[bact.all$V4=="926wF"] <- "V6-V8"
bact.all$V10[bact.all$V4=="1392R"] <- "V6-V8"

p <- ggplot(na.omit(bact.all), aes(x=V4, y=V9, fill=V5)) +
geom_boxplot() + labs(x="Primer", y=ylabelBACT, fill="Mismatches") +
facet_grid(. ~ V10, scales="free")
p<-p+scale_y_continuous(breaks=seq(0,1,0.1))
p<-p+ theme_bw() + theme(plot.title = element_text(hjust = 0.5), panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
                        panel.grid.major.y = element_line(linetype="dashed", colour = "grey", size=0.5),
                        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), 
                        axis.title = element_text(size = rel(1.5)), #axis.line = element_line(colour = "black"),
                        axis.text = element_text(size = rel(1.3)), strip.text = element_text(size = rel(1.3)),
                        legend.text = element_text(size = rel(1.3)), legend.title = element_text(size = rel(1.3),hjust = 0.5),
                        legend.box.background = element_rect(colour = "black", fill=NA, size=0.5) )
ggsave(file=snakemake@output[[2]], plot=p, width=11, height=8)

#Cyano + plastids
cyano.all <- all.summary[with(all.summary, V3=="BACT-CYANO" & V6>=10 & V5!="6-mismatch" & V4!="V4F" & V4!="V4R" & V4!="V4RB" & V4!="1389F" & V4!="1510R"), ]
cyano.all$V4=factor(cyano.all$V4, levels=c("27F","338R","341F","785R","515Y","806RB","926R","V4F","V4R","V4RB","926wF","1392R"))

cyano.all$V10[cyano.all$V4=="515Y"] <- "V4 / V4-V5"
cyano.all$V10[cyano.all$V4=="806RB"] <- "V4 / V4-V5"
cyano.all$V10[cyano.all$V4=="926R"] <- "V4 / V4-V5"

cyano.all$V10[cyano.all$V4=="27F"] <- "V1-V2"
cyano.all$V10[cyano.all$V4=="338R"] <- "V1-V2"

cyano.all$V10[cyano.all$V4=="341F"] <- "V3-V4"
cyano.all$V10[cyano.all$V4=="785R"] <- "V3-V4"

cyano.all$V10[cyano.all$V4=="926wF"] <- "V6-V8"
cyano.all$V10[cyano.all$V4=="1392R"] <- "V6-V8"

cyano.all$V10[cyano.all$V4=="1389F"] <- "V9"
cyano.all$V10[cyano.all$V4=="1510R"] <- "V9"

p <- ggplot(na.omit(cyano.all), aes(x=V4, y=V9, fill=V5)) +
geom_boxplot() + labs(x="Primer", y=ylabelCYANO, fill="Mismatches") +
facet_grid(. ~ V10, scales="free")
p<-p+scale_y_continuous(breaks=seq(0,1,0.1))
p<-p+ theme_bw() + theme(plot.title = element_text(hjust = 0.5), panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
                        panel.grid.major.y = element_line(linetype="dashed", colour = "grey", size=0.5),
                        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), 
                        axis.title = element_text(size = rel(1.5)), #axis.line = element_line(colour = "black"),
                        axis.text = element_text(size = rel(1.3)), strip.text = element_text(size = rel(1.3)),
                        legend.text = element_text(size = rel(1.3)), legend.title = element_text(size = rel(1.3),hjust = 0.5),
                        legend.box.background = element_rect(colour = "black", fill=NA, size=0.5) )
ggsave(file=snakemake@output[[3]], plot=p, width=10, height=8)

#Archaea
arch.all <- all.summary[with(all.summary, V3=="ARCH" & V6>=10 & V9!="NA" & V5!="6-mismatch" & V4!="1389F" & V4!="27F" & V4!="338R" & V4!="1510R" & V4!="V4F" & V4!="V4R" & V4!="V4RB"), ]
arch.all$V4=factor(arch.all$V4, levels=c("341F","785R","515Y","806RB","926R","926wF","1392R", "807F","1050R"))

arch.all$V10[arch.all$V4=="515Y"] <- "V4 / V4-V5"
arch.all$V10[arch.all$V4=="806RB"] <- "V4 / V4-V5"
arch.all$V10[arch.all$V4=="926R"] <- "V4 / V4-V5"

arch.all$V10[arch.all$V4=="341F"] <- "V3-V4"
arch.all$V10[arch.all$V4=="785R"] <- "V3-V4"

arch.all$V10[arch.all$V4=="926wF"] <- "V6-V8"
arch.all$V10[arch.all$V4=="1392R"] <- "V6-V8"

p <- ggplot(na.omit(arch.all), aes(x=V4, y=V9, fill=V5)) +
geom_boxplot() + labs(x="Primer", y=ylabelARCH, fill="Mismatches") +
facet_grid(. ~ V10, scales="free")
p<-p+scale_y_continuous(breaks=seq(0,1,0.1))
p<-p+ theme_bw() + theme(plot.title = element_text(hjust = 0.5), panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
                        panel.grid.major.y = element_line(linetype="dashed", colour = "grey", size=0.5),
                        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), 
                        axis.title = element_text(size = rel(1.5)), #axis.line = element_line(colour = "black"),
                        axis.text = element_text(size = rel(1.3)), strip.text = element_text(size = rel(1.3)),
                        legend.text = element_text(size = rel(1.3)), legend.title = element_text(size = rel(1.3),hjust = 0.5),
                        legend.box.background = element_rect(colour = "black", fill=NA, size=0.5) )
ggsave(file=snakemake@output[[4]], width=10, height=8)
