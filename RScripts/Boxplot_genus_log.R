# The script can be used to plot boxplots of specific features or taxa using a phyloseq object
# Load packages
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library(reshape2)
library(ggpubr)
library(RColorBrewer)
library(ggpubr)
library(ape)
library(vegan)

#set working directory
setwd("path")

# Load phyloseq object
ps <- readRDS(file = 'ps.rds')
ps

# Glom to genus level
ps.genera <- tax_glom(ps, taxrank="Genus")

meltp3<- psmelt(ps.genera)

sampleCompare <- list(c('NCGS','IBS'))

ggplot(meltp3,aes(x=Samples, y=log10(Abundance), fill=Samples))+ geom_boxplot()+ theme_bw()+ facet_wrap(~Genus, scale="free_x")+ stat_compare_means(comparisons = sampleCompare, method = "wilcox.test", size=5, label="p.signif")+ theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black", vjust = 0.2, hjust = 0.5, face = "bold"), axis.title.y = element_text(size = 15, face = "bold"),strip.text = element_text(size = 15, colour = "black", face = "bold"), legend.title = element_text(size = 18, colour = "black", face = "bold"), legend.text = element_text(size = 16, face = "bold", colour = "black"), axis.text.y = element_text(size = 16, colour = "black", face = "bold"))+ scale_y_continuous()+ scale_fill_manual(values = c("#3399FF", "orange", "#996600")) + geom_jitter(aes())
