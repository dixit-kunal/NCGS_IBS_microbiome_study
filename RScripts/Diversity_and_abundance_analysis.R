# This script can be used for Alpha diversity, Beta diversity analysis and to generate bar plots for top genera associated with NCGS and IBS

# Load required packages
library("ggplot2"); packageVersion("ggplot2")
library(reshape2)
library(ggpubr)
library(RColorBrewer)
library(ggpubr)
library(ape)
library(vegan)

# Set working  directory
setwd("path")

# Load phyloseq object
ps <- readRDS("ps.rds")

# Load Metadata
samdf <- read.table(file = "metadata.txt", sep = "\t", header = T, row.names = 1, check.names = F)
row.names(samdf)
sample_data(ps) <- samdf
ps
sample_names(ps)

# Subset samples absed on sample types
ps.baseline <- subset_samples(physeq = ps, Time_point == 'Baseline')
ps.stool <- subset_samples(physeq = ps.baseline, Sample_Type == 'Stool')
ps.duo <- subset_samples(physeq = ps.baseline, Sample_Type == 'Duodenal')
ps.sig <- subset_samples(physeq = ps.baseline, Sample_Type == 'Sigmoid')
ps.ibs.st <- subset_samples(physeq = ps.stool, Samples == 'IBS(AGA_Negative)')
ps.ibs.duo <- subset_samples(physeq = ps.duo, Samples == 'IBS(AGA_Negative)')
ps.ibs.sig <- subset_samples(physeq = ps.sig, Samples == 'IBS(AGA_Negative)')
ps.ncgs.st <- subset_samples(physeq = ps.stool, Samples == 'NCGS')
ps.ncgs.duo <- subset_samples(physeq = ps.duo, Samples == 'NCGS')
ps.ncgs.sig <- subset_samples(physeq = ps.sig, Samples == 'NCGS')
ps.healthy<- subset_samples(physeq = ps.baseline, Sample_Category == 'Healthy_Control')
ps.ibs.ncgs.st <- subset_samples(physeq = ps.stool, Samples == 'IBS(AGA_Negative)' | Samples == 'NCGS')
ps.ibs.ncgs.duo <- subset_samples(physeq = ps.duo, Samples == 'IBS(AGA_Negative)' | Samples == 'NCGS')
ps.ibs.ncgs.sig <- subset_samples(physeq = ps.sig, Samples == 'IBS(AGA_Negative)' | Samples == 'NCGS')
ps.ibs.ncgs.h <- subset_samples(physeq = ps.stool, Sample_Category == 'IBS(AGA_Negative)' | Sample_Category == 'NCGS' | Sample_Category == 'Healthy_Control')
ps.aga.neg <- subset_samples(physeq = ps.stool, AGA_Serology == 'Negative')
ps.subtype <- subset_samples(physeq = ps.aga.neg, Diagnosis == 'IBS-D' | Diagnosis == 'IBS-C' | Diagnosis == 'IBS-M' | Diagnosis == 'IBS-U')
ps.ncgs <- subset_samples(physeq = ps, Sample_Category == 'NCGS')
ps.ncgs.stool <- subset_samples(physeq = ps.ncgs, Sample_Type == 'Whole_gut')
ps.ncgs.duo <- subset_samples(physeq = ps.ncgs, Sample_Type == 'Duodenal')
ps.ncgs.sig <- subset_samples(physeq = ps.ncgs, Sample_Type == 'Sigmoid')

############## Alpha diversity ##############
sampleCompare <- list(c('Baseline','Post_GFD'))
plot_richness(ps.subtype, x="Diagnosis", measures=c("Shannon", "Simpson","Observed", fill="Sample_Category"))+geom_boxplot(lwd=1, aes(fill=treatment))+ scale_fill_manual(values=c("#00CC33", "#FF3333", "#3399FF", "orange", "purple"))+ stat_compare_means(comparisons = sampleCompare, method = "wilcox.test", size=6.5)+ geom_point(size=3, alpha=2.5)+ theme_bw()+ theme(axis.text.x = element_text(angle = 90, size = 15, colour = "black", vjust = 0.2, hjust = 0.5, face = "bold"), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold") ,strip.text = element_text(size = 18, colour = "black", face = "bold"), legend.title = element_text(size = 20, colour = "black", face = "bold"), legend.text = element_text(size = 22, face = "bold", colour = "black"), axis.text.y = element_text(size = 16, colour = "black", face = "bold"))

#################### Beta diversity ###############
########## NMDS-Bray ###############
ps.prop <- transform_sample_counts(ps.ncgs.stool, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="Time_point", title="Bray NMDS",label = "Sample_ID") + geom_point(size=2, alpha=0.5)+ stat_ellipse(geom = "polygon", linetype = 2,alpha=0.2, aes(fill=Time_point))

###################### Top20 genera abundance bar plot ###############
# Define colourscheme
kdcol = c("#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD","#0072B2",
          "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", "#F0E442",
          "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "steelblue", "#CBD588", "#FB9A99")
# Plot
ps.biopsy.pg <- merge_samples(x = ps.biopsy.pg,group = 'Group')
sample_data(ps.biopsy.pg)[,'Group'] <- sample_names(ps.biopsy.pg)
temp.taxa <- as.data.frame(tax_table(ps.biopsy.pg))
temp.taxa$Genus <- gsub(pattern = '_.*',replacement = '',x = temp.taxa$Genus,perl = T)
tax_table(ps.biopsy.pg) <- tax_table(as.matrix(temp.taxa))
ps.genus <- tax_glom(physeq = ps.biopsy.pg, taxrank = "Genus")
ps.genus.na <- subset_taxa(ps.genus, Genus!="NA")
top20 <- names(sort(taxa_sums(ps.genus.na), decreasing=TRUE))[1:20]
ps.top20 <- prune_taxa(top20, ps.genus.na)
ps.top20 <- transform_sample_counts(ps.top20, function(OTU) OTU/sum(OTU))
nb.cols <- length(taxa_names(ps.top20))
nb.cols
my_color = colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
ps.genus.melt <- psmelt(ps.top20)
ps.genus.filt <- ps.genus.melt[,"Abundance" > 0]
ggplot(ps.genus.filt, aes(x = Sample, y = Abundance, fill = Genus)) + geom_bar(stat = "identity", colour= "black")+ theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black", vjust = 0.2, hjust = 0.5, face = "bold"), axis.title.y = element_text(size = 20, face = "bold"),strip.text = element_text(size = 15, colour = "black", face = "bold"), legend.title = element_text(size = 20, colour = "black", face = "bold"), legend.text = element_text(size = 22, face = "bold", colour = "black"), axis.text.y = element_text(size = 16, colour = "black", face = "bold"))+ scale_y_continuous()+ labs(x="", y="Relative_Abundance(%)", fill="")+ scale_fill_manual(values = kdcol)

# Save data
save.image("diversity_analysis.RData")


