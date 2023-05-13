# The script can be used for differential taxa analysis with Phyloseq object as a staring point.

# Load necessary environments
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library(reshape2)
library(ggpubr)
library(RColorBrewer)
library(ggpubr)
library(ape)
library(vegan)
library("DESeq2")

# Set working directory
setwd("C:\\Users\\ncmr\\Desktop\\Kunal_D\\NCGS\\Run_wise_analysis\\Phyloseq\\Analysis_Feb2022\\merged\\phyloseqs")

# Load Phyloseq
ps <- readRDS("merged.ps.rds")
samdf <- read.table(file = "merged.metadata.txt", sep = "\t", header = T, row.names = 1, check.names = F)
row.names(samdf)
sample_data(ps) <- samdf
ps
sample_names(ps)

# Subset phyloseq based on sample groups
ps.baseline <- subset_samples(physeq = ps, Time_point == 'Baseline')
ps.pg <- subset_samples(physeq = ps, Time_point == 'Post_GFD')
ps.stool <- subset_samples(physeq = ps.baseline, Sample_Type == 'Stool')
ps.duo <- subset_samples(physeq = ps.baseline, Sample_Type == 'Duodenal')
ps.sig <- subset_samples(physeq = ps.baseline, Sample_Type == 'Sigmoid')
ps.ibs.ncgs.stool <- subset_samples(physeq = ps.stool, Samples == 'IBS(AGA_Negative)' | Samples == 'NCGS')
ps.ibs.ncgs.duo <- subset_samples(physeq = ps.duo, Samples == 'IBS(AGA_Negative)' | Samples == 'NCGS')
ps.ibs.ncgs.sig <- subset_samples(physeq = ps.sig, Samples == 'IBS(AGA_Negative)' | Samples == 'NCGS')
ps.ibs.stool.all <- subset_samples(physeq = ps.stool, Sample_Category == 'IBS(AGA_Negative)')
ps.ibs.duo.all <- subset_samples(physeq = ps.duo, Sample_Category == 'IBS(AGA_Negative)')
ps.ibs.sig.all <- subset_samples(physeq = ps.sig, Samples == 'IBS(AGA_Negative)')
ps.healthy <- subset_samples(physeq = ps.stool, Sample_Category == 'Healthy_Control')
ps.ibs.ncgs.h <- subset_samples(physeq = ps.stool, Sample_Category == 'IBS(AGA_Negative)' | Sample_Category == 'NCGS' | Sample_Category == 'Healthy_Control')
ps.allibs.ncgs <- subset_samples(physeq = ps.stool, Sample_Category == 'IBS(AGA_Negative)' | Sample_Category == 'NCGS')
ps.ibs.h <- subset_samples(physeq = ps.stool, Sample_Category == 'IBS(AGA_Negative)' | Sample_Category == 'Healthy_Control')
ps.ncgs <- subset_samples(physeq = ps.baseline, Sample_Category == 'NCGS')
ps.ncgs.pg <- subset_samples(physeq = ps.pg, Sample_Category == 'NCGS')
ps.ncgs.stool <- subset_samples(physeq = ps.ncgs, Sample_Type == 'Stool')
ps.ncgs.duo <- subset_samples(physeq = ps.ncgs, Sample_Type == 'Duodenal')
ps.ncgs.sig <- subset_samples(physeq = ps.ncgs, Sample_Type == 'Sigmoid')


# Add +1 to remove 0 values
otu_added_1 <- as.matrix(ps@otu_table, counts)+1

#make new phyloseq by using the generated sequence file
taxa<- tax_table(ps)
ps_added_1<- phyloseq(otu_table(otu_added_1), tax_table(taxa))

samdf <- read.table(file = "metadata.txt", sep = "\t", header = T, row.names = 1, check.names = F)
row.names(samdf)
sample_data(ps_added_1) <- samdf
ps_added_1

# make abundance table
ps.glom <- tax_glom(physeq = ps_added_1,taxrank = 'Genus')

#perform DESeq2
diagdds = phyloseq_to_deseq2(ps_added_1, ~ Time_point)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")


res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps_added_1)[rownames(sigtab), ], "matrix"))
sigtab$Timepoint <- ifelse(sigtab$log2FoldChange>0,'Post_GFD','Baseline')
head(sigtab)

dim(sigtab)

#plot
library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))


#Plot_2
ggplot(sigtab, aes(y=Genus, x=log2FoldChange, color= Timepoint)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=5) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+ theme(axis.text.x = element_text(angle = 90, size = 15, colour = "black", vjust = 0.2, hjust = 0.5, face = "bold"), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold") ,strip.text = element_text(size = 18, colour = "black", face = "bold"), legend.title = element_text(size = 20, colour = "black", face = "bold"), legend.text = element_text(size = 22, face = "bold", colour = "black"), axis.text.y = element_text(size = 16, colour = "black", face = "bold"))

#Plot_3
ggplot(sigtab, aes(y=Phylum, x=log2FoldChange, color= Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_boxplot(size=1) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+ theme(axis.text.x = element_text(angle = 90, size = 15, colour = "black", vjust = 0.2, hjust = 0.5, face = "bold"), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold") ,strip.text = element_text(size = 18, colour = "black", face = "bold"), legend.title = element_text(size = 20, colour = "black", face = "bold"), legend.text = element_text(size = 22, face = "bold", colour = "black"), axis.text.y = element_text(size = 16, colour = "black", face = "bold"))

# Save results
write.csv(sigtab, "data.csv")
