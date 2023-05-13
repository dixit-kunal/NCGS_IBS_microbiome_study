# The script can be used to plot ordination plots (NMDS) for core microbiome and test for significance using PERMANOVA. Here we use prevelence of 0.67 to define core microbiome

# Load necessary packages
library(DESeq2)
library(phyloseq)
library("ggplot2"); packageVersion("ggplot2")
library(reshape2)
library(ggpubr)
library(RColorBrewer)
library(ggpubr)
library(ape)
library(vegan)
library(rlang)
library(microbiome)

# set working directory
setwd("path")

# Load phyloseq object
ps <- readRDS("ps.rds")

# Load metadata
samdf <- read.table(file = "merged.txt", sep = "\t", header = T, row.names = 1, check.names = F)

row.names(samdf)
sample_data(ps) <- samdf
ps
################################################ IBS vs NCGS ##################################################################

ps.ibs.stool <- subset_samples(physeq = ps,Sample_Description == 'IBS_stool')
ps.ncgs.stool <- subset_samples(physeq = ps, Sample_Description == 'NCGS_stool')

ps.ibs.duo <- subset_samples(physeq = ps, Sample_Description == 'IBS_small_intestine')
ps.ncgs.duo <- subset_samples(physeq = ps, Sample_Description == 'NCGS_small_intestine')

ps.ibs.sig <- subset_samples(physeq = ps, Sample_Description == 'IBS_large_intestine')
ps.ncgs.sig <- subset_samples(physeq = ps, Sample_Description == 'NCGS_large_intestine')

############## core microbiota ##############
############### Stool ################
ps.ibs.stool.core <- core(ps.ibs.stool, detection = 0.0000, prevalence = .67)
ps.ncgs.stool.core <- core(ps.ncgs.stool, detection = 0.0000, prevalence = .67)
merged_stool_phyloseq <- merge_phyloseq(ps.ibs.stool.core,ps.ncgs.stool.core)

ps.prop <- transform_sample_counts(merged_stool_phyloseq, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="Sample_Description", title="IBS - NCGS Stool Core Microbiome") + geom_point(size=4, alpha=2.5, shape= 15)+ stat_ellipse(geom = "polygon", linetype = 1,alpha=0) + scale_color_manual(values = c("#3399FF", "orange"))+ theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black", vjust = 0.2, hjust = 0.5, face = "bold"), axis.title.y = element_text(size = 15, face = "bold"), axis.title.x = element_text(size = 15, face = "bold"), legend.title = element_text(size = 20, colour = "black", face = "bold"), legend.text = element_text(size = 16, face = "bold", colour = "black"), axis.text.y = element_text(size = 16, colour = "black", face = "bold"))+ scale_y_continuous() + theme_bw()


metadata.stool <- as(sample_data(merged_stool_phyloseq), "data.frame")
asv.tab.stool <- as.data.frame(otu_table(merged_stool_phyloseq))
permanova_stool <- adonis2(phyloseq::distance(merged_stool_phyloseq, method="bray") ~ Sample_Description,
                           data = metadata.stool)



############### Duo ################

ps.ibs.duo.core <- core(ps.ibs.duo, detection = 0.0000, prevalence = .67)
ps.ncgs.duo.core <- core(ps.ncgs.duo, detection = 0.0000, prevalence = .67)
merged_duo_phyloseq <- merge_phyloseq(ps.ibs.duo.core,ps.ncgs.duo.core)

ps.prop <- transform_sample_counts(merged_duo_phyloseq, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="Sample_Description", title="IBS - NCGS Small intesine Core Microbiome") + geom_point(size=4, alpha=2.5, shape= 19)+ stat_ellipse(geom = "polygon", linetype = 1,alpha=0) + scale_color_manual(values = c("#3399FF", "orange"))+ theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black", vjust = 0.2, hjust = 0.5, face = "bold"), axis.title.y = element_text(size = 15, face = "bold"), axis.title.x = element_text(size = 15, face = "bold"), legend.title = element_text(size = 20, colour = "black", face = "bold"), legend.text = element_text(size = 16, face = "bold", colour = "black"), axis.text.y = element_text(size = 16, colour = "black", face = "bold"))+ scale_y_continuous() + theme_bw()


metadata.duo <- as(sample_data(merged_duo_phyloseq), "data.frame")
asv.tab.duo <- as.data.frame(otu_table(merged_duo_phyloseq))
permanova_duo<- adonis2(phyloseq::distance(merged_duo_phyloseq, method="bray") ~ Sample_Description,
                        data = metadata.duo)



############### Sig ################

ps.ibs.sig.core <- core(ps.ibs.sig, detection = 0.0000, prevalence = .67)
ps.ncgs.sig.core <- core(ps.ncgs.sig, detection = 0.0000, prevalence = .67)
merged_sig_phyloseq <- merge_phyloseq(ps.ibs.sig.core,ps.ncgs.sig.core)


ps.prop <- transform_sample_counts(merged_sig_phyloseq, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="Sample_Description", title="IBS - NCGS Large intesine Core Microbiome") + geom_point(size=4, alpha=2.5, shape= 17)+ stat_ellipse(geom = "polygon", linetype = 1,alpha=0) + scale_color_manual(values = c("#3399FF", "orange"))+ theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black", vjust = 0.2, hjust = 0.5, face = "bold"), axis.title.y = element_text(size = 15, face = "bold"), axis.title.x = element_text(size = 15, face = "bold"), legend.title = element_text(size = 20, colour = "black", face = "bold"), legend.text = element_text(size = 16, face = "bold", colour = "black"), axis.text.y = element_text(size = 16, colour = "black", face = "bold"))+ scale_y_continuous() + theme_bw()


metadata.sig <- as(sample_data(merged_sig_phyloseq), "data.frame")
asv.tab.sig <- as.data.frame(otu_table(merged_sig_phyloseq))
permanova_sig <- adonis2(phyloseq::distance(merged_sig_phyloseq, method="bray") ~ Sample_Description,
                         data = metadata.sig)

################## Save permanova ##########################

setwd("D:\\Kunal\\metaphlan_&_humann3_analysis\\core_ordination")

write.csv(permanova_stool, "permanova.ibs.ncgs.stool.csv")
write.csv(permanova_duo, "permanova.ibs.ncgs.duo.csv")
write.csv(permanova_sig, "permanova.ibs.ncgs.sig.csv")




################################## NCGS vs NCGS_PG #################################################

ps.ncgs <- subset_samples(physeq = ps, Samples == 'NCGS')

ps.ncgs.ncgspg.stool <- subset_samples(physeq = ps.ncgs, Sample_Type == 'stool')
#ps.ncgs.stool.genus <- tax_glom(physeq = ps.ncgs.stool, taxrank = "Genus")
#metadata.ncgs.stool.genus <- sample_data(ps.ncgs.stool.genus)


ps.ncgs.ncgspg.duo <- subset_samples(physeq = ps.ncgs, Sample_Type == 'Duodenal')
#ps.ncgs.duo.genus <- tax_glom(physeq = ps.ncgs.duo, taxrank = "Genus")
#metadata.ncgs.duo.genus <- sample_data(ps.ncgs.duo.genus)
#otu.ncgs.ibs.duo.genus <- as_data_frame(otu_table(ps.ncgs.duo.genus))


ps.ncgs.ncgspg.sig <- subset_samples(physeq = ps.ncgs, Sample_Type == 'Sigmoid')
#ps.ncgs.sig.genus <- tax_glom(physeq = ps.ncgs.sig, taxrank = "Genus")
#metadata.ncgs.sig.genus <- sample_data(ps.ncgs.sig.genus)


#########################
ps.ncgspg.st <- subset_samples(physeq = ps.ncgs.ncgspg.stool,Time_point == 'Post_GFD')
ps.ncgs.st <- subset_samples(physeq = ps.ncgs.ncgspg.stool, Time_point == 'Baseline')

ps.ncgspg.duodenal <- subset_samples(physeq = ps.ncgs.ncgspg.duo, Time_point == 'Post_GFD')
ps.ncgs.duodenal <- subset_samples(physeq = ps.ncgs.ncgspg.duo, Time_point == 'Baseline')

ps.ncgspg.sigmoid <- subset_samples(physeq = ps.ncgs.ncgspg.sig, Time_point == 'Post_GFD')
ps.ncgs.sigmoid <- subset_samples(physeq = ps.ncgs.ncgspg.sig, Time_point == 'Baseline')

############## core microbiota ##############

############### Stool ################
ps.ncgspg.st.core <- core(ps.ncgspg.st, detection = 0.0000, prevalence = .67)
ps.ncgs.st.core <- core(ps.ncgs.st, detection = 0.0000, prevalence = .67)
merged_ncgs_ncgspg_st_phyloseq <- merge_phyloseq(ps.ncgspg.st.core,ps.ncgs.st.core)

ps.prop <- transform_sample_counts(merged_ncgs_ncgspg_st_phyloseq, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="Time_point", title="NCGS - NCGS_PG Stool Core Microbiome") + geom_point(size=4, alpha=2.5, shape= 15)+ stat_ellipse(geom = "polygon", linetype = 1,alpha=0) + scale_color_manual(values = c("orange", "#CC6600", "#3399FF"))+ theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black", vjust = 0.2, hjust = 0.5, face = "bold"), axis.title.y = element_text(size = 15, face = "bold"), axis.title.x = element_text(size = 15, face = "bold"), legend.title = element_text(size = 20, colour = "black", face = "bold"), legend.text = element_text(size = 16, face = "bold", colour = "black"), axis.text.y = element_text(size = 16, colour = "black", face = "bold"))+ scale_y_continuous() + theme_bw()


metadata.ncgs.ncgspg.stool <- as(sample_data(merged_ncgs_ncgspg_st_phyloseq), "data.frame")
asv.tab.stool <- as.data.frame(otu_table(merged_ncgs_ncgspg_st_phyloseq))
permanova_ncgs_ncgspg_stool <- adonis2(phyloseq::distance(merged_ncgs_ncgspg_st_phyloseq, method="bray") ~ Time_point,
                           data = metadata.ncgs.ncgspg.stool)


############### Duodenal ################
ps.ncgspg.duodenal.core <- core(ps.ncgspg.duodenal, detection = 0.0000, prevalence = .67)
ps.ncgs.duodenal.core <- core(ps.ncgs.duodenal, detection = 0.0000, prevalence = .67)
merged_ncgs_ncgspg_duodenal_phyloseq <- merge_phyloseq(ps.ncgspg.duodenal.core,ps.ncgs.duodenal.core)

ps.prop <- transform_sample_counts(merged_ncgs_ncgspg_duodenal_phyloseq, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="Time_point", title="NCGS - NCGS_PG Small intesine Core Microbiome") + geom_point(size=4, alpha=2.5, shape= 19)+ stat_ellipse(geom = "polygon", linetype = 1,alpha=0) + scale_color_manual(values = c("orange", "#CC6600", "#3399FF"))+ theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black", vjust = 0.2, hjust = 0.5, face = "bold"), axis.title.y = element_text(size = 15, face = "bold"), axis.title.x = element_text(size = 15, face = "bold"), legend.title = element_text(size = 20, colour = "black", face = "bold"), legend.text = element_text(size = 16, face = "bold", colour = "black"), axis.text.y = element_text(size = 16, colour = "black", face = "bold"))+ scale_y_continuous() + theme_bw()


metadata.ncgs.ncgspg.duodenal <- as(sample_data(merged_ncgs_ncgspg_duodenal_phyloseq), "data.frame")
asv.tab.duodenal <- as.data.frame(otu_table(merged_ncgs_ncgspg_duodenal_phyloseq))
permanova_ncgs_ncgspg_duodenal <- adonis2(phyloseq::distance(merged_ncgs_ncgspg_duodenal_phyloseq, method="bray") ~ Time_point, data = metadata.ncgs.ncgspg.duodenal)


############### Sigmoid ################
ps.ncgspg.sigmoid.core <- core(ps.ncgspg.sigmoid, detection = 0.0000, prevalence = .67)
ps.ncgs.sigmoid.core <- core(ps.ncgs.sigmoid, detection = 0.0000, prevalence = .67)
merged_ncgs_ncgspg_sigmoid_phyloseq <- merge_phyloseq(ps.ncgspg.sigmoid.core,ps.ncgs.sigmoid.core)

ps.prop <- transform_sample_counts(merged_ncgs_ncgspg_sigmoid_phyloseq, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="Time_point", title="NCGS - NCGS_PG Large intesine Core Microbiome") + geom_point(size=4, alpha=2.5, shape= 17)+ stat_ellipse(geom = "polygon", linetype = 1,alpha=0) + scale_color_manual(values = c("orange", "#CC6600", "#3399FF"))+ theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black", vjust = 0.2, hjust = 0.5, face = "bold"), axis.title.y = element_text(size = 15, face = "bold"), axis.title.x = element_text(size = 15, face = "bold"), legend.title = element_text(size = 20, colour = "black", face = "bold"), legend.text = element_text(size = 16, face = "bold", colour = "black"), axis.text.y = element_text(size = 16, colour = "black", face = "bold"))+ scale_y_continuous() + theme_bw()


metadata.ncgs.ncgspg.sigmoid <- as(sample_data(merged_ncgs_ncgspg_sigmoid_phyloseq), "data.frame")
asv.tab.sigmoid <- as.data.frame(otu_table(merged_ncgs_ncgspg_sigmoid_phyloseq))
permanova_ncgs_ncgspg_sigmoid <- adonis2(phyloseq::distance(merged_ncgs_ncgspg_sigmoid_phyloseq, method="bray") ~ Time_point, data = metadata.ncgs.ncgspg.sigmoid)


################## Save permanova ##########################

setwd("D:\\Kunal\\metaphlan_&_humann3_analysis\\core_ordination")

write.csv(permanova_ncgs_ncgspg_stool, "permanova.ncgs.ncgspg.stool.csv")
write.csv(permanova_ncgs_ncgspg_duodenal, "permanova.ncgs.ncgspg.duo.csv")
write.csv(permanova_ncgs_ncgspg_sigmoid, "permanova.ncgs.ncgspg.sig.csv")
