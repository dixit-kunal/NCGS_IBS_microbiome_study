# The script can be used to run signifance testing for taxa abundance in three sample groups i.e. IBS patients, NCGS patients pre and post gluten free diet

# starting with a clean environment
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects
gc() #free up memory and report the memory usage.

# remove things that mask others
remove(list = conflicts(detail = TRUE)$.GlobalEnv)

# Dependencies
#install.packages("tidyverse", Ncpus = 8)
library(tidyverse)
library(data.table)
#install.packages("ggpubr", Ncpus = 8)
library(ggpubr)
#install.packages("rstatix", Ncpus = 8)
library(rstatix)
library(gtools)
library(readxl)

# assigning rename to use the function from dplyr
rename <- dplyr::rename

# Set working directory
setwd("path")

# Running in a loop per taxa #
#################################
# Using dataframe woth pathway with contributing species and abundance
# https://stackoverflow.com/questions/70386061/one-way-anova-analysis-by-row
kunal_long <- df %>%
  #  pivot_longer(cols = 2:ncol(.)) %>%
  pivot_wider(names_from = Species, values_from = value) %>%
  separate(group, into = c("groups"), sep = "\\d") %>%
  select(-variable)
kunal_long

# https://www.biostars.org/p/295214/
# group as factor
kunal_long$groups <- factor(kunal_long$groups, levels = c("ibs", "ncgs", "pg_ncgs"))
# replacing column names for ease
columns <- colnames(kunal_long)[2:ncol(kunal_long)] %>%
  as.data.frame() %>%
  rename("old_names" = ".") %>%
  mutate(new_names = paste0("taxa", 1:nrow(.)))

# replacing with new column names
colnames(kunal_long)[2:ncol(kunal_long)] <- columns$new_names
kunal_long$groups <- factor(kunal_long$groups, levels = c("ibs", "ncgs", "pg_ncgs"))

# since anova needs all three groups to have values, dropping columns that only have variables in two
sapply(lapply(kunal_long, unique), length)
length(sapply(lapply(kunal_long, unique), length)<2)
# dropping columns with less than two variables
todrop <- names(kunal_long)[which(sapply(lapply(kunal_long, unique), length)<4)]
df_clean = kunal_long[,!(names(kunal_long) %in% todrop)] %>% as.data.frame()
df_clean$groups <- kunal_long$groups
df_clean$groups <- as.factor(df_clean$groups)
df_clean <- df_clean %>% relocate(., "groups" = "groups")
df_clean_new <- df_clean %>% replace(is.na(.),0)

# store all formulae in a list
formulae <- lapply(colnames(df_clean_new)[2:ncol(df_clean_new)], function(x) as.formula(paste0(x, " ~ groups")))

# go through list and run aov()
res <- lapply(formulae, function(x) summary(aov(x, data = df_clean_new)))
names(res) <- format(formulae)
res

# extract just p-values
p <- unlist(lapply(res, function(x) x[[1]]$"Pr(>F)"[1]))
p

aov_pvalues <- data.frame(
  condition = sub(' ~ groups', '', names(p)),
  pvalue = p)

aov_merged_pathwaynames <- left_join(aov_pvalues, columns, by = c("condition" = "new_names"))
aov_pvalues_005 <- aov_merged_pathwaynames %>% filter(pvalue < 0.05)

aov_pvalues_005$adj.p <- p.adjust(aov_pvalues_005$pvalue, method = "fdr", n=length(aov_pvalues_005$pvalue))

write.csv(aov_pvalues_005, "D:\\Kunal\\metaphlan_&_humann3_analysis\\species\\anova_diff_taxa.csv")

# Extracting significant taxa
taxa <- df2 %>% filter(Species %in% aov_pvalues_005$old_names) %>% rename("count" = "value")
taxa

# Getting stats for significant taxa
library(ggpubr)
spp_stat <- compare_means(count ~ group, data = taxa, method = "wilcox.test", group.by = "Species") %>% 
  filter(p.signif != "ns")

write.csv(spp_stat, "D:\\Kunal\\metaphlan_&_humann3_analysis\\species\\annova_diff\\species_statistics.csv")


# plotting significant taxa
taxa %>% 
  ggplot(aes(x = group, y = log10(count))) +
  geom_boxplot(aes(fill = group)) +
  #geom_point(aes(color = group)) +
  geom_jitter(aes()) +
  facet_wrap(.~Species) + 
  scale_fill_manual(values = c("#3399FF", "orange", "#996600")) + 
  scale_color_manual(values = c("#3399FF", "orange", "#996600")) + theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black", vjust = 0.2, hjust = 0.5, face = "bold"), axis.title.x = element_text(size = 20, face = "bold"), axis.title.y = element_text(size = 20, face = "bold"),strip.text = element_text(size = 18, colour = "black", face = "bold"), legend.title = element_text(size = 20, colour = "black", face = "bold"), legend.text = element_text(size = 22, face = "bold", colour = "black"), axis.text.y = element_text(size = 16, colour = "black", face = "bold"))+ scale_y_continuous()

#+geom_signif(
    comparisons = list(c("ibs", "ncgs"), c("ibs", "pg_ncgs"), c("ncgs", "pg_ncgs")),
    map_signif_level = TRUE
    #label.y = 1
  )

             
