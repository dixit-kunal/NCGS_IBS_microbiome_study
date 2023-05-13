# The script can be used to run signifance testing for pathway testing in three sample groups i.e. IBS patients, NCGS patients pre and post gluten free diet

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

# Running in a loop per pathway #
#################################
# set working directory
setwd("path")

# Using dataframe woth pathway with contributing species and abundance 
# https://stackoverflow.com/questions/70386061/one-way-anova-analysis-by-row
kunal_long <- df %>%
  #  pivot_longer(cols = 2:ncol(.)) %>%
  pivot_wider(names_from = pathway, values_from = count) %>%
  separate(group, into = c("groups"), sep = "\\d") %>%
  select(-sample)
kunal_long

# https://www.biostars.org/p/295214/
# group as factor
kunal_long$groups <- factor(kunal_long$groups, levels = c("IBS", "NCGS", "NCGS_PG"))
# replacing column names for ease
columns <- colnames(kunal_long)[2:ncol(kunal_long)] %>%
  as.data.frame() %>%
  rename("old_names" = ".") %>%
  mutate(new_names = paste0("path", 1:nrow(.)))

# replacing with new column names
colnames(kunal_long)[2:ncol(kunal_long)] <- columns$new_names
kunal_long$groups <- factor(kunal_long$groups, levels = c("IBS", "NCGS", "NCGS_PG"))

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

#write.csv(aov_pvalues_005, "anova_diff_species_pathabund.csv")

# Extracting significant path
path <- df %>% filter(pathway %in% aov_pvalues_005$old_names)
path

# Getting stats for significant path
library(ggpubr)
spp_stat <- compare_means(count ~ group, data = path, method = "wilcox.test", group.by = "pathway") %>% 
  filter(p.signif != "ns")

#write.csv(spp_stat, "species_pathabund_statistics.csv")

# plotting significant path
path %>% 
  #filter(pathway == "PWY-7013: (S)-propane-1,2-diol degradation") %>% 
  ggplot(aes(x = group, y = count)) +
  geom_boxplot(aes(fill = group)) +
  #geom_point(aes(color = group)) +
  geom_jitter(aes()) +
  facet_wrap(.~pathway) + 
  scale_fill_manual(values = c("#3399FF", "orange", "#996600")) + 
  scale_color_manual(values = c("#3399FF", "orange", "#996600")) + theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black", vjust = 0.2, hjust = 0.5, face = "bold"), axis.title.x = element_text(size = 20, face = "bold"), axis.title.y = element_text(size = 20, face = "bold"),strip.text = element_text(size = 18, colour = "black", face = "bold"), legend.title = element_text(size = 20, colour = "black", face = "bold"), legend.text = element_text(size = 22, face = "bold", colour = "black"), axis.text.y = element_text(size = 16, colour = "black", face = "bold"))+ scale_y_continuous()


# setting up unique pathnames
aov_pvalues_005_new <- aov_pvalues_005 %>% 
  mutate(path_names = old_names) %>% 
  separate(path_names, c("path_names_v2", "taxa"), sep = "\\|g__") %>% 
  mutate(path_names_v2 = as.factor(path_names_v2))

unique(levels(aov_pvalues_005_new$path_names_v2))
