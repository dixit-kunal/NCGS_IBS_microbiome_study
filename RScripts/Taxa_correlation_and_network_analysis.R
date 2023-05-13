# The script can be used to calculate taxa to taxa correlation values based on Spearman's correlation. This was further used to generate microbial co-occurrence networks and network statistics to illustrate signature taxa in our study. 

# Set working directory
setwd('Path')

# Load necessary6y packages
library(phyloseq)
library(ggplot2)
library(Hmisc)

# Load phyloseq object
ps <- readRDS(file = 'ps.rds')
ps


temp.taxa <- as.data.frame(tax_table(ps))
temp.taxa$Genus <- gsub(pattern = '_.*',replacement = '',x = temp.taxa$Genus,perl = T)
tax_table(ps) <- tax_table(as.matrix(temp.taxa))

# Glom to genus level
ps.glom <- tax_glom(ps, taxrank="Genus")

asv.table <- as.data.frame(otu_table(ps.glom))
asv.table <- t(asv.table)

# Extract sequences
asv.names <- paste0('ASV',1:dim(asv.table)[1])
asv.headers <- paste0('>',asv.names,'\n')
asv.seq <- rownames(asv.table)
write(x = paste0(asv.headers,asv.seq),file = 'asv_seqs.fna')

# otu table
rownames(asv.table) <- asv.names
write.table(x = asv.table,file = 'asv_table.txt',quote = F,sep = '\t',col.names = NA)

# taxonomy
tax.tab <- as.data.frame(tax_table(ps.glom))
rownames(tax.tab) <- asv.names
write.table(x = tax.tab,file = 'tax_table.txt',quote = F,sep = '\t',col.names = NA)

# Sample data
sample.data <- as.data.frame(sample_data(ps.glom))
write.table(x = sample.data,file = 'metadata_table.txt',quote = F,sep = '\t',col.names = NA)


#Correlation analysis
corr <- read.csv("path\\corr.csv", row.names=1) # the file was curated for format manually and reloaded in R environment
corr=t(corr)
flattenCorrMatrix <- function(cormat, pmat) {
     ut <- upper.tri(cormat)
     data.frame(
          row = rownames(cormat)[row(cormat)[ut]],
          column = rownames(cormat)[col(cormat)[ut]],
          cor  =(cormat)[ut],
          p = pmat[ut]
     )
}

res2 <- rcorr(as.matrix(corr[,1:717])) # as the number of interactions was 717
correlation=flattenCorrMatrix(res2$r, res2$P)

# Download correlation file
write.csv(correlation, "correlation.csv")

# Further, the significant positive interactions from the file "correlation.csv" was used ahead to generate bacterial networks






##################
# Kunal Networks #
##################
setwd('path')

# Dependencies
library(tidyverse)
library(igraph)
library(viridis)

# Data
el <- read.csv("correlation_ibs_stool.csv")
el

grph <- graph.data.frame(el, directed = FALSE)
grph

# Adding correlation as edge weights
E(grph)$weight <- el$cor
E(grph)$weight

V(grph)
vertices <- length(V(grph))
vertices
# 41

E(grph)
edges <- length(E(grph))
edges
# 75

density <- edges/((vertices*(vertices-1))/2)
density
# [1]  0.09146341

# Plotting the graph
V(grph)$size <- (degree(grph) + 1)
V(grph)$size <- V(grph)$size
# V(grph)$size<-V(grph)$size/4
V(grph)$name


#V(grph)$label.cex <- degree(grph)/5
plot(grph,
     vertex.size=V(grph)$size,
     #     vertex.label=NA,
     edge.width=E(grph)$weight,
     layout=layout_on_sphere(grph),
     edge.curved=0.5)

l <- layout.fruchterman.reingold(grph, niter=5000, area=vcount(grph)^4*10)

plot(grph, layout=l, 
     #     edge.arrow.size=0.5, 
     vertex.label.cex=0.9, 
     vertex.label.family="Helvetica",
     vertex.label.font=2,
     vertex.shape="circle", 
     vertex.size=V(grph)$size, 
     vertex.label.color="black", 
     edge.width=E(grph)$weight,
     edge.curved=0.3)


# Edge Color by Weight
#Custom colorblind pallette, see: https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
customvermillion<-rgb(213/255,94/255,0/255)
custombluegreen<-rgb(0/255,158/255,115/255)
customblue<-rgb(0/255,114/255,178/255)
customskyblue<-rgb(86/255,180/255,233/255)
customreddishpurple<-rgb(204/255,121/255,167/255)
E(grph)$color <- custombluegreen
E(grph)$color[E(grph)$weight<0] <- customreddishpurple
E(grph)$width <- abs(E(grph)$weight)*2
plot(grph,
     vertex.label=NA,
     layout=layout.circle(grph))

# Edge Weights
plot(density((E(grph)$weight)),xlab="Edge Weight",main="")
# pdf(file = "VALSOREY_second_interactions_boxplot.pdf", width = 8, height = 8)
boxplot(abs(E(grph)$weight)~(E(grph)$weight>0),
        xlab="Positive Interaction?",
        ylab="Strength of Interaction")

#Remove edges with very low weight
hist(E(grph)$weight)
weight_threshold <- 0
grph <- delete.edges(grph,which(abs(E(grph)$weight)<weight_threshold))

#Remove negative edges
grph.pos <- delete.edges(grph,which(E(grph)$weight<0))
plot(grph.pos,
     vertex.label=NA,
     edge.color="black",
     layout=layout_with_fr(grph.pos))


#Remove unconnected vertices
grph.pos <- delete.vertices(grph.pos, which(degree(grph.pos)<1))
set.seed(42)
plot(grph.pos,
     vertex.label=NA,
     edge.color="black",
     layout=layout_with_fr(grph.pos),
     edge.curved = 0.5)

# degree distribution
dd.grph.pos <- degree.distribution(grph.pos)
plot(0:(length(dd.grph.pos)-1), dd.grph.pos, type='b',
     ylab="Frequency", xlab="Degree", main="Degree Distributions")


# Degree is the number of edges that a node (vertex) has. We see that most nodes are connected to few other nodes when consideringly only positive interactions above our weight threshold.
grph.pos_deg<-degree(grph.pos, v=V(grph.pos), mode="all")
fine = 500 # this will adjust the resolving power.

#this gives you the colors you want for every point
graphCol = viridis(fine)[as.numeric(cut(grph.pos_deg,breaks = fine))]

# now plot
set.seed(42)
plot(grph.pos, vertex.color=graphCol,
     edge.color="black",
     vertex.label=NA,
     layout=layout_on_sphere(grph.pos),
     edge.curved = 0.5)

#    layout=layout_with_fr(grph.pos),
# vertex.label=if(1.5^(V(grph.pos)$size)>3){V(grph.pos)$name} else {
#   V(grph.pos)$name==""
#   })

grph.pos_bw<-betweenness(grph.pos, directed=F)

#this gives you the colors you want for every point
graphCol = viridis(fine)[as.numeric(cut(grph.pos_bw,breaks = fine))]

# now plot
set.seed(42)
plot(grph.pos, vertex.color=graphCol,
     vertex.label=NA,
     edge.color="black",
     layout=layout_with_fr(grph.pos), edge.curved = 0.5)


## Clustering
grph.pos.greedy <- cluster_fast_greedy(grph.pos, weights=E(grph.pos)$weight)
modularity(grph.pos.greedy)
sizes(grph.pos.greedy)

colourCount = length(unique(grph.pos.greedy$membership)) # this will adjust the resolving power.

cluster_col = rainbow(colourCount)[as.numeric(cut(grph.pos.greedy$membership,breaks = colourCount))]

# now plot
set.seed(42)
plot(grph.pos, vertex.color=cluster_col,
     vertex.label=NA,
     edge.color="black",
     layout=layout_on_sphere(grph.pos),
     edge.curved = 0.5)

grph.pos.louvain <- cluster_louvain(grph.pos, weights=E(grph.pos)$weight)
modularity(grph.pos.louvain)
sizes(grph.pos.louvain)

colourCount = length(unique(grph.pos.louvain$membership)) # this will adjust the resolving power.

cluster_col = rainbow(colourCount)[as.numeric(cut(grph.pos.louvain$membership,breaks = colourCount))]

# now plot
set.seed(42)
plot(grph.pos, vertex.color=cluster_col,
     vertex.label= NA,
     edge.color="Black",
     vertex.label.cex=1.1,
     vertex.label.font=2,
     edge.width=1.5,
     layout=layout_with_fr(grph.pos),
     vertex.label.color="black",
     edge.curved = 0.5)

# circular plot
set.seed(42)
plot(grph.pos, vertex.color=cluster_col,
     #vertex.label=NA,
     edge.color="Black",
     vertex.label.cex=1.1,
     vertex.label.font=2,
     edge.width=1.5,
     layout=layout.circle(grph.pos),
     vertex.label.color="black",
     edge.curved = 0)


# Getting cluster memberships
V(grph.pos)$cluster=grph.pos.louvain$membership
vertex_attr(grph.pos, index = V(grph.pos))

ids <- which(sizes(grph.pos.louvain)<=5)
grph.pos.main.communities <- delete_vertices(grph.pos,which(V(grph.pos)$cluster %in% ids))

nodes <- V(grph.pos.main.communities)$name
nodes

hist(degree(grph.pos.main.communities))

# now plot
set.seed(42)
plot(grph.pos, vertex.color=cluster_col,
     vertex.label= ifelse(V(grph.pos)$name %in% nodes, V(grph.pos)$name, ""),
     edge.color="Black",
     vertex.label.cex=1.1,
     vertex.label.font=2,
     edge.width=1.5,
     layout=layout_with_fr(grph.pos),
     vertex.label.color="black",
     edge.curved = 0.5)


list_ncgs_duo <- c("Bacillus", "Alcaligenes", "Leucobacter", "Providencia", "C39", "Aeromicrobium")

list_ncgs_stool <- c("Succinivibrio", "Lysinibacillus", "Roseburia CAG-352", "Desulfovibrio", "[Eubacterium] ruminantium group", "Gemella", "Anaerostipes", "Staphylococcus", "Lachnospiraceae AC2044 group", "Methanobrevibacter", "Senegalimassilia", "Intestinimonas", "Weissella", "Rothia", "[Eubacterium] siraeum group", "Caproiciproducens", "Pseudobutyrivibrio", "GCA-900066575", "Mogibacterium", "Treponema", "Lachnospiraceae UCG-001", "Victivallis", "UBA1819", "Candidatus Soleaferrea", "Catenisphaera", "Barnesiella", "Elusimicrobium"
)

list_ncgs_sig <- c("Bacillus", "Lysinibacillus", "Acinetobacter", "Rhodococcus", "Achromobacter", "Ochrobactrum", "Desulfovibrio", "Gordonia", "Delftia", "Agromyces", "CAG-56", "Methanobrevibacter", "Brevibacterium", "Caproiciproducens", "Afipia", "Pseudobutyrivibrio", "Deinococcus", "Herbaspirillum", "Pelomonas", "Vagococcus", "Chryseobacterium", "Peptococcus", "Curvibacter", "UCG-009", "Thermomonas", "Reyranella")

list_ibs_duo <- c("Prevotella", "Agathobacter", "Lactobacillus", "Faecalibacterium", "Shewanella", "Bacteroides", "Roseburia", "Dialister", "Bifidobacterium", "Catenibacterium", "Blautia", "Dorea", "Subdoligranulum", "[Ruminococcus] torques group", "Megamonas", "Sutterella", "Collinsella", "Lachnospiraceae NK4A136 group", "Clostridium sensu stricto 1", "Ralstonia", "Lachnospira", "Zhihengliuella", "Lachnospiraceae UCG-004", "Actinobacillus", "UCG-002", "Brevibacillus", "Romboutsia", "Fusicatenibacter", "Holdemanella", "Lachnospiraceae UCG-010", "Butyricicoccus", "[Eubacterium] ruminantium group", "Lachnoclostridium", "Ruminococcus", "Anaerostipes", "Lachnospiraceae ND3007 group", "Phascolarctobacterium", "Alloprevotella", "Granulicatella", "[Eubacterium] hallii group", "Prevotellaceae UCG-003", "Pseudoxanthomonas", "Intestinibacter", "Howardella", "[Eubacterium] xylanophilum group", "Allisonella", "Colidextribacter", "Odoribacter", "CAG-56", "Acidovorax", "Erysipelotrichaceae UCG-003", "Actinomyces", "Coprococcus", "Stenotrophomonas", "Senegalimassilia", "Parabacteroides", "Sarcina", "Lachnospiraceae FCS020 group", "Weissella", "Christensenellaceae R-7 group", "UCG-003", "Mannheimia", "Stomatobaculum", "Solobacterium", "[Ruminococcus] gauvreauii group", "Terrisporobacter", "Bilophila", "[Eubacterium] siraeum group", "Fournierella", "Afipia", "Chitinophaga", "UCG-005", "Mesorhizobium", "Caulobacter", "Alishewanella", "Oribacterium", "[Eubacterium] nodatum group", "Incertae Sedis", "Family XIII AD3011 group", "Slackia", "Enterococcus", "Flavonifractor", "Alistipes", "Turicibacter", "Sphingobium", "Aquabacterium", "F0058", "Peptostreptococcus", "Marvinbryantia", "Variovorax", "Mogibacterium", "Cupriavidus", "Tepidimonas", "Oxalobacter", "Johnsonella", "Psychrobacillus", "Eikenella", "Psychrobacter", "Empedobacter", "Knoellia", "Catonella", "Treponema", "Butyrivibrio", "Adlercreutzia", "Skermanella", "Lachnospiraceae UCG-001", "Nocardioides", "Zoogloea", "Gordonibacter", "Williamsia", "Conservatibacter", "[Renibacterium] salmoninarum group", "[Eubacterium] oxidoreducens group", "SN8", "Lachnospiraceae UCG-008", "Lachnospiraceae NC2004 group")

list_ibs_sig <- c("Psychrobacter", "Veillonella", "Zhihengliuella", "Paraprevotella", "Flavonifractor", "Turicibacter", "[Ruminococcus] gnavus group", "Lautropia", "UCG-008", "Centipeda", "Dubosiella", "Gemmobacter", "Ruminobacter", "Lachnospiraceae UCG-008", "SN8", "Azotobacter", "[Renibacterium] salmoninarum group", "Pseudochrobactrum")

list_ibs_stool = c("Megasphaera", "Aeromonas", "Fusobacterium", "Megamonas", "Aggregatibacter", "[Ruminococcus] gnavus group", "Pseudocitrobacter")

list_ncgs_pg_duo <- c("Stenotrophomonas", "Lysinibacillus", "Shewanella", "Solibacillus", "Fusobacterium", "Ralstonia", "Actinobacillus", "Porphyromonas", "Pseudoxanthomonas", "Delftia", "Brachybacterium", , "Paenochrobactrum", "Pseudochrobactrum", "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "Agromyces", "Acidovorax", "Leptotrichia", "Capnocytophaga", "Mycoplasma", "Afipia", "Hydrogenophaga", "Mesorhizobium", "Caulobacter", "Leucobacter", "Deinococcus", "Herbaspirillum", "Lachnoanaerobaculum", "Sphingobium", "Chryseobacterium", "F0058", "Lentimicrobium", "Variovorax", "Mogibacterium", "Streptobacillus", "Lautropia", "Psychroglaciecola", "Tepidimonas", "Curvibacter", "Azotobacter", "Bergeyella", "Cloacibacterium", "Filifactor", "Centipeda", "Cellulosimicrobium", "Ciceribacter", "Tannerella", "Butyrivibrio", "Phenylobacterium")

list_ncgs_pg_sig <- c("Campylobacter", "Pseudarthrobacter", "Pseudoxanthomonas", "Dechloromonas", "Pseudochrobactrum", "Hydrogenophaga", "Sphingobium", "Porphyromonas", "Caulobacter", "Methylotenera", "Microbacterium", "UCG-008", "Phreatobacter", "Pseudaminobacter", "Massilia", "Ottowia", "Undibacterium", "Dietzia", "Rhodobacter", "Phenylobacterium")


# now plot
set.seed(42)
plot(grph.pos, vertex.color=cluster_col,
     vertex.label= ifelse(V(grph.pos)$name %in% list_ibs_stool, V(grph.pos)$name, ""),
     edge.color="Black",
     vertex.label.cex=1.4,
     vertex.label.font=2,
     edge.width=1.5,
     layout=layout_with_fr(grph.pos),
     vertex.label.color="black",
     edge.curved = 0.5)



cluster_id <- V(grph.pos.main.communities)$cluster
nodes.ibs.duo<-as.data.frame(cbind(nodes, cluster_id))
nodes.ibs.duo

######## Network stats

# Number of nodes
grph$nodes <- length(V(grph))

# Number of edges
grph$edges <- length(E(grph))

# Density of the network
# Density == Actual connections/Potential connections
# Potential connections == (n*(n-1))/2
grph$density <- length(E(grph))/((length(V(grph))*(length(V(grph))-1))/2)

# Global clustering coefficient
grph$transitivity <- transitivity(grph)

# Average clustering coefficient
grph$avg_transitivity <- transitivity(grph, type = "average")

# Modularity
grph$modularity <- modularity(cluster_fast_greedy(grph.pos, weights=E(grph.pos)$weight))

# Diameter
grph$diameter <- diameter(grph.pos)

# Assortativity
grph$assortativity <- assortativity_degree(grph)

# Mean distance
grph$mean_distance <- mean_distance(grph.pos)

# Cliques
grph$cliques <- length(max_cliques(grph))

# Cluster sizes
grph$clusters <- length(sizes(cluster_fast_greedy(grph.pos)))



ibs.sig.0.05.stats <- data.frame(Stat = c("nodes", "edges", "density", "transitivity", "avg_transitivity", "modularity", "diameter", "assortativity", "mean_distance", "cliques", "clusters"),
                            Value = c(grph$nodes, grph$edges, grph$density, grph$transitivity, grph$avg_transitivity, grph$modularity, grph$diameter, grph$assortativity, grph$mean_distance, grph$cliques, grph$clusters))



######### join stat tables & save ###########

table <- full_join(x = healthy.stats, y = ibs.stats, by = "Stat", all = TRUE) %>% 
  full_join(x = ., y = ncgs.baseline.stats, by = "Stat", all = TRUE) %>% full_join(x = ., y = ncgs.pg.stats, by = "Stat", all = TRUE)

colnames(table) <- c("Stats","Healthy","IBS","NCGS","NCGS_PG")
write.csv(table, "stats_all_positive.csv")
