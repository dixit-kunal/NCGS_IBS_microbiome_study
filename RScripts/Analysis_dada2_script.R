# DADA2 Pipeline for analyzing 16S rRNA gene metabarcoding data
#This script includes Quality profile inspection, Filter and trim, Learning the Error Rates, Merging paired end reads, Construct sequence table, Remove chimeras, Assign taxonomy, Generate Phyloseq object

# Install DADA2 using devtools

install.packages("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2", ref="v1.16") # change the ref argument to get other versions

# Open file path
library(dada2); packageVersion("dada2")
path <- "/home/mcc/kunal_dixit/NCGS/sequencing"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.namesR <- sapply(strsplit(basename(fnRs), "_"), `[`, 1)
sample.names

# Plot Quality Profiles
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])


# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,210), trimLeft = c(19,20), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE; Primer removal using trimleft; Adjust truncLen as per data quality

head(out)

write.csv(out, "out.csv")

# Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Plotting
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Sample inference and merger of paired-end reads
save.image("healthy.RData")
mergers <- vector("list", length(sample.names))
denoisedF <- vector()
denoisedR <- vector()
numMerged <- vector()
getN <- function(x) sum(getUniques(x))
names(mergers) <- sample.names

# dereplication and merging
print("dereplication and merging")
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR,maxMismatch=2,minOverlap=10)
  mergers[[sam]] <- merger
  #count number of processed reads
  denoisedF <- c(denoisedF,getN(ddF))
  denoisedR <- c(denoisedR,getN(ddR))
  numMerged <- c(numMerged,getN(merger))
}
rm(derepF); rm(derepR)

track <- cbind(out,denoisedF,denoisedR,numMerged)

# Construct Sequence Table
print("making sequence table")
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "seqtab.rds")
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab <- seqtab[,nchar(colnames(seqtab)) %in% 252:254]

print("Remove chimera")
# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE)
track <- cbind(track,rowSums(seqtab.nochim))
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab2)

print("Assigning taxonomy")
taxa_g <- assignTaxonomy(seqtab.nochim, "/opt/ancillary/databases/dada2/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
save.image("run5.RData")

print("Assigning species")
taxa <- addSpecies(taxa_g, "/opt/ancillary/databases/dada2/silva_species_assignment_v138.fa.gz")

write.csv(seqtab.nochim, "seqtab.nonchim.csv")
write.csv(taxa, "taxa.csv")

# Install Phyloseq package in R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("phyloseq")

library("phyloseq"); packageVersion("phyloseq")

# Generate Phuloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(taxa))
seqtab.nochim
ps

# Save R Data
save.image("ncgs.RData")
# Save Phyloseq Object
saveRDS(ps.run5, file = "ps.rds")

# Use this Phyloseq Object to generate plots in R
