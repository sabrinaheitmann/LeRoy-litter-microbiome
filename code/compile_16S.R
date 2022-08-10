#V4
#Compile MiSeq library and sample data into single phyloseq object 
#Also add taxonomy and remove non-target sequences

library(tidyverse)
library(magrittr)
library(foreach)
library(dada2)
library(Biostrings)
library(phyloseq)

#read in sample meta data
meta <- read.csv("data/meta/LeRoy_litter_metadata_3-14-22.csv",as.is=T,row.names=1) #row.names = column with sample names that match OTU table

#read in denoised sequence table
seqTab <- readRDS("output/dada/16S/bact/seqTab.rds")
#row.names(seqTab) %<>% str_remove("\\.$") # This needs to be fixed in denoise.R!

#extract sequences from OTU table 
seqs <- getSequences(seqTab) %>% dada2::rc() %>% DNAStringSet

#Creat data frame for summary of processing
compile.summary <- data.frame(SampID=rownames(seqTab),
                              denoised=rowSums(seqTab))
#rename sequences
OTU.names <- paste0("OTU.",1:ncol(seqTab))
colnames(seqTab) <- OTU.names
names(seqs) <- OTU.names

#Create compile/16S.scratch directory
dir.create("output/compile/16S.scratch")
#write temp file
writeXStringSet(seqs,file="output/compile/16S.scratch/tmp.fasta",width=1000)

#change out to Salix genome if needed
#remove host contamination with Bowtie2 - uses a Ptrichocarpa v3 genome assembly pre-processed with bowtie2-build
#system("/home/busbylab/miniconda3/bin/bowtie2 --threads 10 -x data/Ptri_genome/CopyOfPtri.v.3.db -f output/compile/16S.scratch/tmp.fasta --un output/compile/16S.scratch/noHost.fa --al output/compile/16S.scratch/host.fa --very-sensitive-local")
#seqs.noHost <- readDNAStringSet("output/compile/16S.scratch/noHost.fa")
#compile.summary$noHost <- seqTab %>% .[,names(seqs.noHost)] %>% rowSums()

#convert to phyloseq object
#changed to tmp.fasta from noHost.fa
seqs.noHost <- readDNAStringSet("output/compile/16S.scratch/tmp.fasta") 

otuTab <- seqTab %>% as.data.frame() %>% otu_table(taxa_are_rows = F)
phy <- phyloseq(otuTab, refseq(seqs))

#collapse sequences with >= 99% similarity
cluster <- function(phy.in,method="UPGMA",dissimilarity=0.01){
  require(DECIPHER)
  clusters <- DistanceMatrix(refseq(phy.in), includeTerminalGaps = T, processors=NULL) %>%
    IdClusters(method=method, cutoff=dissimilarity, processors=NULL) %>%
    rownames_to_column("OTUid")
  for(i in unique(clusters$cluster)){
    foo <- clusters$OTUid[which(clusters$cluster==i)] 
    if(length(foo)>1){phy.in %<>% merge_taxa(foo)}
  }
  return(phy.in)
}  
phy <- phy %>% cluster  

#Re-assign taxonomy with bac-only database
#the rc of each sequence was used if it is a better match to the reference sequences than the forward sequence
taxa <- assignTaxonomy(refseq(phy),"data/taxonomy_db/silva_nr_v132_train_set.fa.gz", multithread = T, minBoot=50, tryRC = TRUE)
rownames(taxa) %<>% names

#Add taxonomy to phyloseq object
tax_table(phy) <- taxa %>% as.matrix %>% tax_table

#convert meta to phyloseq sample data
meta %<>% sample_data

#add meta data to phyloseq object
phy %<>% merge_phyloseq(meta)

#removed Eukaryota
phy <- subset_taxa(phy, !Kingdom=="Eukaryota")

#Save phyloseq object
saveRDS(phy,"output/phy_16S/phy_16S.rds")

#Save components for possible manual inspection
otu_table(phy) %>% write.csv("output/phy_16S/16S.OTU.table.csv")
tax_table(phy) %>% write.csv("output/phy_16S/16S.taxonomy.table.csv")
sample_data(phy) %>% write.csv("output/phy_16S/16S.sample.data.csv")
