#ITS2
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
seqTab <- readRDS("output/dada/ITS2/seqTab.rds")
row.names(seqTab) %<>% str_remove("\\.$") # This needs to be fixed in denoise.R!

#extract sequences from OTU table 
seqs <- getSequences(seqTab) %>% dada2::rc() %>% DNAStringSet

#Creat data frame for summary of processing
compile.summary <- data.frame(SampID=rownames(seqTab),
                              denoised=rowSums(seqTab))
#rename sequences
OTU.names <- paste0("OTU.",1:ncol(seqTab))
colnames(seqTab) <- OTU.names
names(seqs) <- OTU.names

### Assign taxonomy ###
#download a version of the UNITE database if you dont have it
if(!file.exists("data/taxonomy_db")){
  download.file("https://files.plutof.ut.ee/public/orig/EB/0C/EB0CCB3A871B77EA75E472D13926271076904A588D2E1C1EA5AFCF7397D48378.zip",
                "data/taxonomy_db.zip",
                method="wget")
  unzip("data/taxonomy_db.zip",exdir="data/taxonomy_db")
  file.remove("data/taxonomy_db.zip")}

#Extract ITS region with ITSx
#Create compile/ITS2.scratch directory
dir.create("output/compile/ITS2.scratch")
#write temp file
writeXStringSet(seqs,file="output/compile/ITS2.scratch/tmp.fasta",width=1000)

#system("/home/busbylab/miniconda3/bin/bowtie2 --threads 10 -x data/Ptri_genome/CopyOfPtri.v.3.db -f output/compile/ITS2.scratch/tmp.fasta --un output/compile/ITS2.scratch/noHost.fa --al output/compile/ITS2.scratch/host.fa --very-sensitive-local")

#changed to tmp.fasta from noHost.fa
seqs.noHost <- readDNAStringSet("output/compile/ITS2.scratch/tmp.fasta")

#summary of denoised and nohost sequence counts
compile.summary$noHost <- seqTab %>% .[,names(seqs.noHost)] %>% rowSums()

# ITSx parameters (these will depend on your data)
#changed to tmp.fasta from noHost.fa
ITSx.flags <- paste("-i output/compile/ITS2.scratch/tmp.fasta",
                    "-t 'fungi'",
                    "--preserve T",
                    "--complement F",
                    "--summary T",
                    "--cpu 10",
                    "-o output/compile/ITS2.scratch/ITSx",
                    "--only_full T",
                    "-E 1e-2")

system2("ITSx", args = ITSx.flags)

#remove OTUs from seqTab that are not in ITSx output (and longer than 75 bp)
seqs.ITS <- readDNAStringSet("output/compile/ITS2.scratch/ITSx.ITS2.fasta") %>% .[.@ranges@width > 75]
#remove OTUs from seqTab that are not in ITSx output
seqTab %<>% .[,names(seqs.ITS)] 
#remove OTUs from seqTab that are not in ITSx output
colnames(seqTab) <- seqs.ITS %>% as.character 
seqTab %<>% collapseNoMismatch()
compile.summary$ITSx.filtered <- rowSums(seqTab)

###############
### OUTPUTs ###
###############

#set final OTU names
final.names <- paste0("OTU.",1:ncol(seqTab))

final.seqs <- getSequences(seqTab) %>% DNAStringSet()
names(final.seqs) <- final.names
colnames(seqTab) <- final.names
#rownames(tax) <- final.names does this make a difference if hashed?

#make phyloseq object
phy <- phyloseq(otu_table(seqTab,taxa_are_rows = F), 
                sample_data(meta),
                refseq(final.seqs))

cluster <- function(phy.in,method="single",dissimilarity=0.01){
  require(DECIPHER)
  clusters <- DistanceMatrix(refseq(phy.in), includeTerminalGaps = T, processors=NULL) %>%
    IdClusters(method=method, cutoff=dissimilarity, processors=NULL) 
  clusters[,1] %<>% as.character()
  tax_table(phy.in) <- clusters %>% as.matrix %>% tax_table
  phy.in %<>% speedyseq::tax_glom("cluster")
  tax_table(phy.in) <- NULL
  return(phy.in)
}  
phy %<>% cluster  

#assign taxonomy with DADA2 implementation of RDP classifier (needs 12+ gb ram)
tax <- assignTaxonomy(refseq(phy),"data/taxonomy_db/sh_general_release_dynamic_02.02.2019.fasta",multithread=T)
rownames(tax) %<>% names

#Add taxonomy to phyloseq object
tax_table(phy) <- tax %>% as.matrix %>% tax_table

#Save phyloseq object
saveRDS(phy,"output/phy_ITS2/phy_ITS2.rds")

#Save components for possible manual inspection
otu_table(phy) %>% write.csv("output/phy_ITS2/ITS2.OTU.table.csv")
tax_table(phy) %>% write.csv("output/phy_ITS2/ITS2.taxonomy.table.csv")
sample_data(phy) %>% write.csv("output/phy_ITS2/ITS2.sample.data.csv")
