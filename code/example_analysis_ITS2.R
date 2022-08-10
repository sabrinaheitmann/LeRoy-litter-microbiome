library(tidyverse)
library(magrittr)
library(phyloseq)
library(vegan)

### Read in phyloseq object
(phy <- readRDS("output/phy_ITS2/phy_ITS2.rds"))

#' ### Removed NLPA and NLPW bc no Weevil info (which causes the perMANOVA not to run)
phy = subset_samples(phy, sample != "NLPA" & sample != "NLPW")

#' ### Check meta data to confirm controls are removed
meta <- phy %>% sample_data %>% data.frame()

minDepth <- 250
data.frame(SeqDepth=sort(sample_sums(phy)), Study=sample_data(phy)$Study) %>%
  mutate(cutoff=SeqDepth>minDepth) %>%
  ggplot(aes(x=Study, y=SeqDepth)) +
  geom_violin() +
  geom_point(aes(color=cutoff),position=position_jitter(width=0.1)) +
  theme_classic()

#' ### Remove samples below sequencing depth cutoff
(phy %<>% prune_samples(sample_sums(.)>minDepth,.))

#' ### Remove low abundance OTUs by prevelance
# Only keep OTUs present in at least 1% of samples
(phy %>% filter_taxa(., function(x) {sum(x>0) > 0.01*nsamples(.)}, TRUE))

#' ### Convert to proportional abundance
phy %<>% transform_sample_counts(function(x){x*min(sample_sums(.)/sum(x))})

#' ## Community analyses

#' ### Plot an ordination (PCoA)
phy.ord <- phy %>% ordinate("MDS","bray")
plot_ordination(phy,phy.ord,color="Weevil") + theme_classic()
#Save image to file
ggsave("output/figs/Weevil_ITS2_PCoA.pdf")

#' ### Run a perMANOVA
phy.dist <- phy %>% phyloseq::distance("bray")  
adonis(phy.dist~sample_data(phy)$Weevil)
