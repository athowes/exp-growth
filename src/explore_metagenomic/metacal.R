#' Following https://mikemc.github.io/metacal/articles/tutorial.html

library(phyloseq)
library(tidyverse)
library(ggbeeswarm)
library(cowplot)
library(metacal)

theme_set(theme_cowplot())

#' Black cottonwood was inoculated with fungal rust pathogen Melampsora x columbiana and
#' 8 other species of follar fungi. Fungi relative abundances were measured using ITS
#' amplicon sequencing. Also created a set of DNA mocks of these fungi, which were measured
#' along with the experimental samples.
#'
#' Approach: use the mocks to estimate the bias of the 9 focal species, then use this to work
#' out the species abundance in the full set of samples.

data_path <- here::here("src", "explore_metagenomic", "data", "leopold2020host")

if(!dir.exists(data_path)) {
  dir.create(data_path, recursive = TRUE)
  download.file(
    "https://zenodo.org/record/3872145/files/dleopold/Populus_priorityEffects-v1.2.zip",
    file.path(data_path, "Populus_priorityEffects-v1.2.zip")
  )
  unzip(
    file.path(data_path, "Populus_priorityEffects-v1.2.zip"),
    exdir = data_path
  )
}

ps <- file.path(
  data_path,
  "dleopold-Populus_priorityEffects-8594f7c/output/compiled/phy.rds"
) %>%
  readRDS() %>%
  print()

#' OTU: operational taxonomic unit
#' "In biology, a taxon (back-formation from taxonomy; plural taxa) is a group of one or
#'  more populations of an organism or organisms seen by taxonomists to form a unit."

#' close_etls: Close the elements of x to proportions. Could just say "normalise" (?) as
#' it's just x / sum(x, na.rm = na.rm)

mock_actual <- file.path(
  data_path,
  "dleopold-Populus_priorityEffects-8594f7c/data/MockCommunities.csv"
) %>%
  read.csv(row.names = 1) %>%
  select(-Sym4) %>%
  as("matrix") %>%
  otu_table(taxa_are_rows = FALSE) %>%
  #' Reciprocals of the dilution factors
  transform_sample_counts(function(x) close_elts(1 / x))

#' "actual" means they are the (somewhat) correct answers
#' "mock" means they are control samples

mock_taxa <- taxa_names(mock_actual)
mock_taxa #' These are the taxa names
#' * Melampsora is a genus of Basidiomycota fungi -- they're a plant pathogen
#' * Dioszegia is a genus of fungi in the family Bulleribasidiaceae
#' * ... not going to look up every fungi genus on Wikipedia

#' Observed abundances are read counts
otu_table(ps) %>%
  #' This is calling prune_taxa(mock_taxa, otu_table(ps))
  #' I think it's filtering for those taxa included in mock_taxa
  prune_taxa(mock_taxa, .) %>%
  #' Displays a corner section of a rectangular data set
  corner()

sam <- sample_data(ps) %>%
  as("data.frame") %>%
  as_tibble(rownames = "Sample")

glimpse(sam)

#' The control samples are the n = 10 with Samp_type Mock
#' What are the n = 9 with Samp_type Single?
count(sam, Samp_type)

#' Taxnonomy table
tax <- tax_table(ps) %>%
  as("matrix") %>%
  as_tibble(rownames = "Taxon")

#' It looks a bit strange to me that you would need the "k__" in front of the names
#' within "Kingdom" or the "p__" in front of the names within "Phylum" (as in it's a
#' duplication of information)
head(tax)

#' These are the only taxa which we have calibration information about from the set
#' of mock experiments which were ran
tax %>%
  filter(Taxon %in% mock_taxa) %>%
  select(Taxon, Family:Species)

#' There are in total 219 taxa, of which 9 are assigned to the mock taxa
ntaxa(ps)

#' We are now going to inspect the mock community measurements
#' All of the taxa are in every sample, in proportions which vary from 0.01 to 0.3
mock_actual %>%
  psmelt() %>%
  mutate(
    across(OTU, factor, levels = mock_taxa),
    across(Sample, factor, levels = sample_names(mock_actual))
  ) %>%
  ggplot(aes(x = Sample, y = OTU, fill = Abundance)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log10", breaks = c(0.02, 0.05, 0.1, 0.2))

ps_mock <- ps %>%
  subset_samples(Samp_type == "Mock") %>%
  prune_taxa(mock_taxa, .)

ps_mock
