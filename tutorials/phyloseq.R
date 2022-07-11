#' Following https://vaulot.github.io/tutorials/Phyloseq_tutorial.html

library(dplyr)
library(readxl)
library(ggplot2)
library(phyloseq)

#' Read in data
otu <- read_excel("data/CARBOM data.xlsx", sheet = "OTU matrix")
tax <- read_excel("data/CARBOM data.xlsx", sheet = "Taxonomy table")
samples <- read_excel("data/CARBOM data.xlsx", sheet = "Samples")

otu <- tibble::column_to_rownames(otu, "otu")
tax <- tibble::column_to_rownames(tax, "otu")
samples <- tibble::column_to_rownames(samples, "sample")

#' Create phyloseq class object
carbom <- phyloseq::phyloseq(
  phyloseq::otu_table(as.matrix(otu), taxa_are_rows = TRUE),
  phyloseq::tax_table(as.matrix(tax)),
  phyloseq::sample_data(samples)
)
carbom

sample_names(carbom)
rank_names(carbom)

carbom <- subset_samples(carbom, Select_18S_nifH == "Yes")
carbom #' Only lost one sample here by the looks of things

#' Keep only photopynthetic taxa (ok?)
carbom <- carbom %>%
  subset_taxa(
    Division %in% c("Chlorophyta", "Dinophyta", "Cryptophyta", "Haptophyta", "Ochrophyta", "Cercozoa"),
    !(Class %in% c("Syndiniales", "Sarcomonadea"))
  )
carbom

#' Normalise the number of reads in each sample using median sequencing depth
total <- median(sample_sums(carbom))
standf <- function(x, t = total) round(t * (x / sum(x)))
carbom <- transform_sample_counts(carbom, standf)

plot_bar(carbom, fill = "Division") +
  geom_bar(aes(col = Division, fill = Division), stat = "identity", position = "stack")

#' "Regroup together pico vs nano samples" (I don't know what this means really)
carbom_fraction <- merge_samples(carbom, "fraction")

plot_bar(carbom_fraction, fill = "Division") +
  geom_bar(aes(col = Division, fill = Division), stat = "identity", position = "stack")

#' Keep only Chlorophyta and colour according to genus. Facet by pico vs nano and surface vs deep samples
carbom_chloro <- subset_taxa(carbom, Division == "Chlorophyta")

plot_bar(carbom_chloro, x = "Genus", fill = "Genus", facet_grid = level ~ fraction) +
  geom_bar(aes(col = Genus, fill = Genus), stat = "identity", position = "stack")

#' Heatmaps
#' What is method NMDS?
#' What is distance bray?
plot_heatmap(carbom, method = "NMDS", distance = "bray")

#' Apparently the above is "cluttered" so they suggest filtering out those OTUs which represent
#' less than 20% of the reads in at least one sample
carbom_abund <- filter_taxa(carbom, function(x) sum(x > total * 0.20) > 0, TRUE)
carbom_abund

otu_table(carbom_abund)

plot_heatmap(carbom_abund, method = "NMDS", distance = "bray")

#' It's possible to use other distances and methods
#' This one is called Jaccard distance (A + B - 2 * J / (A + B - J))
plot_heatmap(
  carbom_abund, method = "MDS", distance = "(A+B-2*J)/(A+B-J)",
  taxa.label = "Class", taxa.order = "Class", trans = NULL
)

#' Many different distances are built-in
(dist_methods <- unlist(distanceMethodList))

#' You can also construct your own
#' J = sum(x * y)
#' A = sum(x^2)
#' B = sum(y^2)

#' Alpha diversity
plot_richness(carbom, measures = c("Chao1", "Shannon"))
plot_richness(carbom, measure = c("Chao1", "Shannon"), x = "level", col = "fraction")

#' Ordination
carbom_ord <- ordinate(carbom, "NMDS", "bray")

plot_ordination(carbom, carbom_ord, type = "taxa", col = "Class", shape = "Division", title = "OTUs")

#' Facet by taxonomic division
plot_ordination(carbom, carbom_ord, type = "taxa", col = "Class", title = "OTUs", label = "Class") +
  facet_wrap(~ Division, 3)

plot_ordination(carbom, carbom_ord, type = "samples", col = "fraction", shape = "level", title = "Samples") +
  geom_point(size = 3)

plot_ordination(carbom, carbom_ord, type = "split", col = "Class", shape = "level", title = "biplot", label = "station") +
  geom_point(size = 3)

#' Network analysis
plot_net(carbom, distance = "(A+B-2*J)/(A+B)", type = "taxa", maxdist = 0.7, col = "Class", point_label = "Genus")
plot_net(carbom_abund, distance = "(A+B-2*J)/(A+B)", type = "taxa", maxdist = 0.8, col = "Class", point_label="Genus")
