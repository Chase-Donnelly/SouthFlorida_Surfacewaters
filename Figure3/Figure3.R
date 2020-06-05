# this code is used to create figure 3
# data for this code is found in figure 2 
##load packages and set the directory 
library(ggplot2)
library(phyloseq)
library(ape)
library(vegan)
library(ggplot2)
library(tidyverse)
library(ggsci)
library(RColorBrewer)
# set our color pallete 
nb.cols <- 13
mypal = pal_npg("nrc", alpha = 0.7)(9) 
mycolors <- colorRampPalette(mypal)(nb.cols)
#read in otu table 
otu_table=read.table(file= "feature-table_Rfix.txt", header=TRUE, sep ="\t", row.names = 1)
otu_table=as.matrix(otu_table)
##read in taxonomy 
#make sure these are seperated columns for kpcofgs
taxonomy=read.table(file = "taxonomy_R_sep.txt", sep = "\t", header = T, row.names = 1)
head(taxonomy)
taxonomy=as.matrix(taxonomy)
##add metadata
metadata=read.table("cd_metadata_water_final_filled.txt", header=T, sep = "\t", row.names = 1)
##load tree
phy_tree=read_tree("tree-unrooted.nwk")
##import as phyloseq objects
OTU= otu_table(otu_table,taxa_are_rows=TRUE)
TAX=tax_table(taxonomy)
META=sample_data(metadata)
##check that you OTU names are consistent across objects 
taxa_names(TAX)
taxa_names(OTU)
taxa_names(phy_tree)
##merge into one phyloseq object
physeq = phyloseq(OTU,TAX,META,phy_tree)
physeq
##check rank names of taxonomy
rank_names(physeq)
##prune taxa from the OTU table that are in zero samples (these are in other samples on the run)
merge=prune_taxa(taxa_sums(physeq)>0,physeq)
merge
##create for taxa above RA of 1% 
merge99 = transform_sample_counts(merge, function(x){x/sum(x)})
otu_table(merge99)[otu_table(merge99)<.01] <- 0 
merge99 = prune_taxa(taxa_sums(merge99)>0,merge99)
merge99 = transform_sample_counts(merge99, function(x){x*100})
otu_table(merge99) = floor(otu_table(merge99))
merge99
#create a normalized data set for lowest reads 
# Normalize to 24381 reads per sample (proportions) rounding down
mnorm = transform_sample_counts(physeq, function(x) {24381*x/sum(x)})
otu_table(mnorm) = floor(otu_table(mnorm)) 
mnorm = prune_taxa(taxa_sums(mnorm)>0,mnorm)
mnorm

#networking now to make the network 
set.seed(711L)
ig <- make_network(merge, dist.fun="bray", max.dist=0.7)
net <- plot_network(ig, merge, color="SampleLocation", shape="SiteType", line_weight=0.4, label = NULL)
# set the theme for the graph 
mynamestheme <- theme(plot.title = element_text(family = "Helvetica", face = "bold", size = (15), hjust = 0.5), 
                      legend.title = element_text(colour = "steelblue",  face = "bold.italic", family = "Helvetica"), 
                      legend.text = element_text(face = "italic", colour="steelblue4",family = "Helvetica"), 
                      axis.title = element_text(family = "Helvetica", size = (10), colour = "steelblue4", face = "bold"),
                      axis.text = element_text(family = "Courier", colour = "cornflowerblue", size = (10), face = "bold"))
# finnally print our graph
print(net + mynamestheme  + scale_color_manual(values = mycolors) + ggtitle("Network Analysis: Bray Curtis Disimilarity (0.7 max dist)")) 


      