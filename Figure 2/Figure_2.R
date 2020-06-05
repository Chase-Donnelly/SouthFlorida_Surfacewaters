# this code is used to create figure 2 
##load packages and set the directory 
library(ggplot2)
library(phyloseq)
library(ape)
library(vegan)
###now to import to phyloseq
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
##now continue analysis in phyloseq
## check reads of samples
sample_sums(physeq)[1:10]
## basic stats for read of samples
mean(sample_sums(physeq))
min(sample_sums(physeq))
max(sample_sums(physeq))
sd(sample_sums(physeq))

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

##Alpha diversity 
q = plot_richness(merge, x="SampleLocation", color = "SiteType", measures = c("Shannon","InvSimpson"))
q + geom_boxplot(data = q$data, aes(x=SampleLocation, y=value, color=NULL, fill=NULL ),alpha=0.1) ##+ geom_point(size =3, alpha=0.7)

#now to run an anova on our results, first we save the alpha values 
alpha.diversity <- estimate_richness(merge, measures = c("Shannon","InvSimpson"))
head(alpha.diversity)
write.table(alpha.diversity, "alpha.txt")
#next to run the anova for shannon
data <- cbind(sample_data(merge), alpha.diversity)
merge.anova <- aov(Shannon ~ SampleLocation, data)
summary(merge.anova)
#
data <- cbind(sample_data(merge), alpha.diversity)
merge.anova <- aov(InvSimpson ~ SiteType, data)
summary(merge.anova)















