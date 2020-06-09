# This is an attempt to make a CCA plot with ggvegan package 
#First load the needed libraries 
library(vegan)
library(base)
library(ggvegan)
library(tidyverse)
library(ggsci)
library(RColorBrewer)
#start with loading SV data 
dat <- read.table(file= "feature-table_Rfix.txt", header = TRUE,sep ="\t", row.names = 1)
#transpose the data to rows 
dat <- as.data.frame(t(dat))
#import the metadata and view 
metadata <- read.table(file= "cd_metadata_water_final_filled_nosalt.txt", header=T,sep = "\t", row.names=1)
#now we need to make it so we only have the data for the specific rows we are looking at 
common.rownames <- intersect(rownames(dat),rownames(metadata))
dat <- dat[common.rownames,]
metadata <- metadata[common.rownames,]
all.equal(rownames(dat),rownames(metadata))
#reduce noise (get rid of sigle and doubletons)
otu.abund<-which(colSums(dat)>2)
dat.dom<-dat[,otu.abund]
#reduce otus that occur in small amount of samples
dat.pa<-decostand(dat.dom, method ="pa")
dat.otus.05per<-which(colSums(dat.pa) > (0.05*nrow(dat.pa)))
dat.05per<-dat.dom[,dat.otus.05per]
#transform data for relative abundance 
dat.ra<-decostand(dat.05per, method = "total")
# CCA Time
set.seed(42);ord <- cca(dat.ra ~Month+Latitude+Longitude+WaterTemp+AirTemp+Rain3days+Rain1month+TP+DO+CONDUCTIVITY+pH+TN+TEMP+AMMONIA+TSS, data =metadata)
ord
# check the numbers!
vif.cca(ord)
#make sure they add up to more than ten or you may need to remove if its over 20 def remove
#step 2, zero variables
set.seed(42);lwr<- cca(dat.ra~1, data=metadata)
lwr
#forward selecting model, this is to go through and remove unneeded variables
set.seed(42);mods.all<- ordiR2step(lwr, scope = formula(ord))
mods.all
vif.cca(mods.all)
# 
mods.all <- cca(dat.ra ~Month+Latitude+TSS+TN+TP+DO, data = metadata)
mods.all
# now to take the arrows by fortifying
#ford <- fortify(mods.all, axes = 1:2)  # fortify the ordination
#take <- c('CCA1', 'CCA2')  # which columns contain the scores we want
#arrows <- subset(ford, Score == 'biplot')  # take only biplot arrow scores
## multiplier for arrows to scale them to the plot range
#mul <- ggvegan:::arrowMul(arrows[, take],
                         # subset(ford, select = take, Score == 'sites'))
#arrows[, take] <- arrows[, take] * mul  # scale biplot arrows

#now try plotting attempt 2 
# set some parameters 
#ef.all <- envfit(mods.all,metadata[,c("CONDUCTIVITY","TP","DO","Latitude")])
site_type <- metadata %>% 
  select(SampleLocation, SiteType)
pal <- c("lightsalmon1", "gold1", "palegreen4")
clean_background <- theme(plot.background = element_rect("white"),
                          panel.background = element_rect("white"),
                          panel.grid = element_line("white"),
                          axis.line = element_line("gray25"),
                          axis.text = element_text(size = 12, color = "gray25"),
                          plot.title = element_text(family = "Helvetica", face = "bold", size = (15), hjust = 0.5), 
                          legend.title = element_text(colour = "steelblue",  face = "bold.italic", family = "Helvetica"), 
                          legend.text = element_text(face = "italic", colour="steelblue4",family = "Helvetica"), 
                          axis.title = element_text(family = "Helvetica", size = (10), colour = "steelblue4", face = "bold"),)

#
ccaplot <- plot(mods.all)
ccavectors <- as.data.frame(ccaplot$biplot * 5.15)

site_data <- as.data.frame(ccaplot$sites) %>% 
  bind_cols(., site_type)

species_data <- as.data.frame(ccaplot$species)

plot_cca <- ggplot(site_data) +
  geom_point(aes(x = CCA1, y = CCA2, color = SiteType), shape = 19, size = 2, alpha = 0.8) +
  scale_color_npg() +
  geom_segment(data = ccavectors, aes(x = 0, y = 0, xend = CCA1, yend = CCA2), arrow = arrow(length = unit(0.2, "cm"))) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  scale_x_continuous(limits = c(-6, 6)) +
  geom_text(data = ccavectors, aes(x = CCA1, y = CCA2, label = rownames(ccavectors)), nudge_x = 0.3, nudge_y = 0.3) +
  clean_background +
  labs(title = "Canonical Correspondence Analysis")
print(plot_cca)  





