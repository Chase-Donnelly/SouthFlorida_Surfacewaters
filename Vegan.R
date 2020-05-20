#start with loading data
dat <- read.table(file= "feature-table_Rfix.txt", header = TRUE,sep ="\t", row.names = 1)
#look at imported file 
View(dat)
#transpose the data to rows 
t.dat <- as.data.frame(t(dat))
##view first rows 
t.dat[1:5,1:5]
dat <-t.dat
#import the metadata and view 
metadata <- read.table(file= "cd_metadata_water_chem_test.txt", header=T,sep = "\t", row.names=1)
#view to check 
View(metadata)
#now we need to make it so we only have the data for the specific rows we are looking at 
common.rownames <- intersect(rownames(dat),rownames(metadata))
dat <- dat[common.rownames,]
metadata <- metadata[common.rownames,]
all.equal(rownames(dat),rownames(metadata))
#reduce noise (get rid of sigle and doubletons)
otu.abund<-which(colSums(dat)>2)
dat.dom<-dat[,otu.abund]
#reduce otus that occur in small amount of samples
library(vegan)
library(base)
dat.pa<-decostand(dat.dom, method ="pa")
dat.otus.05per<-which(colSums(dat.pa) > (0.05*nrow(dat.pa)))
dat.05per<-dat.dom[,dat.otus.05per]
#transform data for relative abundance 
dat.ra<-decostand(dat.05per, method = "total")
#print to excell sheet for kronos if using 
dat.rat <- as.data.frame(t(dat.ra))
View(dat.rat)
write.table(dat.rat, file = "Water_Kronos", sep="\t",row.names = T) 
# look at bray curtis dissimilarity 
dat.bc.dist<-vegdist(dat.ra, method = "bray")
#adonis 
adonis(dat.bc.dist~SiteType*SPCONDUCTIVITYFIELD, data = metadata)
#run a pcoa for adonis results 
dat.betadisp<-betadisper(dat.bc.dist,metadata$SiteType)
boxplot(dat.betadisp)
plot(dat.betadisp) ##(use this plot for slides)

#run pairwise adonis 
library(RVAideMemoire)
pairwise.perm.manova(dat.bc.dist,metadata$SiteType)
#simper 
dat.simp<-simper(dat.ra, metadata$SiteType, permutations = 999)##change to 999 after intial run 
sink("Simper_by_siteTYPE.csv")
summary(dat.simp)
sink()
##look at the file and you can see what otus are causeing the difference between the sites, look up the otu and see if that is interesting 


##CCA time 
set.seed(42);env.cca<-cca(dat.ra~Rain3days+Rain1month+PHOSPHATETOTALASP+TEMP+DISSOLVEDOXYGEN+SPCONDUCTIVITYFIELD+PHFIELD+TOTALNITROGEN+AMMONIA.N, data =metadata) #+SODIUM+POTASSIUM+CALCIUM+MAGNESIUM+CHLORIDE+SULFATE+
env.cca
vif.cca(env.cca)
#make sure they add up to more than ten or you may need to remove if its over 20 def remove
#step 2, zero variables
set.seed(42);lwr<- cca(dat.ra~1, data=metadata)
lwr
#forward selecting model 
set.seed(42);mods.all<- ordiR2step(lwr, scope = formula(env.cca))
mods.all
vif.cca(mods.all)

R2.adj.all<-RsquareAdj(mods.all)
R2.adj.all

mods.all$anova
#repeat this for different sites to see if the variance is different for each site (to do this just change the metadata file)

## try ploting this 
cca.p <- plot(mods.all,type = "none")
points(cca.p, "sites", col= as.numeric(metadata$SiteType), pch = as.numeric(metadata$SiteType))

ef.all<- envfit(cca.p,metadata[,c("SPCONDUCTIVITYFIELD","TEMP","PHOSPHATETOTALASP","Rain1month","DISSOLVEDOXYGEN")])
plot(ef.all)

legend("topright",legend = as.character(paste(" ", unique(metadata$SiteType))), pch= as.numeric(unique(metadata$SiteType)))

ordiellipse(cca.p, metadata$SiteType, label = T, conf = 0.95)

#looking into ndms chart 
set.seed(42)
comm.bc.mds<-metaMDS(dat.ra, distance="bray")
mds.fig<-ordiplot(comm.bc.mds, display="sites")
ordiellipse(mds.fig, metadata$SiteType, label = T, conf = 0.95)
#Adding environmental and trait data to ordinations
plot(envfit(comm.bc.mds, metadata[,c(5:8)]))
##this is how you can adjust the x and y axis 
mds.fig<-ordiplot(comm.bc.mds, display="sites", xlim=c(-2,3), ylim = c(-2,3))  #### adjust xlim and ylim
##adjust colors
points(mds.fig,"sites", pch = 15, col = "black", select = metadata$SiteType == "Ag")
points(mds.fig,"sites", pch = 16, col = "red", select = metadata$SiteType == "Urban")
points(mds.fig,"sites", pch = 17, col = "green3", select = metadata$SiteType == "Control")
points(mds.fig,"sites", pch = 17, col = "blue", select = metadata$SiteType == "UrbanSalt")
#legend
legend("topright",legend=as.character(paste(" ",unique(metadata$SiteType))), cex = 0.99,pch=19,col=1:length(unique(metadata$SiteType)))
ordiellipse(mds.fig, metadata$SiteType, label = F, conf = 0.95, lty = 2)
palette()

##Tukey Test
comm.bc.dist<-vegdist(dat, method = "bray")
comp.variation<-betadisper(comm.bc.dist, metadata$SiteType)
TukeyHSD(comp.variation)
##
comp.variation2<-betadisper(comm.bc.dist, metadata$Date) 
boxplot(comp.variation2)
##
comm<- decostand(dat,method = "total")
diversity<- diversity(comm,index="invsimpson")
diversity
div.table <- as.data.frame(diversity)          
cor.table <- cbind(metadata,div.table)
plot(diversity~metadata$SiteType)
sit.aov <- aov(diversity~metadata$SiteType)
summary(sit.aov)
TukeyHSD(sit.aov)
