#this is the code used to create figure 1 A & B 
#import data, remember to set your directory first
data <- read.csv('taxa_barplots_final_corrected.csv', header=T, sep = ";", row.names=1)
#view the data to check for any errors during import
View(data)
##start libraries 
library(ggplot2)
library(tidyverse)
library(ggsci)
library(RColorBrewer)
# create an extended color palete, i am using ggsci color palettes and extending the number of colors to allow for use in the charts
nb.cols <- 87
mypal = pal_npg("nrc", alpha = 0.7)(9) 
mycolors <- colorRampPalette(mypal)(nb.cols)
# here we set our main plot using ggplot2 
phylum <- ggplot() + geom_bar(aes(y= Percentage,
                        x = Site_Type, 
                        fill = Phylum),
                    data = data,
                    stat = 'identity',
                    position = "fill")
# create a theme, this is used to add styles to the plot
mynamestheme <- theme(plot.title = element_text(family = "Helvetica", face = "bold", size = (15), hjust = 0.5), 
                      legend.title = element_text(colour = "steelblue",  face = "bold.italic", family = "Helvetica"), 
                      legend.text = element_text(face = "italic", colour="steelblue4",family = "Helvetica"), 
                      axis.title = element_text(family = "Helvetica", size = (10), colour = "steelblue4", face = "bold"),
                      axis.text = element_text(family = "Courier", colour = "cornflowerblue", size = (10), face = "bold"))
# now we are going to print the chart as well as adding the color, theme and titles/labels 
print(phylum + mynamestheme + scale_fill_manual(values = mycolors) + ggtitle("Comparing Top 0.1% Phylum \nPer Sampling Location") +labs(y="Proportion of Top Phylum", x="Sampling Location"))

# now we create the plot for order 
order <- ggplot() + geom_bar(aes(y= Percentage,
                                  x = Site_Type, 
                                  fill = Order),
                              data = data,
                              stat = 'identity',
                              position = "fill")
#our theme and colors are already set so we can just print the graph with new titles 
print(order + mynamestheme + scale_fill_manual(values = mycolors) + ggtitle("Comparing Top 0.1% Order \nPer Sampling Location") +labs(y="Proportion of Top Orders", x="Sampling Location") + aes(position = "fill"))
#now for family
family <- ggplot() + geom_bar(aes(y= Percentage,
                                           x = Site_Type, 
                                           fill = Family),
                                       data = data,
                                       stat = 'identity',
                                       position = "fill")
# we adjust the theme alittle so we need to run again 
mynamestheme <- theme(plot.title = element_text(family = "Helvetica", face = "bold", size = (15), hjust = 0.5), 
                      legend.title = element_text(colour = "steelblue",  face = "bold.italic", family = "Helvetica"), 
                      legend.text = element_text(face = "italic", colour="steelblue4",family = "Helvetica", size = (4)), 
                      axis.title = element_text(family = "Helvetica", size = (10), colour = "steelblue4"),
                      axis.text = element_text(family = "Courier", colour = "cornflowerblue", size = (10)))
# and print 
print(family + mynamestheme + scale_fill_manual(values = mycolors) + ggtitle("Comparing Top 0.1% Family \nPer Sampling Location") +labs(y="Proportion of Top Families", x="Sampling Location") + aes(position = "fill"))













