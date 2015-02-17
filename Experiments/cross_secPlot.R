### 
# Creates the cross-section plots
###

library(ggplot2)
library(dplyr)
library(extrafont) #for latex fonts
library(reshape2)

setwd("/Users/vishalgupta/Documents/Research/BayesDRO/AmbiguitySets/Experiments/Results")
#font = "CM Roman"
font = "Times New Roman"

#merge the relevant pieces
dat = rbind(read.csv("secMom50.csv", header=TRUE), 
            read.csv("bernvar50.csv", header=TRUE), 
            read.csv("pstar50.csv", header=TRUE))

# dat = rbind(read.csv("secMom500.csv", header=TRUE), 
#             read.csv("bernvar500.csv", header=TRUE), 
#             read.csv("pstar500.csv", header=TRUE))


g <- ggplot(aes(x=p1, y=p2, linetype=Type, color=Type), data=dat) + 
  geom_path() + 
  xlab(expression(p[1])) + ylab(expression(p[2])) +
  theme_bw(base_size=12) + 
  theme(legend.title=element_blank(), 
        legend.direction="horizontal",
        text=element_text(family=font))

mylabs <- c(expression(chi^2), "KL", "Opt")
g<- g + scale_color_discrete(labels=mylabs) + 
    scale_linetype_discrete(labels=mylabs) + 
  ylim(0, .45) + xlim(0, .4) +
  theme(legend.position= c(.5, .9))

ggsave("../../TexDocuments/Figures/comparisonPlot500.pdf", 
       g, width=3.25, height=3.25, units="in")



####
#skewed version
dat = rbind(read.csv("secMom50_skew.csv", header=TRUE), 
            read.csv("bernvar50_skew.csv", header=TRUE), 
            read.csv("pstar50_skew.csv", header=TRUE))

g <- ggplot(aes(x=p1, y=p2, linetype=Type, color=Type), data=dat) + 
  geom_path() + 
  xlab(expression(p[1])) + ylab(expression(p[2])) +
  theme_bw(base_size=12) + 
  theme(legend.title=element_blank(), 
        legend.direction="horizontal",
        text=element_text(family=font))

mylabs <- c(expression(chi^2), "KL", "Opt")
g<- g + scale_color_discrete(labels=mylabs) + 
  scale_linetype_discrete(labels=mylabs) + 
#  ylim(0, .45) + xlim(0, .35) +
  theme(legend.position= c(.5, .9))

ggsave("../../TexDocuments/Figures/comparisonPlot50_skew.pdf", 
       g, width=3.25, height=3.25, units="in")



