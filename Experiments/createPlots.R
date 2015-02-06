###  
# Creates the conv in N Plots for
###
library(ggplot2)
library(dplyr)
library(extrafont) #for latex fonts
library(reshape2)

setwd("/Users/vishalgupta/Documents/Research/BayesDRO/AmbiguitySets/Experiments/")
#font = "CM Roman"
font = "Times New Roman"

#########
#first create a plot without asymptotic values
dat = read.csv(file="plotConvInN.csv", header=TRUE)

g = ggplot(aes(x=N, y=Ratio, color=Type), 
       data=filter(dat, N>0 & N < 1e3)) + 
  geom_line(aes(linetype=factor(Direction))) + 
  geom_point(aes(shape=Type), size=2) +
  scale_x_log10() + 
  scale_linetype_discrete(guide=FALSE) + 
  theme_bw(base_size=12) + 
  theme(legend.title=element_blank(), 
        legend.position=c(.7, .8), 
        text=element_text(family=font), 
        legend.direction="horizontal")

my.labs <- c(expression(chi^2), expression(chi[C]^2), expression(KL), expression(KL[C]))
g<- g + scale_color_discrete(breaks=c("ChiSq", "ChiSq_C", "KL", "KL_C"), 
                     labels=my.labs) + 
  scale_shape_discrete(breaks=c("ChiSq", "ChiSq_C", "KL", "KL_C"), 
                       labels=my.labs)

g<- g + guides(color=guide_legend(ncol=2))

ggsave("../../TexDocuments/Figures/convInN.pdf", 
       g, width=3.25, height=3.25, units="in")

##########
# The epsilon plot
######
setwd("/Users/vishalgupta/Documents/Research/BayesDRO/AmbiguitySets/Experiments/")
dat_ = read.csv(file="eps_plot.csv", header=TRUE)
dat = melt(dat, id.vars=c("eps"), variable.name="Type")

g <- ggplot(aes(x=eps, y=value, color=Type), 
           data=dat) + 
  geom_line(aes(linetype=Type)) + 
  geom_point(aes(shape=Type), size=2) +
  xlab(expression(epsilon)) + ylab("") +
  scale_linetype_discrete(guide=FALSE) + 
  theme_bw(base_size=12) + 
  theme(legend.title=element_blank(), 
        legend.position=c(.7, .8), 
        text=element_text(family=font))

my.labs <- c(expression(KL), expression(chi^2))
g<- g + scale_color_discrete(breaks=c("KL", "ChiSq"), 
                             labels=my.labs) + 
  scale_shape_discrete(breaks=c("KL", "ChiSq"), 
                       labels=my.labs)

ggsave("../../TexDocuments/Figures/eps_plot.pdf", 
       g, width=3.25, height=3.25, units="in")



##########
# Plots comparing epsilons
dat = read.csv("comp_eps_kl.csv", header=TRUE)
dat$d = factor(dat$d)

g<- ggplot(aes(x=Cov_eps, y=Real_eps, 
           group=d, color=d, linetype=d, shape=d), 
       data=dat) + 
  geom_point() + 
  geom_line() 

g<- g + theme_bw(base_size=12) + 
  xlab(expression(paste("Desired ",epsilon))) + 
  ylab(expression(paste("Achieved ", epsilon))) + 
  theme(legend.position=c(.2, .8), 
        legend.title=element_blank(), 
        text=element_text(family=font))

ggsave("../../TexDocuments/Figures/comp_eps_kl.pdf", 
       g, width=3.25, height=3.25, units="in")


