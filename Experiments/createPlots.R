###  
# Creates the Conv in N Plots for Paper
###
library(ggplot2)
library(plyr)
library(dplyr)
library(extrafont) #for latex fonts
library(reshape2)

#VG Kill this before submitting?
setwd("/Users/vishalgupta/Documents/Research/BayesDRO/AmbiguitySets/Experiments/Results/")
font = "CM Roman"
font = "Times New Roman"

#########
#first create a plot without asymptotic values
dat = read.csv(file="new_conv_N_plot.csv", header=TRUE)
dat = read.csv(file="../temp_vals.csv", header=TRUE)

g = ggplot(aes(x=N, y=Ratio, color=Type), 
       data=filter(dat, N>0 & N < 1e4)) + 
  geom_line(aes(linetype=factor(Direction))) + 
  geom_point(aes(shape=Type), size=2) +
  scale_x_log10(breaks=c(10, 100, 1000, 10000)) +
  scale_linetype_discrete(guide=FALSE) + 
  theme_minimal(base_size=12) + 
  theme(legend.title=element_blank(), 
        legend.position=c(.7, .8), 
        text=element_text(family=font), 
        legend.direction="horizontal") + 
  ylab("")

my.labs <- c(expression(KL), expression(chi^2), expression(KL[C]), expression(chi[C]^2))
my.breaks <- c("KL", "ChiSq", "KL_C", "ChiSq_C")


g<- g + scale_color_discrete(breaks=my.breaks, 
                     labels=my.labs) + 
  scale_shape_discrete(breaks=my.breaks, 
                       labels=my.labs)

g<- g + guides(color=guide_legend(ncol=2))

ggsave("../../../TexDocuments/Figures/convInN.pdf", 
       g, width=3.25, height=3.25, units="in")

##########
# The epsilon plot
######
setwd("/Users/vishalgupta/Documents/Research/BayesDRO/AmbiguitySets/Experiments/Results/")
dat_ = read.csv(file="new_eps_plot.csv", header=TRUE)
dat = melt(dat_, id.vars=c("eps"), variable.name="Type")

g <- ggplot(aes(x=eps, y=value, color=Type), 
           data=dat) + 
  geom_line(aes(linetype=Type)) + 
  geom_point(aes(shape=Type), size=2) +
  xlab(expression(epsilon)) + ylab("") +
  scale_linetype_discrete(guide=FALSE) + 
  theme_minimal(base_size=12) + 
  theme(legend.title=element_blank(), 
        legend.position=c(.5, .95), 
        text=element_text(family=font), 
        legend.direction="horizontal")

my.labs <- c(expression(KL), expression(chi^2), 
             expression(phi-"Div, d=5"), 
             expression(phi-"Div, d=10"), 
             expression(phi-"Div, d=20") )
my.breaks <- c("KL", "ChiSq", "KL_C5", "KL_C10", "KL_C20")
g<- 
  g + scale_color_discrete(breaks=my.breaks, 
                             labels=my.labs) + 
  scale_shape_discrete(breaks=my.breaks, 
                       labels=my.labs) + 
  guides(col=guide_legend(nrow=2))

ggsave("../../../TexDocuments/Figures/eps_plot.pdf", 
       g, width=3.25, height=3.25, units="in")



##########
# Plots comparing epsilons
dat = read.csv("comp_eps_plot.csv", header=TRUE)
dat$d = factor(dat$d)

#d = 3 is the same as the Good KL
dat$d = revalue(dat$d, c("3"="KL"))

g<- ggplot(aes(x=Cov_eps, y=Real_eps, 
           group=d, color=d, shape=d), 
       data=dat) + 
  geom_point(size=3) + 
  geom_line() 

g<- g + theme_minimal(base_size=12) + 
  xlab(expression(paste("Desired ",epsilon))) + 
  ylab(expression(paste("Achieved ", epsilon))) + 
  theme(legend.position=c(.15, .8), 
        legend.title=element_blank(), 
        text=element_text(family=font))

my.labs <- c(expression(KL), 
             expression(chi^2),
             expression("d=5"), 
             expression("d=10") )
my.breaks <- c("KL", "Chisq", "5", "10")

g<- 
  g + scale_color_discrete(breaks=my.breaks, 
                           labels=my.labs) + 
  scale_shape_discrete(breaks=my.breaks, 
                       labels=my.labs) 
  
g<- g + geom_abline(linetype="dotted")

ggsave("../../../TexDocuments/Figures/comp_eps_kl.pdf", 
       g, width=3.25, height=3.25, units="in")


