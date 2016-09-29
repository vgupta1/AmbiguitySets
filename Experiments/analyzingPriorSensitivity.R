library(ggplot2)
library(plyr)
library(dplyr)
library(extrafont) #for latex fonts
library(reshape2)
font = "Times New Roman"

setwd("~/Documents/Research/BayesDRO/AmbiguitySets/Experiments/")

#############
# The increasingly incorrect prior
##############
dat = read.csv("Results/randWrongPriorScale_300_72.csv", header=TRUE)
budget = 3
dat.sum<- dat %>% group_by(Scale) %>%
  summarise(InPerf=mean(inReturn), 
            outPerf=mean(outReturn), 
            outCVaR=mean(CVaR), 
            exceed=mean(CVaR > budget),
            PerfUp=quantile(outReturn, .9), 
            PerfDown=quantile(outReturn, .1),
            CVaRUp= quantile(CVaR, .90), 
            CVaRDown=quantile(CVaR, .1))

g <- ggplot(aes(x=Scale, y=outPerf), 
            data=dat.sum) + 
  geom_line() + 
  geom_point(size=2) + 
  geom_errorbar(aes(ymin=PerfDown, ymax=PerfUp)) + 
  theme_minimal(base_size=10) + 
  theme(legend.title=element_blank(), 
        legend.position="top", 
        text=element_text(family=font)) + 
  ylab("Return (%)") + 
  xlab(expression(tau[0]))

ggsave("../../TexDocuments/Figures/scalePriorReturn.pdf", 
       g, width=3.25, height=3.25, units="in")

g<- ggplot(aes(x=Scale, y=outCVaR), data=dat.sum) +
  geom_line() + 
  geom_point(size=2) +
  geom_errorbar(aes(ymin=CVaRDown, ymax=CVaRUp)) + 
  theme_minimal(base_size=10) + 
  theme(legend.title=element_blank(), 
        text=element_text(family=font), 
        legend.position="top") + 
  ylab("CVaR (%)") +
  geom_hline(yintercept=3, linetype="dashed") + 
  xlab(expression(tau[0]))

ggsave("../../TexDocuments/Figures/scalePriorRisk.pdf", 
       g, width=3.25, height=3.25, units="in")

####
# Attenuation in N
####
dat = read.csv(file="Results/wrongPrior_IncrN_72_3.csv")
dat.sum<- dat %>% group_by(N, Method) %>%
  summarise(InPerf=mean(inReturn), 
            outPerf=mean(outReturn), 
            outCVaR=mean(CVaR), 
            exceed=mean(CVaR > budget),
            PerfUp=quantile(outReturn, .9), 
            PerfDown=quantile(outReturn, .1),
            CVaRUp= quantile(CVaR, .90), 
            CVaRDown=quantile(CVaR, .1))

g<- ggplot(aes(x=N, outPerf, ymin=PerfDown, ymax=PerfUp), 
       data = dat.sum) + 
  geom_point() + geom_errorbar() + geom_line() +
  theme_minimal(base_size=10) + 
  xlab("N") + ylab("Return (%)")
ggsave("../../TexDocuments/Figures/scalePriorReturnIncreasingN.pdf", 
       g, width=3.25, height=3.25, units="in")


g<- ggplot(aes(x=N, outCVaR, ymin=CVaRDown, ymax=CVaRUp), 
       data = dat.sum) + 
  geom_point() + geom_errorbar() + geom_line() +
  theme_minimal(base_size=10) + 
  xlab("N") + ylab("CVaR (%)") + 
  geom_hline(yintercept=3, linetype= "dotted")
ggsave("../../TexDocuments/Figures/scalePriorRiskIncreasingN.pdf", 
       g, width=3.25, height=3.25, units="in")



####
# Randomly Generated Priors and Performance
####
dat = read.csv("Results/random_wrong_priors_72_3.csv", header=TRUE)

head(dat)

filter(dat, Method == "ChiSq", N == 300) %>%
  ggplot(aes(x=sqrt(MSE), y=outReturn, group=factor(N), color=factor(N)), data=.) + 
  geom_smooth() + geom_point()


dat.sum <- dat %>% group_by(iPrior, N, relProb, MSE, Method) %>%
  summarize(Return = mean(outReturn), CVaR =mean(outCVaR))

filter(dat.sum, Method == "ChiSq", N == 300) %>%
  ggplot(aes(x=MSE, y=Return, group=factor(N), color=factor(N)), data=.) + 
  geom_smooth() + geom_point()

