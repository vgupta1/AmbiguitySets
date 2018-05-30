library(ggplot2)
library(plyr)
library(dplyr)
library(extrafont) #for latex fonts
library(reshape2)
font = "Times"

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
  xlab(expression(hat(tau)[0]))

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
  xlab(expression(hat(tau)[0]))

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
dat = read.csv("Results/random_wrong_priors2_72_3.csv")
dat = read.csv("Results/random_wrong_priors_1000_72_3.csv")

head(dat)

##test if you can collapse things by prior
t <- dat %>% group_by(iPrior, N, strength, Method) %>%
  summarize(Return = mean(outReturn), 
            stdErrReturn = sd(outReturn)/ sqrt(n()), 
            CVaR= mean(outCVaR), 
            stdErrCVar = sd(outCVaR)/ sqrt(n()),
            Num = n())

#std errors for return and CvAr are smaller than digits shown?
max(t$stdErrReturn)
max(t$stdErrCVar)


### summarize data
d1 <- dat %>% group_by(iPrior, Method) %>%
  summarize(Return = mean(outReturn),
            CVaR= mean(outCVaR), 
            N = n())
d2 <- dat %>% select(iPrior, Method, MSE, portMSE, AssetMSE)
d <- join(d1, d2, type="left")


d %>%
  ggplot(aes(x=AssetMSE, y=CVaR, group=1), data=.) + 
  geom_smooth() + geom_point()

## Plots look like garbage
## Just make a table.
dat.sum<- dat %>% group_by(Method, strength) %>%
  summarize(avgReturn = mean(outReturn), 
            lowReturn = quantile(outReturn, .1), 
            highReturn = quantile(outReturn, .9),
            minReturn = min(outReturn), 
            maxReturn = max(outReturn), 
            avgCVaR = mean(outCVaR), 
            lowCVaR = quantile(outCVaR, .1),
            highCVaR = quantile(outCVaR, .9),            
            minCVaR = min(outCVaR), 
            maxCVaR = max(outCVaR))  %>%
            mutate(strength = 100*strength)
dat.sum %>% select(-minReturn, -maxReturn, -minCVaR, -maxCVaR)

