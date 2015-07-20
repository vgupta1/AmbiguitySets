###
# Portfolio Exps Plots with Real data for paper
####
library(ggplot2)
library(plyr)
library(dplyr)
library(extrafont) #for latex fonts
library(reshape2)
font = "Times New Roman"


### Messing around with the data
setwd("./Documents/Research/BayesDRO/AmbiguitySets/Experiments/")
mkt = read.csv(file = "12_Industry_Portfolios_clean.csv")
mkt.melt = melt(mkt, id.vars = "Date")

mkt.melt %>% group_by(variable) %>% 
            summarise(avg=mean(value), std=sd(value), var=quantile(value, .1)) %>%
            arrange(avg)

mkt %>% select(-Date) %>% cor()

### Increasing N plots
dat = read.csv("Results/portExp2_72_3.csv")
budget = 3

dat.sum<- dat %>% group_by(N, Method) %>%
  summarise(InPerf=mean(inReturn), 
            outPerf=mean(outReturn), 
            outCVaR=mean(CVaR), 
            exceed=mean(CVaR > budget),
            PerfUp=quantile(outReturn, .9), 
            PerfDown=quantile(outReturn, .1),
            CVaRUp= quantile(CVaR, .90), 
            CVaRDown=quantile(CVaR, .1))

pos = position_dodge(width=100)
my.labs <- c(expression(KL), expression(chi^2), expression(KL[C]), expression(chi[C]^2), expression(SAA))
my.breaks <- c("KL", "Chisq", "KLCov", "ChisqCov", "SAA")

g <- ggplot(aes(x=N, color=Method, y=outPerf, 
           shape=Method), 
       data=dat.sum) + 
    geom_line(aes(linetype=Method), position=pos) + 
  geom_point(position=pos, size=2) + 
    geom_errorbar(aes(ymin=PerfDown, ymax=PerfUp), position=pos) + 
  theme_bw(base_size=10) + 
  theme(legend.title=element_blank()) + 
  ylab("Return (%)")

g <- g + scale_color_discrete(breaks=my.breaks, 
                             labels=my.labs) + 
  scale_shape_discrete(breaks=my.breaks, 
                       labels=my.labs) +
  scale_linetype_discrete(breaks = my.breaks, labels=my.labs)

g<- g + theme(legend.position="top", 
              text=element_text(family=font))

ggsave("../../TexDocuments/Figures/portConvInNRet.pdf", 
       g, width=3.25, height=3.25, units="in")


g<- ggplot(aes(x=N, color=Method, shape=Method, y=outCVaR), data=dat.sum) +
  geom_line(aes(linetype=Method), position=pos) + 
  geom_point(position=pos, size=2) +
  geom_errorbar(aes(ymin=CVaRDown, ymax=CVaRUp), position=pos) + 
  theme_bw(base_size=10) + 
  theme(legend.title=element_blank(), 
        text=element_text(family=font), 
        legend.position="top") + 
  ylab("CVaR (%)") +
  geom_hline(yintercept=3, linetype="dashed")

g <- g + scale_color_discrete(breaks=my.breaks, 
                              labels=my.labs) + 
  scale_shape_discrete(breaks=my.breaks, 
                       labels=my.labs) +
  scale_linetype_discrete(breaks = my.breaks, labels=my.labs)


ggsave("../../TexDocuments/Figures/portConvInNRisk.pdf", 
       g, width=3.25, height=3.25, units="in")

dat.sum %>% filter(N > 750) %>%
  select(N, Method, exceed)

dat.sum %>% filter(Method %in% c("KL", "Chisq"))

########
# Increasing d plots
dat = dat = read.csv("Results/incrDExp_300_3.csv")
dat.sum<- dat %>% group_by(d, Method) %>%
  summarise(InPerf=mean(inReturn), 
            outPerf=mean(outReturn), 
            outCVaR=mean(CVaR), 
            exceed=mean(CVaR > budget),
            PerfUp=quantile(outReturn, .9), 
            PerfDown=quantile(outReturn, .1),
            CVaRUp= quantile(CVaR, .90), 
            CVaRDown=quantile(CVaR, .1))

pos = position_dodge(width=10)
g <- ggplot(aes(x=d, color=Method, y=outPerf, 
                shape=Method), 
            data=dat.sum) + 
  geom_line(aes(linetype=Method), position=pos) + 
  geom_point(position=pos, size=2) + 
  geom_errorbar(aes(ymin=PerfDown, ymax=PerfUp), position=pos) + 
  theme_bw(base_size=10) + 
  theme(legend.title=element_blank(), 
        legend.position="top", 
        text=element_text(family=font)) + 
  ylab("Return (%)")

g <- g + scale_color_discrete(breaks=my.breaks, 
                              labels=my.labs) + 
  scale_shape_discrete(breaks=my.breaks, 
                       labels=my.labs) +
  scale_linetype_discrete(breaks = my.breaks, labels=my.labs)

ggsave("../../TexDocuments/Figures/portConvInDRet.pdf", 
       g, width=3.25, height=3.25, units="in")


g<- ggplot(aes(x=d, color=Method, shape=Method, y=outCVaR), data=dat.sum) +
  geom_line(aes(linetype=Method), position=pos) + 
  geom_point(position=pos, size=2) +
  geom_errorbar(aes(ymin=CVaRDown, ymax=CVaRUp), position=pos) + 
  theme_bw(base_size=10) + 
  theme(legend.title=element_blank(), 
        text=element_text(family=font), 
        legend.position="top") + 
  ylab("CVaR (%)") +
  geom_hline(yintercept=3, linetype="dashed")

g <- g + scale_color_discrete(breaks=my.breaks, 
                              labels=my.labs) + 
  scale_shape_discrete(breaks=my.breaks, 
                       labels=my.labs) +
  scale_linetype_discrete(breaks = my.breaks, labels=my.labs)


ggsave("../../TexDocuments/Figures/portConvInDRisk.pdf", 
       g, width=3.25, height=3.25, units="in")

#################
# analyzing wrong priors
dat = read.csv("Results/randWrongPrior2_20.csv")

dat.sum<- dat %>% group_by(N, Method) %>%
  summarise(InPerf=mean(inReturn), 
            outPerf=mean(outReturn), 
            outCVaR=mean(CVaR), 
            exceed=mean(CVaR > budget),
            PerfUp=quantile(outReturn, .9), 
            PerfDown=quantile(outReturn, .1),
            CVaRUp= quantile(CVaR, .90), 
            CVaRDown=quantile(CVaR, .1))

pos = position_dodge(width=30)
g <- ggplot(aes(x=N, color=Method, y=outPerf, shape=Method), 
            data=dat.sum) + 
  geom_line(aes(linetype=Method), position=pos) + 
  geom_point(position=pos, size=2) + 
  geom_errorbar(aes(ymin=PerfDown, ymax=PerfUp), position=pos) + 
  theme_bw(base_size=10) + 
  theme(legend.title=element_blank(), 
        legend.position="top", 
        text=element_text(family=font)) + 
  ylab("Return (%)")

ggsave("../../TexDocuments/Figures/priorConvRet.pdf", 
       g, width=3.25, height=3.25, units="in")

g<- ggplot(aes(x=N, color=Method, shape=Method, y=outCVaR), data=dat.sum) +
  geom_line(aes(linetype=Method), position=pos) + 
  geom_point(position=pos, size=2) +
  geom_errorbar(aes(ymin=CVaRDown, ymax=CVaRUp), position=pos) + 
  theme_bw(base_size=10) + 
  theme(legend.title=element_blank(), 
        text=element_text(family=font), 
        legend.position="top") + 
  ylab("CVaR (%)") +
  geom_hline(yintercept=3, linetype="dashed")


ggsave("../../TexDocuments/Figures/priorConvRisk.pdf", 
       g, width=3.25, height=3.25, units="in")



#############
# the full sim
#####
myfun <- function(t){quantile(t, .05)}

dat = read.csv("Results/real_sim_36_4.csv")
dat %>% group_by(Method) %>%
  select(-ix, -inReturn) %>%
  summarise_each(funs(mean, sd, myfun))

## dump for julia to do the cvar calc
dcast(dat, ix ~Method, value=outReturn) %>% 
  write.csv("Results/temp_out.csv")

dat %>% filter(Method %in% c("Chisq", "KL", "SAA")) %>%
ggplot(aes(x=outReturn, group=Method, fill=Method)) + 
  geom_histogram(aes(y=..density..)) + 
  facet_grid(Method~.) + 
  ylab("") + 
  xlab("Return (%)") +
  theme_bw(base_size=10) + 
  theme(legend.position="none")

dat %>% filter(Method %in% c("Chisq", "KL", "SAA")) %>%
ggplot(aes(x=inReturn, y=outReturn, color=Method)) + 
  geom_point() + 
  geom_abline()

dat %>% mutate(comp = inReturn > outReturn) %>%
  group_by(Method) %>%
  summarise(mean(comp))

dat %>% mutate(diff = outReturn - inReturn) %>%
  filter(Method %in% c("Chisq", "KL", "SAA")) %>%
  ggplot(aes(x=diff, group=Method, fill=Method)) + 
  geom_histogram() + 
  facet_grid(Method~.) + 
  theme_bw(base_size=10) + 
  theme(legend.position="none")

### plot smoothed returns just to see what's happening
#convert to actual date objects
dat <- dat %>% 
  mutate(Date = mkt$Date[ix-36 + 1062-200], 
         Date2 = as.Date(paste(Date, "01", sep=""), "%Y%m%d"))
  
g<- dat %>% group_by(Method) %>%
  filter(Method %in% c("KL", "Chisq", "SAA")) %>%
  mutate(cumReturn = cumprod(1 + .01 * outReturn)) %>%
  ggplot(aes(x=Date2, y=cumReturn, 
           color=Method, group=Method)) + 
  geom_line(aes(linetype=Method)) + theme_bw(base_size=10) + 
  theme(legend.position = c(.2, .8), 
        text=element_text(family=font), 
        legend.title=element_blank()) +
  ylab("Wealth") + xlab("")

my.labs <- c(expression(chi^2), expression(KL),  expression(SAA))
my.breaks <- c("Chisq", "KL", "SAA")

g <- g + scale_color_discrete(breaks=my.breaks, 
                                labels=my.labs) + 
  scale_linetype_discrete(breaks = my.breaks, labels=my.labs)

ggsave("../../TexDocuments/Figures/returnsHistorical.pdf", 
       g, width=6, height=3.25, units="in")


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
  theme_bw(base_size=10) + 
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
  theme_bw(base_size=10) + 
  theme(legend.title=element_blank(), 
        text=element_text(family=font), 
        legend.position="top") + 
  ylab("CVaR (%)") +
  geom_hline(yintercept=3, linetype="dashed") + 
  xlab(expression(tau[0]))

ggsave("../../TexDocuments/Figures/scalePriorRisk.pdf", 
       g, width=3.25, height=3.25, units="in")


