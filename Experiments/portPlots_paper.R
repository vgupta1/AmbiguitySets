###
# Portfolio Exps Plots with Real data for paper
####
library(ggplot2)
library(plyr)
library(dplyr)
library(extrafont) #for latex fonts
library(reshape2)
font = "Times New Roman"
font = "CM Roman"

### Messing around with the data
setwd("~/Documents/Research/BayesDRO/AmbiguitySets/Experiments/")
mkt = read.csv(file = "Data/12_Industry_Portfolios_clean.csv")

### Increasing N plots
dat = read.csv("Results/portExp2b_72_3.csv")  #uses the empirical mkt
dat = read.csv("Results/portExp3_3.csv")  #uses the cluster mkt
budget = 3

dat.sum<- dat %>% group_by(N, Method) %>%
  summarise(InPerf=mean(inReturn), 
            outPerf=mean(outReturn), 
            outCVaR=mean(CVaR), 
            exceed=mean(CVaR > budget),
            PerfUp=quantile(outReturn, .9), 
            PerfDown=quantile(outReturn, .1),
            CVaRUp= quantile(CVaR, .90), 
            CVaRDown=quantile(CVaR, .1), 
            NumZero = mean( X_Norm <= .001), 
            SDPorts = sd(X_Norm), 
            numRuns = n())

pos = position_dodge(width=50)
my.labs <- c(expression(KL), expression(chi^2), expression(KL[C]), expression(chi[C]^2), expression(SAA), 
             expression(Naive), expression(MinVar))
my.breaks <- c("KL", "Chisq", "KLCov", "ChisqCov", "SAA", "Naive", "MinVar")

g <- ggplot(aes(x=N, color=Method, y=outPerf, 
           shape=Method), data=dat.sum) + 
    geom_line(aes(linetype=Method), position=pos) + 
  geom_point(position=pos, size=2) + 
    geom_errorbar(aes(ymin=PerfDown, ymax=PerfUp), position=pos) + 
  theme_minimal(base_size=10) + 
  theme(legend.title=element_blank(), 
        legend.direction="horizontal", 
        legend.position=c(.5, 1)) + 
  ylab("Return (%)")

g <- g + scale_color_discrete(breaks=my.breaks, 
                             labels=my.labs) + 
  scale_shape_discrete(breaks=my.breaks, 
                       labels=my.labs) +
  scale_linetype_discrete(breaks = my.breaks, labels=my.labs)

# UNCOMMENT for the empirical mkt
# g<- g + theme(legend.direction="horizontal", 
#               legend.position=c(.7, .25), 
#               text=element_text(family=font)) + 
#         guides(col=guide_legend(ncol=2))
# 
# ggsave("../../TexDocuments/Figures/portConvInNRet_a.pdf", 
#        g, width=3.25, height=3.25, units="in")

# UNCOMMNENT for ClUSTER MKT
g<- g + theme(legend.direction="horizontal", 
              legend.position=c(.5, .95), 
              text=element_text(family=font)) + 
        guides(col=guide_legend(nrow=2))
ggsave("../../TexDocuments/Figures/portConvInNRet_b.pdf", 
       g, width=3.25, height=3.25, units="in")



g<- ggplot(aes(x=N, color=Method, shape=Method, y=outCVaR), data=dat.sum) +
  geom_line(aes(linetype=Method), position=pos) + 
  geom_point(position=pos, size=2) +
  geom_errorbar(aes(ymin=CVaRDown, ymax=CVaRUp), position=pos) + 
  theme_minimal(base_size=10) + 
  theme(legend.title=element_blank(), 
        text=element_text(family=font), 
        legend.direction="horizontal", 
        legend.position=c(.7, .25)) + 
  ylab("CVaR (%)") +
  geom_hline(yintercept=3, linetype="dashed")

g <- g + scale_color_discrete(breaks=my.breaks, 
                              labels=my.labs) + 
  scale_shape_discrete(breaks=my.breaks, 
                       labels=my.labs) +
  scale_linetype_discrete(breaks = my.breaks, 
                          labels=my.labs)  

# # UNCOMMNENT for Empirical  MKT
# g <- g + guides(col=guide_legend(ncol=2))
# 
# ggsave("../../TexDocuments/Figures/portConvInNRisk_a.pdf", 
#        g, width=3.25, height=3.25, units="in")

# UNCOMMNENT for ClUSTER MKT
g <- g + theme(legend.position = c(.5, .95))+
      guides(col=guide_legend(nrow=2))
ggsave("../../TexDocuments/Figures/portConvInNRisk_b.pdf", 
       g, width=3.25, height=3.25, units="in")

dat.sum %>% filter(N > 750) %>%
  select(N, Method, exceed)

dat.sum %>% filter(Method %in% c("KL", "Chisq"))

###isolate the number of times something is zero
ggplot(aes(x=N, y=NumZero, fill=Method,group=Method), data=dat.sum) +
  geom_bar(stat = "identity", position="dodge")

dat.sum %>% 
  select(N, Method, NumZero) %>%
  dcast(N ~ Method)

dat.sum %>% 
  select(N, Method, SDPorts) %>%
  dcast(N ~ Method)


d1 = dat.sum %>% filter(Method %in% c("ChisqCov", "KLCov")) %>%
    select(N, Method, NumZero) %>%
  dcast(N ~ Method)
d2 = dat.sum %>% 
  select(N, Method, SDPorts) %>%
  dcast(N ~ Method)


########
# Increasing d plots
dat = read.csv("Results/incrDExp_a_300_3.csv")
dat = read.csv("Results/incrDExp_a_750_3.csv")

dat.sum<- dat %>% group_by(d, Method) %>%
  summarise(InPerf=mean(inReturn), 
            outPerf=mean(outReturn), 
            outCVaR=mean(CVaR), 
            exceed=mean(CVaR > budget),
            PerfUp=quantile(outReturn, .9), 
            PerfDown=quantile(outReturn, .1),
            CVaRUp= quantile(CVaR, .90), 
            CVaRDown=quantile(CVaR, .1))

convInDPlot <- function(dat, yval, yvaldown, yvalup,  ylabel){
  pos = position_dodge(width=10)
  g <- ggplot(aes_string(x="d", color="Method", y=yval, 
                  shape="Method", ymin=yvaldown, ymax=yvalup), 
              data=dat) + 
    geom_line(aes(linetype=Method), position=pos) + 
    geom_point(position=pos, size=2) + 
    geom_errorbar(position=pos) + 
    theme_minimal(base_size=10) + 
    theme(legend.title=element_blank(), 
          legend.position="top", 
          text=element_text(family=font)) + 
    ylab(ylabel) 
  g <- g + scale_color_discrete(breaks=my.breaks, 
                                labels=my.labs) + 
    scale_shape_discrete(breaks=my.breaks, 
                         labels=my.labs) +
    scale_linetype_discrete(breaks = my.breaks, labels=my.labs)
  return(g)
}

#first do one with everyone for the appendix.  These an be big
g <- convInDPlot(dat.sum, "outPerf", "PerfDown", "PerfUp", "Return (%)")
ggsave("../../TexDocuments/Figures/portConvInDRet_300_full.pdf", 
       g, width=6.5, height=3.25, units="in")

g <- convInDPlot(dat.sum, "outCVaR", "CVaRDown", "CVaRUp", "CVaR (%)")
ggsave("../../TexDocuments/Figures/portConvInDRisk_300_full.pdf", 
       g, width=6.5, height=3.25, units="in")

#now limit down to interesting things for main text.  Should be smaller.
g <- dat.sum %>% filter(!Method %in% c("Naive", "MinVar")) %>%
  convInDPlot(., "outPerf", "PerfDown", "PerfUp", "Return (%)")
g <- g + theme(legend.position = c(.7, .8)) + 
    guides(col=guide_legend(ncol=2))
ggsave("../../TexDocuments/Figures/portConvInDRet_300.pdf", 
       g, width=3.25, height=3.25, units="in")

g <- dat.sum %>% filter(!Method %in% c("Naive", "MinVar")) %>%
  convInDPlot(., "outCVaR", "CVaRDown", "CVaRUp", "CVaR (%)")
g <- g + theme(legend.position = c(.75, .3)) + 
  guides(col=guide_legend(ncol=2))
ggsave("../../TexDocuments/Figures/portConvInDRisk_300.pdf", 
       g, width=3.25, height=3.25, units="in")

### Finally do a big one with 700 points for the appendix
#first do one with everyone for the appendix.  These an be big
g <- convInDPlot(dat.sum, "outPerf", "PerfDown", "PerfUp", "Return (%)")
ggsave("../../TexDocuments/Figures/portConvInDRet_700_full.pdf", 
       g, width=6.5, height=3.25, units="in")

g <- convInDPlot(dat.sum, "outCVaR", "CVaRDown", "CVaRUp", "CVaR (%)")
ggsave("../../TexDocuments/Figures/portConvInDRisk_700_full.pdf", 
       g, width=6.5, height=3.25, units="in")

#############
# the full sim
#####
myfun <- function(t){-quantile(t, .1)}

cvar <- function(t){
    var = quantile(t, .1)
    -mean(t[ t <= var ])}

dat = read.csv("Results/real_sim_36_6.csv")

#Add actual date objects for plotting
dat$dt = mkt$Date[dat$ix-36 + 1062-200]
dat <- dat %>% 
  mutate(dt2 = as.Date(paste(dt, "01", sep=""), "%Y%m%d"))

### Comptue the rolling CVaR by method
dat <- dat %>% group_by(Method) %>%
      mutate( CVaR72 = rollapply(outReturn, width=72, cvar, 
                             fill=NA, align="right"))

## First compute the summary table for the document:
dat %>% group_by(Method) %>% #filter(!Method %in% c("MinVarUnsc", "NaiveUnsc")) %>%
  summarise(AvgReturn = mean(outReturn), 
            StdDev = sd(outReturn), 
            VaR = -quantile(outReturn, .1), 
            CVaR = cvar(outReturn),
            RetFeas = mean(outReturn > inReturn), 
            Ratio = AvgReturn/CVaR) %>%
  arrange(Ratio)

wealthPlot <- function(d){
  g<- d %>% group_by(Method) %>%
    filter(!Method %in% c("NaiveUnsc", "MinVarUnsc")) %>%
    mutate(cumReturn = cumprod(1 + .01 * outReturn)) %>%
    ggplot(aes(x=dt2, y=cumReturn, 
               color=Method, group=Method)) + 
    geom_line(aes(linetype=Method)) + theme_minimal(base_size=10) + 
    theme(legend.position = c(.2, .8), 
          text=element_text(family=font), 
          legend.title=element_blank(),
          legend.direction="horizontal") +
    ylab("Wealth") + xlab("")  
  g <- g + scale_color_discrete(breaks=my.breaks, 
                                labels=my.labs) + 
    scale_linetype_discrete(breaks = my.breaks, labels=my.labs) + 
    guides(col=guide_legend(ncol=2))
}

# Plot the wealth process
#SAVE a big one for teh appendix with everybody
#Since it's color, you can just natural labeling
g <- wealthPlot(dat)
ggsave("../../TexDocuments/Figures/returnsHistorical_Full.pdf", 
       g, width=6.5, height=4, units="in")


#Now do one without MinVar, without KL, and without KLC, 
g <- dat %>% filter(!Method %in% c("KL", "MinVar", "KLCov")) %>%
  wealthPlot()
#can't assume it's color, so add label in PPT
g <- g +theme(legend.position="none")
ggsave("../../TexDocuments/Figures/returnsHistorical_nolabels.pdf", 
       g, width=3.25, height=3.25, units="in")

##Plot rolling CVAR
cvarPlot <- function(dt){
  dt %>% filter(!is.na(CVaR72), 
             !Method %in% c("NaiveUnsc", "MinVarUnsc")) %>%
    ggplot(aes(x=dt2, y=CVaR72, group=Method, 
               color=Method, linetype=Method), 
           data=.) + 
    scale_linetype_discrete(breaks = my.breaks, labels=my.labs) +
    scale_color_discrete(breaks = my.breaks, labels=my.labs) +
    geom_line() +
    geom_hline(yintercept=6, linetype="dotted") + 
    xlab("") + ylab("Rolling CVaR") + 
    theme_minimal(base_size=10) + 
    theme(legend.title=element_blank(), 
          text=element_text(family=font), 
          legend.direction="horizontal", 
          legend.position =c(.5, 1))
}
g <- cvarPlot(dat)

#SAVE a big one for teh appendix with everybody
ggsave("../../TexDocuments/Figures/rollingCVaR_Full.pdf", 
       g, width=6.5, height=4, units="in")


##Make a smaller one.  No labels.  ADded in ppt. 
g <- dat %>% filter(!Method %in% c("KL", "MinVar", "KLCov")) %>%
      cvarPlot()
g <- g +theme(legend.position="none")
ggsave("../../TexDocuments/Figures/CVaRHistorical_nolabels.pdf", 
       g, width=3.25, height=3.25, units="in")

