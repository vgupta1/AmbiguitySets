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
mkt.melt = melt(mkt, id.vars = "Date")

mkt.melt %>% group_by(variable) %>% 
            summarise(avg=mean(value), std=sd(value), var=quantile(value, .1)) %>%
            arrange(avg)

mkt %>% select(-Date) %>% cor()

### Increasing N plots
dat = read.csv("Results/portExp2b_72_3.csv")
dat = read.csv("Results/portExp3_3.csv")
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

# g<- g + theme(legend.direction="horizontal", 
#               legend.position=c(.7, .25), 
#               text=element_text(family=font)) + 
#         guides(col=guide_legend(ncol=2))
# 
# ggsave("../../TexDocuments/Figures/portConvInNRet_a.pdf", 
#        g, width=3.25, height=3.25, units="in")

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

# g <- g + guides(col=guide_legend(ncol=2))
# 
# ggsave("../../TexDocuments/Figures/portConvInNRisk_a.pdf", 
#        g, width=3.25, height=3.25, units="in")

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

convInDPlot <- function(dat, yval, yvaldown, yvalup){
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
    ylab("Return (%)") 
  g <- g + scale_color_discrete(breaks=my.breaks, 
                                labels=my.labs) + 
    scale_shape_discrete(breaks=my.breaks, 
                         labels=my.labs) +
    scale_linetype_discrete(breaks = my.breaks, labels=my.labs)
  return(g)
}

#first do one with everyone for the appendix.  These an be big
g <- convInDPlot(dat.sum, "outPerf", "PerfDown", "PerfUp")
ggsave("../../TexDocuments/Figures/portConvInDRet_300_full.pdf", 
       g, width=6.5, height=3.25, units="in")

g <- convInDPlot(dat.sum, "outCVaR", "CVaRDown", "CVaRUp")
ggsave("../../TexDocuments/Figures/portConvInDRisk_300_full.pdf", 
       g, width=6.5, height=3.25, units="in")

#now limit down to interesting things for main text.  Should be smaller.
g <- dat.sum %>% filter(!Method %in% c("Naive", "MinVar")) %>%
  convInDPlot(., "outPerf", "PerfDown", "PerfUp")
g <- g + theme(legend.position = c(.7, .8)) + 
    guides(col=guide_legend(ncol=2))
ggsave("../../TexDocuments/Figures/portConvInDRet_300.pdf", 
       g, width=3.25, height=3.25, units="in")

g <- dat.sum %>% filter(!Method %in% c("Naive", "MinVar")) %>%
  convInDPlot(., "outCVaR", "CVaRDown", "CVaRUp")
g <- g + theme(legend.position = c(.8, .3)) + 
  guides(col=guide_legend(ncol=2))
ggsave("../../TexDocuments/Figures/portConvInDRisk_300.pdf", 
       g, width=3.25, height=3.25, units="in")






##now do another dropping some of the methods

ggsave("../../TexDocuments/Figures/portConvInDRet_700_a.pdf", 
       g, width=3.25, height=3.25, units="in")




g<- ggplot(aes(x=d, color=Method, shape=Method, y=outCVaR), data=dat.sum) +
  geom_line(aes(linetype=Method), position=pos) + 
  geom_point(position=pos, size=2) +
  geom_errorbar(aes(ymin=CVaRDown, ymax=CVaRUp), position=pos) + 
  theme_minimal(base_size=10) + 
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


ggsave("../../TexDocuments/Figures/portConvInDRisk_700_a.pdf", 
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
  theme_minimal(base_size=10) + 
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
  theme_minimal(base_size=10) + 
  theme(legend.title=element_blank(), 
        text=element_text(family=font), 
        legend.position="top") + 
  ylab("CVaR (%)") +
  geom_hline(yintercept=3, linetype="dashed")


ggsave("../../TexDocuments/Figures/priorConvRisk.pdf", 
       g, width=3.25, height=3.25, units="in")


#########
# The new wrong-prior test
########
dat = read.csv("Results/randWrongPrior_100_72_8675309.csv")
head(dat)

##contour plot doesn't like that dist isn't rounded
#So fit a response surface.
prior.loess <- loess( CVaR ~ Dist + Scale + Dist*Scale, data = dat)
dat.prior <- expand.grid(list(Dist=seq(.05, .09, .0025), Scale=seq(.1, 3, .01)))
dat.prior <- dat.prior %>% mutate(CVaRSooth = as.numeric(predict(prior.loess, dat.prior)))
prior.loess2 <- loess( outReturn ~ Dist + Scale + Dist*Scale, data = dat)
dat.prior <- dat.prior %>% mutate(ReturnSmooth = as.numeric(predict(prior.loess2, dat.prior)))

#Figure out how to get labels on the contours...
##Something seems wrong with the returns....
ggplot(aes(x=Dist, y=Scale, z=ReturnSmooth), data=dat.prior) +
  geom_tile(aes(fill=ReturnSmooth)) + stat_contour()
  
##where do we get violations?
dat %>% mutate(Viols = CVaR > 3) %>%
  group_by(Scale) %>%
  summarize(probViol = mean(Viols) ) %>%
  qplot(x=Scale, y=probViol, data=.) + geom_point()



#############
# the full sim
#####
myfun <- function(t){-quantile(t, .05)}

dat = read.csv("Results/real_sim_36_6.csv")
dat %>% filter(!Method %in% c("NaiveUnsc", "MinVarUnsc")) %>%
  group_by(Method) %>%
  select(-ix, -inReturn) %>%
  summarise_each(funs(mean, sd, myfun)) %>%
  mutate(ratio = mean/myfun) %>%
  arrange(ratio)

## dump for julia to do the cvar calc
dcast(dat, ix ~Method, value.var="outReturn") %>% 
  write.csv("Results/temp_out.csv")

dat %>% filter(Method %in% c("Chisq", "KL", "Naive", "MinVar")) %>%
ggplot(aes(x=outReturn, group=Method, color=Method, linetype=Method)) + 
  geom_density() + 
  ylab("") + 
  xlab("Return (%)") +
  theme_minimal(base_size=10) + 
  theme(legend.position=c(.3, .8), 
        legend.title=element_blank())

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
  theme_minimal(base_size=10) + 
  theme(legend.position="none")

### plot smoothed returns just to see what's happening
#convert to actual date objects
dat <- dat %>% 
  mutate(Date = mkt$Date[ix-36 + 1062-200], 
         Date2 = as.Date(paste(Date, "01", sep=""), "%Y%m%d"))
  
g<- dat %>% group_by(Method) %>%
  filter(!Method %in% c("NaiveUnsc", "MinVarUnsc")) %>%
  mutate(cumReturn = cumprod(1 + .01 * outReturn)) %>%
  ggplot(aes(x=Date2, y=cumReturn, 
           color=Method, group=Method)) + 
  geom_line(aes(linetype=Method)) + theme_minimal(base_size=10) + 
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
# clustering analysis
###
head(mkt)
T = 120
mkt_small = mkt[(1062-T+1):1062, 2:13])
t = hclust(dist(mkt_small)
plot(t)
t_cut <- cutree(t, 36)

mkt_small$cluster = t_cut

g <- mkt_small %>% group_by(cluster) %>%
  summarise(n = n()/T) %>%
  arrange(desc(n)) %>%
  ggplot(aes(x=1:36, y=n), data=.) + 
  geom_bar(stat="identity") +
  geom_hline(yintercept=1/36, linetype="dotted") + 
  xlab("Scenario") + 
  ylab("Prob")  + 
  theme_minimal(base_size=10) + 
  theme(text=element_text(family=font)  )

ggsave("../../TexDocuments/Figures/histOfClusterProbs.pdf", 
       g, width=3.25, height=3.25, units="in")


d1 = mkt_small %>% group_by(cluster) %>%
  summarise_each(funs(mean), -Date)
d2 = mkt_small %>% group_by(cluster) %>% summarise(Prob=n()/T)
d1 = join(d1, d2)

write.csv(d1, file="Data/ClusterScenarios.csv")

