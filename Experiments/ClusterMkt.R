####
# Building Mkt Clusters for analysis
###
library(ggplot2)
library(plyr)
library(dplyr)
library(extrafont) #for latex fonts
library(reshape2)

setwd("~/Documents/Research/BayesDRO/AmbiguitySets/Experiments/")
mkt = read.csv(file = "Data/12_Industry_Portfolios_clean.csv")

head(mkt)
T = 120  #Use 10 years of data
mkt_small = mkt[(1062-T+1):1062, 2:13]
t = hclust(dist(mkt_small))
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
             ylab("Probability")  + 
             theme_minimal(base_size=10) + 
             theme(text=element_text(family=font)  )
           
ggsave("../../TexDocuments/Figures/histOfClusterProbs.pdf", 
                  g, width=3.25, height=3.25, units="in")
           
           
 d1 = mkt_small %>% group_by(cluster) %>%
  summarise_each(funs(mean))
           d2 = mkt_small %>% group_by(cluster) %>% summarise(Prob=n()/T)
           d1 = join(d1, d2)
           write.csv(d1, file="Data/ClusterScenarios.csv")
           

#what do certain clusters look like?
d1.melt = melt(d1, id.vars="cluster", variable.name="Index")
d1.melt$cluster = as.character(d1.melt$cluster)

#Pick out the smallest probs
small_clusters = d1.melt %>% filter(Index == "Prob") %>%
  filter(value < .0084) %>% select(cluster)

large_clusters = d1.melt %>% filter(Index == "Prob") %>%
  filter(value >= .1) %>% select(cluster)

d1.melt %>% filter(Index != "Prob") %>%
  filter(cluster  %in% c(small_clusters$cluster)) %>%
  ggplot(aes(x=Index, y=value, color=cluster, group=cluster), data=.) + 
  geom_point() + geom_line()


d1.melt %>% filter(Index != "Prob") %>%
  filter(cluster  %in% c(large_clusters$cluster)) %>%
  ggplot(aes(x=Index, y=value, color=cluster, group=cluster), data=.) + 
  geom_point() + geom_line()

#Graph a subset of the small clusters for paper
g <- d1.melt %>% filter(Index != "Prob") %>%
  filter(cluster  %in% c("23", "27", "25", "35")) %>%
  ggplot(aes(x=Index, y=value, color=cluster, group=cluster, 
             linetype=cluster, shape=cluster), data=.) + 
  geom_point() + geom_line() + 
  theme_minimal(base_size=10) + 
  ylab("Return (%)") + 
  xlab("") +
  theme(legend.position="none", 
        axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("../../TexDocuments/Figures/SampleClusters.pdf", 
       g, width=3.25, height=3.25, units="in")


