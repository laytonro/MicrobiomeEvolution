library(dplyr)
library(ggplot2)
library(tidyverse)
library("stringr")  


setwd("/Users/laytonrosenfeld/Dropbox/layton_project/Ecological_structure")
d <- read.csv(file = 'indexes.csv')


for (i in 1:nrow(d)){
  a = d$mod_rep_hybrid[i]
  d$mod_rep_hybrid2[i]=str_sub(a, -3, -1)
}

for (i in 1:nrow(d)){
  a = d$mod_rep_hybrid2[i]
  if (a=="mod"){
    d$mod_rep_hybrid2[i]="mod_only"
  }
  
  if (a == "ins") {
    d$mod_rep_hybrid2[i]="twins"
  }
  else if (a == "ges") {
    d$mod_rep_hybrid2[i]="no_within_species_changes"
  }
}


d$mod_rep_binned = d$mod_rep_hybrid

for (i in 1:nrow(d)){
  a = d$mod_rep_hybrid[i]
  if (a=="3_mod"){
    d$mod_rep_binned[i]="2_mod"
  }
  if (a=="3_rep"){
    d$mod_rep_binned[i]="2_rep"
  }
}



l <- c("2_mod", "1_mod", "no_within_species_changes", "1_rep", "2_rep", "twins",
        "unrelated_hosts")

colors <- c("#d6e5ff", "#0d64fc", "#c2bebe", "#ff7c70",
             "#ff2a17", "#77d195", "#e8ff75")


final_labels <- c("M ≥ 2\nR = 0", "M = 1\nR = 0", "M = 0\nR = 0", 
                   "R = 1\nM ≥ 0", "R ≥ 2\nM ≥ 0",  "twins", "unrelated\nhosts"
)


d_binned <- d %>% select(jensen_shannon, mod_rep_binned)
h_shuffled_pairs <- read.csv(file = 'null_hmp_indexes.csv')
h_shuffled_pairs$mod_rep_binned = "unrelated_hosts"

d_binned <- rbind(d_binned, h_shuffled_pairs %>% select(jensen_shannon, mod_rep_binned))


#main boxplot with bins
ggplot(d_binned, aes(x=factor(mod_rep_binned, level=l), y=jensen_shannon))+ 
  geom_boxplot(aes(group=mod_rep_binned, fill= factor(mod_rep_binned, level=l)),
               outlier.shape=NA) + 
  labs(y="Jensen Shannon Index", x = "evolutionary events ~ cohort") + 
  scale_fill_manual(values = colors, name = "", labels = final_labels) + 
  scale_x_discrete(labels = final_labels) + 
  theme(
    legend.spacing.y = unit(1.0, 'cm'),
    legend.key.size = unit(2.0, 'lines')
  ) + 
  labs(subtitle = "R = Strain replacement\nM = Evolutionary modification")


#boxplot of replacements vs jensen_shannon by cohort
ggplot(d, aes(x=replacements, y=jensen_shannon, color = cohort))+ 
  geom_boxplot(aes(group=replacements),outlier.shape=NA) + 
  geom_jitter(height=0, width=.4) + scale_x_continuous(breaks=0:5) +
  facet_wrap( ~ cohort,scale="free_x" )

#boxplot of modificatins vs jensen_shannon by cohort

ggplot(d, aes(x=modifications, y=jensen_shannon, color = cohort))+ 
  geom_boxplot(aes(group=modifications), outlier.shape=NA) + 
  geom_jitter(height=0, width=.4) + scale_x_continuous(breaks=0:5) +
  facet_wrap(~ cohort,scale="free_x" )


#Regression on replacements x JS (for hmp samples)

js_model = lm(jensen_shannon ~ replacements, data = d%>%filter(cohort=="hmp"))
js_model
summary(js_model)

confint(js_model, level=.95)


ggplot(d, aes(x=replacements, y=jensen_shannon)) + 
  geom_point(color='#2980B9', size = 4) + 
  geom_smooth(method="lm", color='#2C3E50')
