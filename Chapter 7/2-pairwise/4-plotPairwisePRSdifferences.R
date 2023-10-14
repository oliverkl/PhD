# R code to produce result figures

library(tidyverse)
library(patchwork)
library(ggpubr)
library(ggridges)

pairwise_database <- readRDS("pairwise_database_optimised.rds")

### Boxplots
##### concordant v discordant

sample_size = pairwise_database %>% 
  group_by(sev_diff) %>% 
  summarize(num=n())

sample_size

mean.degrees <- pairwise_database %>% 
  group_by(sev_diff) %>% 
  summarize(m=mean(prs_diff))

mean.degrees

# Box plots

p0 <- pairwise_database %>% 
  ggplot(aes(x = factor(sev_diff), y = prs_diff)) +
  geom_boxplot(aes(color = factor(sev_diff))) +    
  scale_x_discrete(labels=c("Concordant \npairs\nN=784",
                            "Discordant \npairs\nN=1615")) +
  stat_summary(fun="mean") +
  xlab("") +                          
  ylab("PRS difference") +
  ggtitle("A \n") +
  theme_classic() +
  theme(legend.position = "none") + 
  scale_color_manual(values=c("yellow4","blue4")) +
    theme(
    plot.margin = unit(c(0,0.5,0,0.5), "cm"),
    axis.text.x = element_text(colour="black", size = 9, margin = margin(10, 0, 0), family = "URWHelvetica"),
    axis.text.y = element_text(colour="black", size = 6, family = "URWHelvetica"),
    panel.grid.minor.y = element_blank(),
    axis.title.y = element_text(color = "black", size = 9, margin = margin(0,15,0), family = "URWHelvetica"),
    axis.title.x = element_text(color = "black", size = 9, margin = margin(0,15,0), family = "URWHelvetica"),
    plot.title = element_text(color = "black", size = 9, margin=margin(0,0,20), family = "URWHelvetica"), 
    plot.title.position = "plot"
  ) +
  geom_text(data=mean.degrees,
    aes(y=m+0.45, label=round(m,2)),
    color='black', position=position_dodge(0.8),
    size = 1.5, family = "URWHelvetica"
  )

p0

##### discordant group stratified

sample_size = pairwise_database %>% 
  group_by(diff_group) %>% 
  summarize(num=n())

sample_size

mean.degrees <- pairwise_database %>% 
  group_by(diff_group) %>% 
  summarize(m=mean(prs_diff))

mean.degrees

# Box plots - degree of mismatch phenotype

p1 <- pairwise_database %>% 
  ggplot(aes(x = factor(diff_group), y = prs_diff)) +
  geom_boxplot(aes(color = factor(diff_group))) +    
  scale_x_discrete(labels=c("0-grades\nN=784",
    "1-grade\nN=659",
                            "2-grades\nN=484",
                            "3-grades\nN=436",
                            "4-grades\nN=36")) +
  stat_summary(fun="mean") +
  xlab("Degree of phenotype heterogeneity") +                          
  ylab("PRS difference") +
  ggtitle("B \n") +
  theme_classic() +
  theme(legend.position = "none") + 
  scale_color_manual(values=c("yellow4","blue4", "blue4", "blue4", "blue4")) +
    theme(
    plot.margin = unit(c(0,0.5,0,0.5), "cm"),
    axis.text.x = element_text(colour="black", size = 6, margin = margin(10, 0, 0), family = "URWHelvetica"),
    axis.text.y = element_text(colour="black", size = 6, family = "URWHelvetica"),
    panel.grid.minor.y = element_blank(),
    axis.title.y = element_text(color = "black", size = 9, margin = margin(0,15,0), family = "URWHelvetica"),
    axis.title.x = element_text(color = "black", size = 9, margin = margin(0,15,0), family = "URWHelvetica"),
    plot.title = element_text(color = "black", size = 9, margin=margin(0,0,20), family = "URWHelvetica"), 
    plot.title.position = "plot"
  ) +
  geom_text(data=mean.degrees,
    aes(y=m+0.45, label=round(m,2)),
    color='black', position=position_dodge(0.8),
    size = 1.5, family = "URWHelvetica"
  )

p1

tiff("Figure1.tiff", units="cm", width=17, height=7, res=300)
p0 + p1
dev.off()