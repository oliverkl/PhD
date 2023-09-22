library(tidyverse)

#read in pairwise data
#filter for conconcordant pairs

concordant_database <- readRDS("df.rds") %>% 
  filter(diff_group == 0) 

#Randomly re-assign +/- to model different relative minus the other for each pair

meanObt <- mean(concordant_database$prs_diff)
difference <- concordant_database$prs_diff

nreps <- 10000
set.seed(1086)

resampMeanDiff <- numeric(nreps)
for (i in 1:nreps) {
       signs <- sample( c(1,-1),length(difference), replace = T)
       resamp <- difference * signs
       resampMeanDiff[i] <- mean(resamp)
       }

resampMeanDiff.df <- as.data.frame(resampMeanDiff)

#summarise results
summary(resampMeanDiff.df$resampMeanDiff)
meanObt <- mean(resampMeanDiff.df$resampMeanDiff)

#plot results
p1 <- ggplot(resampMeanDiff.df, aes(x=resampMeanDiff)) +
   geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", bins = 30)+
  geom_density(alpha=.2, fill="brown1") +
    labs(title="Concordant pairs \nTotal cohort",
        x ="Mean PRS difference \n (+/-) signs switched x10000", y = "Density") +
  geom_vline(xintercept=c(meanObt), color="brown3", #to see where our observed mean prs diff sits
             linetype="dashed", show.legend = FALSE) +
    theme_classic() +
    theme(
    plot.margin = unit(c(0,0.5,0,0.5), "cm"),
    axis.text.x = element_text(colour="black", size = 6, margin = margin(10, 0, 0), family = "URWHelvetica"),
    axis.text.y = element_text(colour="black", size = 6, family = "URWHelvetica"),
    panel.grid.minor.y = element_blank(),
    axis.title.y = element_text(color = "black", size = 9, margin = margin(0,15,0), family = "URWHelvetica"),
    axis.title.x = element_text(color = "black", size = 7, margin = margin(0,15,0), family = "URWHelvetica"),
    plot.title = element_text(color = "black", size = 9, margin=margin(0,0,20), family = "URWHelvetica"), 
    plot.title.position = "plot"
  )

p1
