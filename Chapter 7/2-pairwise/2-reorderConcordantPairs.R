# R code to create pairwise comparisons within families

library(tidyverse)

concordant_database <- readRDS("pairwise_database.rds") %>% 
  filter(diff_group == 0) 

discordant_database <- readRDS("pairwise_database.rds") %>% 
  filter(diff_group != 0) 


### Switch signs
# Randomly re-assign +/- to model different relative minus the other for each pair

meanObt <- mean(concordant_database$prs_diff)
difference <- concordant_database$prs_diff

nreps <- 9066
set.seed(1086)

resampMeanDiff <- numeric(nreps)
for (i in 1:nreps) {
       signs <- sample( c(1,-1),length(difference), replace = TRUE)
       resamp <- difference * signs
       resampMeanDiff[i] <- mean(resamp)
       }

# Plot resulting distribution

resampMeanDiff.df <- as.data.frame(resampMeanDiff)

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

# How representative was the mean we by chance were using?
#original mean
concordant_database %>% 
  summarise(mean = mean(prs_diff))

#resampled distribution
summary(resampMeanDiff)

### We now know the optimal PRS difference mean value for each cohort/sub-cohort.
### We need to use this information to identify the optimal pair ordering to 
### achieve a mean PRS difference as close to the simulated mean.

opt_mean <- -0.00058 #update as required from above analysis

# Re-run the simulation and compare output values with optimal mean.

difference <- concordant_database$prs_diff

nreps <- 10000
set.seed(1086)

resampMeanDiff <- numeric(nreps)
for (i in 1:nreps) {
       signs <- sample( c(1,-1),length(difference), replace = TRUE)
       resamp <- difference * signs
       resampMeanDiff[i] <- mean(resamp)
       }

resampMeanDiff.df <- as.data.frame(resampMeanDiff)

opt_run <- which.min(abs(opt_mean - resampMeanDiff.df$resampMeanDiff))

# re-run simulations stopping at the optimal run

nreps2 <- opt_run
set.seed(1086)

resampMeanDiff2 <- numeric(nreps2)
for (i in 1:nreps2) {
       signs <- sample( c(1,-1),length(difference), replace = TRUE)
       resamp2 <- difference * signs
       resampMeanDiff2[i] <- mean(resamp2)
       }

new_prs_diff <- as.data.frame(resamp2)

# check it worked
mean(new_prs_diff$resamp2)

# new_prs_diff$resamp2 becomes the new prs_diff for the concordant group (total cohort)

concordant_database2 <- cbind(concordant_database, new_prs_diff) 

concordant_database3 <- concordant_database2 %>% 
  mutate(prs_diff = resamp2) %>% 
  select(-(resamp2))

glimpse(concordant_database3)

# rejoin with discordant pair data
total_database <- rbind(disconcordant_database, concordant_database3)

total_database %>% 
  group_by(diff_group, known_mutn.x) %>% 
  summarise(mean = mean(prs_diff))

saveRDS(total_database, "pairwise_database_optimised.rds")

