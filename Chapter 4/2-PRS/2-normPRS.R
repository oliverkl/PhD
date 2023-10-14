# R code to curate PRSs

## Scores for GGE, Focal and All-Epilepsy (2018 and 2022) were
## calculated using PRSice.  Scores were calculated based on multiple p-value
## thresholds: 1e-7, 1e-5, 0.01, 0.05, 0.1, 0.5. Also FS sumstats


# Read in PRS results and rename p-value thresholds:

library(tidyverse)

prs.dir <- " "

allEpi_prs_2018 <- read_table(file.path(prs.dir, "allEpi_2018/allEpi_2018-PRS_all_gsa.all.score"),
                                  col_names = TRUE) %>% 
  unite("sample.id", FID:IID, sep = "-", remove = FALSE) %>% 
  rename("allEpi2_0.1" = "0.1",
         "allEpi2_0.5" = "0.5")

allEpi_prs_2022 <- read_table(file.path(prs.dir, "allEpi_2022/allEpi_2022-PRS_merged.all.score"),
                                  col_names = TRUE) %>% 
  unite("sample.id", FID:IID, sep = "-", remove = FALSE) %>% 
  rename("allEpi3_0.1" = "0.1",
         "allEpi3_0.5" = "0.5")

focal_prs_2018 <- read_table(file.path(prs.dir, "Focal_2018/Focal_2018-PRS_all_gsa.all.score"),
                                  col_names = TRUE) %>% 
  unite("sample.id", FID:IID, sep = "-", remove = TRUE) %>%
  rename("Focal2_0.1" = "0.1",
         "Focal2_0.5" = "0.5")

focal_prs_2022 <- read_table(file.path(prs.dir, "Focal_2022/Focal_2022-PRS_merged.all.score"),
                                  col_names = TRUE) %>% 
  unite("sample.id", FID:IID, sep = "-", remove = TRUE) %>%
  rename("Focal3_0.1" = "0.1",
         "Focal3_0.5" = "0.5")

gge_prs_2018 <- read_table(file.path(prs.dir, "GGE_2018/GGE_2018-PRS_all_gsa.all.score"),
                                  col_names = TRUE) %>% 
  unite("sample.id", FID:IID, sep = "-", remove = TRUE) %>% 
  rename("GGE2_0.1" = "0.1",
         "GGE2_0.5" = "0.5")

gge_prs_2022 <- read_table(file.path(prs.dir, "GGE_2022/GGE_2022-PRS_merged.all.score"),
                                  col_names = TRUE) %>% 
  unite("sample.id", FID:IID, sep = "-", remove = TRUE) %>% 
  rename("GGE3_0.1" = "0.1",
         "GGE3_0.5" = "0.5")

fs_prs_2022 <- read_table(file.path(prs.dir, "FS_2022/FS_2022-PRS_merged.all.score"),
                                  col_names = TRUE) %>% 
  unite("sample.id", FID:IID, sep = "-", remove = TRUE) %>% 
  rename("FS_0.1" = "0.1",
         "FS_0.5" = "0.5")

master_prs1 <- left_join(allEpi_prs_2018, focal_prs_2018, by = "sample.id") %>% 
  left_join(gge_prs_2018, by = "sample.id") %>% 
  left_join(allEpi_prs_2022, by = "sample.id") %>% 
  left_join(focal_prs_2022, by = "sample.id") %>% 
  left_join(gge_prs_2022, by = "sample.id") %>% 
  left_join(fs_prs_2022, by = "sample.id") 


# Replace scores with normalised values based on control scores.
## First read in clinical data so that can identify/subset control samples.

master_prs1 <- master_prs1 %>% 
  mutate(status = case_when(str_detect(sample.id, "QSK") ~ "Control",
                                   TRUE ~ "Case"))

table(master_prs1$status)

# Normalise scores:

#All Epilepsy PRS models
master_prs1$allEpi2_0.5.norm <- (master_prs1$allEpi2_0.5 - mean(master_prs1$allEpi2_0.5[master_prs1$status == "Control"])) / sd(master_prs1$allEpi2_0.5[master_prs1$status == "Control"])
master_prs1$allEpi2_0.1.norm <- (master_prs1$allEpi2_0.1 - mean(master_prs1$allEpi2_0.1[master_prs1$status == "Control"])) / sd(master_prs1$allEpi2_0.1[master_prs1$status == "Control"])

master_prs1$allEpi3_0.5.norm <- (master_prs1$allEpi3_0.5 - mean(master_prs1$allEpi3_0.5[master_prs1$status == "Control"])) / sd(master_prs1$allEpi3_0.5[master_prs1$status == "Control"])
master_prs1$allEpi3_0.1.norm <- (master_prs1$allEpi3_0.1 - mean(master_prs1$allEpi3_0.1[master_prs1$status == "Control"])) / sd(master_prs1$allEpi3_0.1[master_prs1$status == "Control"])

## Focal PRS models
master_prs1$Focal2_0.5.norm <- (master_prs1$Focal2_0.5 - mean(master_prs1$Focal2_0.5[master_prs1$status == "Control"])) / sd(master_prs1$Focal2_0.5[master_prs1$status == "Control"])
master_prs1$Focal2_0.1.norm <- (master_prs1$Focal2_0.1 - mean(master_prs1$Focal2_0.1[master_prs1$status == "Control"])) / sd(master_prs1$Focal2_0.1[master_prs1$status == "Control"])

master_prs1$Focal3_0.5.norm <- (master_prs1$Focal3_0.5 - mean(master_prs1$Focal3_0.5[master_prs1$status == "Control"])) / sd(master_prs1$Focal3_0.5[master_prs1$status == "Control"])
master_prs1$Focal3_0.1.norm <- (master_prs1$Focal3_0.1 - mean(master_prs1$Focal3_0.1[master_prs1$status == "Control"])) / sd(master_prs1$Focal3_0.1[master_prs1$status == "Control"])


## GGE PRS models
master_prs1$GGE2_0.5.norm <- (master_prs1$GGE2_0.5 - mean(master_prs1$GGE2_0.5[master_prs1$status == "Control"])) / sd(master_prs1$GGE2_0.5[master_prs1$status == "Control"])
master_prs1$GGE2_0.1.norm <- (master_prs1$GGE2_0.1 - mean(master_prs1$GGE2_0.1[master_prs1$status == "Control"])) / sd(master_prs1$GGE2_0.1[master_prs1$status == "Control"])

master_prs1$GGE3_0.5.norm <- (master_prs1$GGE3_0.5 - mean(master_prs1$GGE3_0.5[master_prs1$status == "Control"])) / sd(master_prs1$GGE3_0.5[master_prs1$status == "Control"])
master_prs1$GGE3_0.1.norm <- (master_prs1$GGE3_0.1 - mean(master_prs1$GGE3_0.1[master_prs1$status == "Control"])) / sd(master_prs1$GGE3_0.1[master_prs1$status == "Control"])

# FS model
master_prs1$FS_0.5.norm <- (master_prs1$FS_0.5 - mean(master_prs1$FS_0.5[master_prs1$status == "Control"])) / sd(master_prs1$FS_0.5[master_prs1$status == "Control"])
master_prs1$FS_0.1.norm <- (master_prs1$FS_0.1 - mean(master_prs1$FS_0.1[master_prs1$status == "Control"])) / sd(master_prs1$FS_0.1[master_prs1$status == "Control"])

## Sense check it worked

#all Epilepsy
master_prs1 %>% group_by(status) %>%
  summarise(`20%`=quantile(allEpi2_0.5.norm, probs=0.2),
            `50%`=quantile(allEpi2_0.5.norm, probs=0.5),
            `80%`=quantile(allEpi2_0.5.norm, probs=0.8),
            avg=mean(allEpi2_0.5.norm),
            n=n()) %>% knitr::kable("markdown")

## Remove the raw PRS columns and keep just the normalised scores.

master_prs2 <- master_prs1 %>% 
  select((c(sample.id, FID.x, IID.x, allEpi2_0.1.norm, allEpi2_0.5.norm,
             Focal2_0.1.norm, Focal2_0.5.norm, GGE2_0.1.norm, GGE2_0.5.norm, allEpi3_0.1.norm, allEpi3_0.5.norm,
             Focal3_0.1.norm, Focal3_0.5.norm, GGE3_0.1.norm, GGE3_0.5.norm, FS_0.1.norm, FS_0.5.norm)))

# Save this set of scores for downstream analysis.

saveRDS(master_prs2, file.path(prs.dir, "PRScoresNormalised_date.rds"))
