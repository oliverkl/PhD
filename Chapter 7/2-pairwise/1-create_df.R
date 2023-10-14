# R code to create pairwise comparisons within families

library(tidyverse)

gefsplus_database <- readRDS("gefsplus_database.rds") 

#select variables
family.pairs <- gefsplus_database %>% 
  filter(sev_scale %in% c("1A", "1B", "2", "3", "4A", "4B", "4C", "5")) %>% 
  select(sample.id, `allEpi3_0.5.norm`, `Focal3_0.5.norm`, `GGE3_0.5.norm`, FID, IID, SEX, BroadPHENO, FS, FS.bin, sev_scale, severity, known_mutn,
         Gene, Gene2, PFN, PFN2, cohort, outcome2)

#match on PFN2 instead of FID as was losing singleton SCN1B samples
family.pairs2 <- family.pairs %>% 
  left_join(family.pairs, by = "PFN2", relationship = "many-to-many") %>% 
  mutate(duplicate = case_when(sample.id.x == sample.id.y ~ "yes",
                               sample.id.x != sample.id.y ~ "no")) %>% 
  filter(duplicate == "no") %>% 
  unite(col = "z", c(IID.x, IID.y), sep = "_", remove = FALSE) %>% 
  unite(col = "y", c(IID.y, IID.x), sep = "_", remove = FALSE) %>% 
  mutate(sev_diff = case_when(severity.x == severity.y ~ "no",
                               severity.x != severity.y ~ "yes")) %>% 
  mutate(low_to_high = case_when(severity.x < severity.y ~ "yes",
                                 severity.x > severity.y ~ "no",
                                 TRUE ~ "equal")) %>% 
  mutate(mild_all = case_when(low_to_high == "yes" ~ allEpi3_0.5.norm.x,
                          low_to_high == "no" ~ allEpi3_0.5.norm.y,
                          low_to_high == "equal" ~ allEpi3_0.5.norm.y)) %>% 
  mutate(severe_all = case_when(low_to_high == "no" ~ allEpi3_0.5.norm.x,
                          low_to_high == "yes" ~ allEpi3_0.5.norm.y,
                          low_to_high == "equal" ~ allEpi3_0.5.norm.x)) %>% 
  mutate(prs_diff_all = severe_all - mild_all) %>% 
  mutate(diff_group = case_when(low_to_high == "yes" ~ (severity.y - severity.x),
                                low_to_high == "no" ~ (severity.x - severity.y),
                                TRUE ~ 0)) %>% 
  mutate(mild_gge = case_when(low_to_high == "yes" ~ GGE3_0.5.norm.x,
                          low_to_high == "no" ~ GGE3_0.5.norm.y,
                          low_to_high == "equal" ~ GGE3_0.5.norm.y)) %>% 
  mutate(severe_gge = case_when(low_to_high == "no" ~ GGE3_0.5.norm.x,
                          low_to_high == "yes" ~ GGE3_0.5.norm.y,
                          low_to_high == "equal" ~ GGE3_0.5.norm.x)) %>% 
  mutate(prs_diff_gge = severe_gge - mild_gge) %>% 
  mutate(mild_focal = case_when(low_to_high == "yes" ~ Focal3_0.5.norm.x,
                          low_to_high == "no" ~ Focal3_0.5.norm.y,
                          low_to_high == "equal" ~ Focal3_0.5.norm.y)) %>% 
  mutate(severe_focal = case_when(low_to_high == "no" ~ Focal3_0.5.norm.x,
                          low_to_high == "yes" ~ Focal3_0.5.norm.y,
                          low_to_high == "equal" ~ Focal3_0.5.norm.x)) %>% 
  mutate(prs_diff_focal = severe_focal - mild_focal) %>% 
  mutate(diff_abs = round(abs(prs_diff_all), 5)) %>% 
  group_by(PFN2, diff_group) %>% distinct(diff_abs, .keep_all= TRUE)

#sense check worked as expected
addmargins(table(family.pairs2$sev_diff, family.pairs2$low_to_high))

family.pairs3 <- family.pairs2 %>% 
  unite(pheno_comp1, c(severity.x,severity.y), sep = "-", remove = FALSE) %>% 
  unite(pheno_comp2, c(severity.y,severity.x), sep = "-", remove = FALSE) %>% 
  mutate(pheno_comp = case_when(low_to_high == "equal" ~ pheno_comp1,
                                low_to_high == "yes" ~ pheno_comp1,
                                low_to_high == "no" ~ pheno_comp2)) %>% 
  mutate(pheno_pairs = case_when(pheno_comp == "1-1" ~ "unaffected pair",
                                 pheno_comp == "2-2" ~ "FS pair",
                                 pheno_comp == "3-3" ~ "FS+ pair",
                                 pheno_comp == "4-4" ~ "Epilepsy pair",
                                 pheno_comp == "5-5" ~ "DEE pair",
                                 pheno_comp == "1-2" ~ "FS-unaffected",
                                 pheno_comp == "1-3" ~ "FS+-unaffected",
                                 pheno_comp == "1-4" ~ "Epilepsy-unaffected",
                                 pheno_comp == "1-5" ~ "DEE-unaffected",
                                 pheno_comp == "2-3" ~ "FS+-FS",
                                 pheno_comp == "2-4" ~ "Epilepsy-FS",
                                 pheno_comp == "2-5" ~ "DEE-FS",
                                 pheno_comp == "3-4" ~ "Epilepsy-FS+",
                                 pheno_comp == "3-5" ~ "DEE-FS+",
                                 pheno_comp == "4-5" ~ "DEE-Epilepsy",
                                 TRUE ~ pheno_comp))

# sense check numbers
addmargins(table(family.pairs3$pheno_comp, family.pairs3$sev_diff))
addmargins(table(family.pairs3$pheno_comp, family.pairs3$diff_group))
addmargins(table(family.pairs3$pheno_comp, family.pairs3$pheno_pairs))

addmargins(table(family.pairs3$diff_group, family.pairs3$sev_diff))

addmargins(table(family.pairs3$sev_diff, family.pairs3$Gene2.x))

# total cohort
family.pairs4 <- family.pairs3 %>% 
  mutate(pheno_diff = case_when(sev_diff == "no" ~ 0,
                                sev_diff == "yes" ~ 1,
                                TRUE ~ 0))

### Add in Kinship values so can include in models

kinship.values <- read_table("/stornext/Bioinf/data/lab_bahlo/users/oliver/projects/PhD/PRS/severity/7-ibd/ibd/king.kin") %>% 
  mutate(pair_id = paste0(ID1, "_", ID2)) %>% 
  select(pair_id, Kinship)

family.pairs5 <- family.pairs4 %>% 
  left_join(kinship.values, by = c("z" = "pair_id")) %>% 
  left_join(kinship.values, by = c("y" = "pair_id")) %>% 
  mutate(Kinship = coalesce(Kinship.x, Kinship.y))

addmargins(table(family.pairs5$pheno_diff, useNA = "ifany"))

family.pairs6 <- family.pairs5 %>% 
  select(-c(Kinship.y, Kinship.x, duplicate, diff_abs)) %>% 
  ungroup()

saveRDS(family.pairs6, "pairwise_database.rds") 