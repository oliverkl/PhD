# R code to compare PRS distributions between cases and controls

library(tidyverse)
library(lme4)
library(lmerTest)
library(GMMAT)
library(ggpubr)
library(ggridges)
library(patchwork)

## Read in required datasets

# clinical
# PCAs
# PRSs

prs.dir <- " "
clin.dir <- " "
pca.dir <- " "

prs.data <- readRDS(file.path(prs.dir, "prs.data.rds"))
clin.data <- readRDS(file.path(clin.dir, "clin.data.rds"))
pca.data <- readRDS(file.path(pca.dir, "pca.data.rds"))

prs_database <- clin.data %>% 
  left_join(prs.data, by = c("IID" = "sample.id")) %>% 
  left_join(pca.data, by = "sample.id")

## Logistic regression ##

# status - case v control (1, 0)

### lme4 ###
model_01 <- lme4::glmer(status ~ PRS + SEX + PC1 + PC2 + PC3 + PC4 + PC5 +
                     (1|FID), family=binomial, data = prs_database)

summary(model_01)

# calc ORs
cc <- confint(model_01, parm = "beta_", method = "Wald")
ctab <- cbind(est = fixef(model_foc01), cc)
rtab <- exp(ctab) #exponentiate to get ORs
print(rtab, digits = 3)

### gmmat ###

model_02 <- glmmkin(status ~ PRS + SEX + PC1 + PC2 + PC3 + PC4 + PC5,
                  data = prs_database, kins = relmat, id = "IID", 
                  family = binomial(link = "logit"))

# calc ORs
coeff <- coefficients(model_02)
se <- diag(sqrt(abs(model_02$cov)))
z <- coeff / se
p <- 2*pnorm(-abs(z))
or <- exp(coeff)
lower <- or * exp(-1.96 * se)
upper <- or * exp(1.96 * se)

model_result <- data.frame(coeff, se, z, or, lower, upper, p)
model_result

## Plots ##

# Ridges plot of score distributions

sd <- sd(prs_database$PRS)
mean <- mean(prs_database$PRS)

Epilepsy <- factor(prs_database$status, levels = c(0 , 1), labels = c("No", "Yes"))

p1 <- ggplot(prs_database, aes(x = PRS, y = Epilepsy, color = Epilepsy, fill = Epilepsy)) + 
  stat_density_ridges(alpha = 0.2, show.legend = TRUE, bandwidth = 0.5,
                      quantile_lines = TRUE, quantile_fun=function(x,...)mean(x)) + 
  labs(title ="A", x = "Polygenic Risk Score\n(Units of SD)", y = "Density") +
  theme_ridges(center_axis_labels = TRUE) + # simpler theme
  scale_fill_manual(values = c("gray30", "coral3")) +
  scale_color_manual(values = c("gray30", "coral3")) +
  scale_y_discrete(expand = expansion(add = c(0, 1.9))) +
  coord_cartesian(xlim = c(-4, 4)) +
  theme(
    plot.margin = unit(c(0,2,0,0), "cm"),
    legend.position = c(0.95, 0.95),
    legend.title = element_text(size = 9, family = "URWHelvetica"),
    legend.text = element_text(size = 6, family = "URWHelvetica"),
    axis.text.x = element_text(colour="black", size = 6, family = "URWHelvetica"),
    axis.text.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.title.y = element_text(color = "black", size = 9, margin = margin(0,0,10), family = "URWHelvetica"),
    axis.title.x = element_text(color = "black", size = 6, margin = margin(-20), family = "URWHelvetica"),
    plot.title = element_text(color = "black", size = 9, margin=margin(0,0,20), face="plain", family = "URWHelvetica"), 
    plot.title.position = "plot"
  ) +
  guides(fill = guide_legend(reverse = TRUE)) + guides(color = guide_legend(reverse = TRUE))

p1

p2 <- prs_database %>%  
  ggerrorplot(x = "status", y = "PRS", desc_stat = "mean_ci", size = 0.8, color = "status",
              order = c("0", "1")) +
  ggtitle("B") +
  labs(y = "PRS mean +/- 95% CI", x = " ") +
  theme_minimal() +
  scale_fill_manual(values = c("gray30", "coral3")) +
  scale_color_manual(values = c("gray30", "coral3")) +
  theme(legend.position = "right", text = element_text(size=16)) +
  scale_x_discrete(labels = c("Unaffected", "Affected")) +
  font("ylab", size = 6, family = "URWHelvetica") +
  font("xlab", size = 6, family = "URWHelvetica") +
  theme(
    plot.margin = unit(c(0,0.5,0,0), "cm"),
    legend.position = "none",
    axis.text.x = element_text(colour="black", size = 6, margin = margin(10, 0, 0), family = "URWHelvetica"),
    axis.text.y = element_text(colour="black", size = 6, family = "URWHelvetica"),
    panel.grid.minor.y = element_blank(),
    axis.title.y = element_text(color = "black", size = 9, margin = margin(0,15,0), family = "URWHelvetica"),
    plot.title = element_text(color = "black", size = 9, margin=margin(0,0,20), family = "URWHelvetica"), 
    plot.title.position = "plot"
  ) +
  geom_signif(comparisons = list(c(2,3)), annotations = "list(italic(p)~'='~0.01)", parse = TRUE, y_position = 0.2, tip_length = 0.004, color = "black", family = "URWHelvetica", textsize = 2.5) 

p1 + ggpar(p2, ylim = c(-0.15, 0.65), ticks = FALSE)

tiff("Fig1_PRS.tiff", units="cm", width=17, height=7, res=300)
p1 + ggpar(p2, ylim = c(-0.15, 0.65), ticks = FALSE)
dev.off()

