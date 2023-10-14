# R code to produce similarity matrix

library(readxl)
library(seriation)
library(tidyverse)

setwd(" ")

table <- read_excel("table_for_correlation_transposed.xlsx")

names <- table[ ,1]
names <- as.matrix(names)

table.mat <- as.matrix(table[ ,-1])

row.names(table.mat) <- names

# Compute dissimilarity matrix
dist_result <- dist(table.mat)
# Seriate objects, reorder rows based on their similarity

object_order <- seriate(dist_result)
# Extract object orders
head(get_order(object_order), 15)

pimage(dist_result, main = "Random order")
pimage(dist_result, order = object_order, main = "Reordered")

# Heamap of the raw data
pimage(scale(table.mat), main = "Random")

# Heatmap of the reordered data
pimage(scale(table.mat), order = c(object_order, NA), main = "Reordered")

# Standardize the data
df_scaled <- scale(table.mat, center = FALSE)

# Produce a heat map with optimally reordered dendrograms
# Hierarchical clustering is used to produce dendrograms
hmap(df_scaled, margin = c(7, 4), cexCol = 1, labRow = FALSE)

## Seriation ##

# Visualize the original data
bertinplot(
    table.mat, 
    options = list(panel = panel.circles)
)

# Seriate rows and columns using the bond energy algorithm (BEA)

orders <- seriate(table.mat, method = "BEA", control = list(rep = 10))
bertinplot(
    table.mat, order = orders, 
    options = list(panel = panel.circles)
)
