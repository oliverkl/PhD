gcta.dir <- " "

eigenvals <- read_tsv(file.path(gcta.dir, "filename.eigenval"),
                      col_names = "eigenval") 

eigenvals1 <- eigenvals[1:10, ]

total <- sum(eigenvals1$eigenval)

eigenvals2 <- eigenvals1 %>% 
  mutate("% variance explained" = (eigenval/total)*100) %>% 
  mutate(component = paste0("PC", row.names(eigenvals1)))

eigenvals2$component <- factor(eigenvals2$component, levels = eigenvals2$component)

eigenvals2 %>%
  ggplot(aes(x=component, y=`% variance explained`, group = 1)) +
  geom_line(linetype = "dashed") +
  geom_point() + theme_classic()