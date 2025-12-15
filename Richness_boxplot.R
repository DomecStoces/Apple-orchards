library(tidyverse)
library(readxl)
library(ggplot2)
library(ggpubr)

Sady_final_sorted <- read_excel("Sady_final_sorted.xlsx", sheet = 1)
long_df <- Sady_final_sorted %>%
  pivot_longer(
    cols = -(Locality:Shrubs),
    names_to = "Species",
    values_to = "Abundance"
  )

rich_per_trap_month <- long_df %>%
  group_by(Trap, Month, Treatment) %>%
  summarise(
    richness = sum(Abundance > 0, na.rm = TRUE),
    .groups = "drop"
  )

rich_summary <- rich_per_sample %>%
  group_by(Treatment) %>%
  summarise(
    mean_rich = mean(richness),
    se_rich   = sd(richness) / sqrt(n()),
    n         = n(),
    .groups = "drop"
  )

p<-ggplot(rich_per_trap_month,
       aes(x = Treatment, y = richness)) +
  geom_boxplot(fill = "white",
               color = "black",
               width = 0.6,
               outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.7,
    color="black"
  ) +
  
  labs(
    x = "Režim hospodaření",
    y = "Počet druhů na past/ za měsíc"
  ) +
  theme_minimal() +
  theme(
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.ticks.length = unit(2, "mm"),
    axis.line = element_line(color = "black"),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 15)
  )
ps<-p +
  stat_compare_means(
    comparisons = list(c("BIO","EXT"), c("BIO","INT"), c("EXT","INT")),
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "p.format"
  )
ps

tiff("Druhova_bohatost.tiff", units = "in", width = 8, height = 10, res = 600)
ps
dev.off()
pdf("Druhova_bohatost.pdf", width = 8, height = 10)
ps
dev.off()
