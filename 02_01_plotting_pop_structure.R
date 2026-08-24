

library(tidyverse)

setwd("D:/Landscape_210")
admixture <- read.csv("Admixture_532_results.csv", header = TRUE)

# ── 1. Reshape ────────────────────────────────────────────────────────────────
admix_long <- admixture %>%
  select(Geno_vcf, Geographic_area, order_barplot,
         A1, A2,
         B1, B2, B3,
         c1, c2, c3, c4,
         d1, d2, d3, d4, d5,
         e1, e2, e3, e4, e5, e6) %>%
  pivot_longer(
    cols      = -c(Geno_vcf, Geographic_area, order_barplot),
    names_to  = "cluster",
    values_to = "Q"
  ) %>%
  mutate(
    K = case_when(
      cluster %in% c("A1", "A2")                          ~ "K=2",
      cluster %in% c("B1", "B2", "B3")                    ~ "K=3",
      cluster %in% c("c1", "c2", "c3", "c4")              ~ "K=4",
      cluster %in% c("d1", "d2", "d3", "d4", "d5")        ~ "K=5",
      cluster %in% c("e1", "e2", "e3", "e4", "e5", "e6")  ~ "K=6"
    ),
    K = factor(K, levels = c("K=2", "K=3", "K=4", "K=5", "K=6"))
  )

# ── 2. Filter + harmonise area names ─────────────────────────────────────────
geo_levels <- c(
  "Turkey", "Balkanic Europe", "North Africa", "Iberian Peninsula",
  "East Europe", "Italy", "Middle East", "Ethiopia",
  "Asia", "Indian Continent", "Central America", "South America", "USA"
)

admix_long <- admix_long %>%
  filter(
    !is.na(Geographic_area),
    Geographic_area != "",
    tolower(trimws(Geographic_area)) != "others"
  ) %>%
  mutate(
    Geographic_area = trimws(Geographic_area),
    Geographic_area = recode(Geographic_area,
                             "Iberian Peninsola"  = "Iberian Peninsula",
                             "Indian continent"   = "Indian Continent",
                             "\nIndian continent" = "Indian Continent",
                             "\tIndian continent" = "Indian Continent"
    ),
    Geographic_area = factor(Geographic_area, levels = geo_levels)
  ) %>%
  filter(!is.na(Geographic_area)) %>%
  # ── Sort: geographic area order first, then order_barplot within group ──
  arrange(Geographic_area, order_barplot)

# ── 3. Build ind_order WITH gaps (FIXED & CLEAN) ─────────────────────────────

gap_size <- 3

# Unique individuals sorted properly
ind_meta <- admix_long %>%
  filter(K == "K=2") %>%
  distinct(Geno_vcf, Geographic_area, order_barplot) %>%
  arrange(Geographic_area, order_barplot) %>%
  group_by(Geographic_area) %>%
  mutate(within_rank = row_number()) %>%
  ungroup()

# Count individuals per group
group_sizes <- ind_meta %>%
  count(Geographic_area, name = "n") %>%
  arrange(Geographic_area)

# Compute starting x for each group (with gaps)
group_sizes <- group_sizes %>%
  mutate(
    start = cumsum(lag(n + gap_size, default = 0)) + 1
  )

# Join start positions back
ind_meta <- ind_meta %>%
  left_join(group_sizes, by = "Geographic_area") %>%
  mutate(
    x_pos = start + within_rank - 1
  )

# Join x_pos back to the long data
admix_long <- admix_long %>%
  left_join(ind_meta %>% select(Geno_vcf, x_pos), by = "Geno_vcf")

# ── 4. Boundary positions for separators and x-axis labels ───────────────────
group_bounds <- ind_meta %>%
  group_by(Geographic_area) %>%
  summarise(
    xmin = min(x_pos),
    xmax = max(x_pos),
    .groups = "drop"
  ) %>%
  mutate(
    xmid      = (xmin + xmax) / 2,
    sep_after = xmax + gap_size / 2   # separator in the middle of the gap
  )

sep_x <- group_bounds$sep_after[-nrow(group_bounds)]

# ── 5. Colour palette ─────────────────────────────────────────────────────────
cluster_colors <- c(
  A1 = "blue",      A2 = "darkorange",
  B1 = "blue",      B2 = "darkorange",  B3 = "#228833",
  c1 = "blue",      c2 = "darkorange",  c3 = "#228833",  c4 = "gold1",
  d1 = "blue",      d2 = "darkorange",  d3 = "#228833",  d4 = "gold1",d5 = "darkgoldenrod4",
  e1 = "blue",      e2 = "darkorange",  e3 = "#228833",  e4 = "gold1",e5 = "darkgoldenrod4",     e6 = "darkorchid4"
)

# ── 6. Plot ───────────────────────────────────────────────────────────────────
p <- ggplot(admix_long, aes(x = x_pos, y = Q, fill = cluster)) +
  
  geom_bar(stat = "identity", width = 1, linewidth = 0) +
  
  # Thin white separator lines in the middle of each gap
  geom_vline(xintercept = sep_x, colour = "grey40", linewidth = 0.4,
             linetype = "solid") +
  
  facet_grid(K ~ ., switch = "y") +
  
  scale_fill_manual(values = cluster_colors) +
  
  scale_x_continuous(
    expand = expansion(mult = c(0.005, 0.005)),
    breaks = group_bounds$xmid,
    labels = group_bounds$Geographic_area
  ) +
  
  scale_y_continuous(
    breaks = c(0, 0.5, 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  labs(x = NULL, y = "Ancestry proportion (Q)") +
  
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle  = 40,
      hjust  = 1,
      vjust  = 1,
      size   = 12,
      face   = "italic",
      colour = "grey20"
    ),
    axis.ticks.x        = element_blank(),
    axis.text.y         = element_text(size = 7.5, colour = "grey40"),
    axis.title.y        = element_text(size = 10, margin = margin(r = 6)),
    strip.text.y.left   = element_text(face = "bold", size = 12, angle = 0,
                                       margin = margin(l = 3, r = 6)),
    strip.background    = element_rect(fill = "grey88", color = NA),
    strip.placement     = "outside",
    panel.grid          = element_blank(),
    panel.border        = element_rect(fill = NA, colour = "grey70",
                                       linewidth = 0.3),
    panel.spacing.y     = unit(0.4, "lines"),
    #plot.margin         = margin(t = 8, r = 15, b = 10, l = 10),
    legend.position     = "none"
  )
p
# ── 7. Export ─────────────────────────────────────────────────────────────────

ggsave("admixture_plot.tiff", plot = p, width = 20, height = 9,
       units = "in", dpi = 600)

message("Done — admixture_plot.png / .pdf")










