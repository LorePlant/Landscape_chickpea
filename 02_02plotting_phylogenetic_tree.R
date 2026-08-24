# ── 0. Install (run once) ────────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("ggtree", "treeio"))
install.packages(c("ape", "phytools", "ggplot2", "dplyr", "RColorBrewer", "ggnewscale"))

setwd("D:/Landscape_210")
library(ape)
library(ggplot2)
library(ggtree)
library(dplyr)

# ── 1. Read tree ──────────────────────────────────────────────────────────────
tree <- read.nexus("Country_mdist_K3_DESI.nj.nex")

# ── 2. Read color assignments from file ───────────────────────────────────────
color_id <- read.table("C:/Users/rocchett/Downloads/colorid_K3_DESI.txt",
                       header = FALSE, col.names = "color")

# ── 3. Build metadata from file ───────────────────────────────────────────────
meta <- data.frame(
  label   = as.character(1:101),
  color   = color_id$color,
  K_group = case_when(
    color_id$color == "mediumblue" ~ "K1",
    color_id$color == "orange"     ~ "K2",
    color_id$color == "green3"     ~ "K3",
    TRUE ~ NA_character_
  ),
  stringsAsFactors = FALSE
)

# Quick check
table(meta$K_group)
head(meta)

# ── 4. Get tree layout coordinates ────────────────────────────────────────────
tree_data <- fortify(tree, layout = "daylight", branch.length = "branch.length")

# ── 5. Recursive function to get all descendants of a node ───────────────────
get_descendants <- function(tree, node) {
  children <- vector("list", max(tree$edge))
  for (i in seq_len(nrow(tree$edge))) {
    p <- tree$edge[i, 1]; ch <- tree$edge[i, 2]
    children[[p]] <- c(children[[p]], ch)
  }
  get_desc_recursive <- function(nd) {
    ch <- children[[nd]]
    if (is.null(ch)) return(nd)
    c(ch, unlist(lapply(ch, get_desc_recursive)))
  }
  get_desc_recursive(node)
}

# ── 6. Propagate K group from tips to internal nodes ─────────────────────────
n_tips    <- length(tree$tip.label)
n_nodes   <- tree$Nnode
all_nodes <- n_tips + n_nodes

tip_lookup <- setNames(meta$K_group, meta$label)

node_group <- character(all_nodes)
for (i in seq_len(n_tips)) {
  node_group[i] <- tip_lookup[tree$tip.label[i]]
}
for (nd in (n_tips + 1):all_nodes) {
  desc     <- get_descendants(tree, nd)
  tip_desc <- desc[desc <= n_tips]
  groups   <- unique(node_group[tip_desc])
  node_group[nd] <- if (length(groups) == 1) groups else "mixed"
}

# ── 7. Build edge data ────────────────────────────────────────────────────────
edge_groups <- node_group[tree$edge[, 2]]

node_coords <- tree_data %>% select(node, x, y)

edge_data <- data.frame(
  child  = tree$edge[, 2],
  parent = tree$edge[, 1],
  group  = edge_groups
) %>%
  left_join(node_coords, by = c("child"  = "node")) %>%
  left_join(node_coords, by = c("parent" = "node"), suffix = c("", "_parent"))

tip_data <- tree_data %>%
  filter(isTip) %>%
  left_join(meta, by = "label")


rotation <- 90 - (-95.1949)  # = 185.19 degrees
# ── Rotation function ─────────────────────────────────────────────────────────
rotate_coords <- function(x, y, angle_deg) {
  angle_rad <- angle_deg * pi / 180
  list(
    x = x * cos(angle_rad) - y * sin(angle_rad),
    y = x * sin(angle_rad) + y * cos(angle_rad)
  )
}

# ── Apply rotation to edge_data ───────────────────────────────────────────────
rotated_child  <- rotate_coords(edge_data$x,        edge_data$y,        rotation)
rotated_parent <- rotate_coords(edge_data$x_parent, edge_data$y_parent, rotation)
edge_data_r <- edge_data %>%
  mutate(x        = rotated_child$x,
         y        = rotated_child$y,
         x_parent = rotated_parent$x,
         y_parent = rotated_parent$y)

# ── Apply rotation to tip_data ────────────────────────────────────────────────
rotated_tips <- rotate_coords(tip_data$x, tip_data$y, rotation)
tip_data_r <- tip_data %>%
  mutate(x = rotated_tips$x,
         y = rotated_tips$y)

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot() +
  geom_segment(data = edge_data_r,
               aes(x = x_parent, y = y_parent, xend = x, yend = y,
                   color = group),
               linewidth = 1, lineend = "round") +
  geom_point(data = tip_data_r,
             aes(x = x, y = y),
             shape = 21, size = 3.5, fill = "darkgoldenrod4")+
             #stroke = 0.4, color = "black") +
  scale_color_manual(
    values = k_colors,
    breaks = c("K1", "K2", "K3"),
    labels = c("K1 (Blue)", "K2 (Orange)", "K3 (Green)"),
    name   = "K group"
  ) +
  coord_equal() +
  theme_void() +
   guides(color = guide_legend(override.aes = list(linewidth = 3)))

print(p)
ggsave("tree_Kgroups_DESI.png", p, width = 10, height = 10, dpi = 300)
ggsave("tree_Kgroups_final.pdf", p, width = 10, height = 10, bg = "black")



