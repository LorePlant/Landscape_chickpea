############################################################
## LANDSCAPE GENOMICS ANALYSIS
## Mantel tests, RDA, and Genotype-Environment Association (GEA)
############################################################

# ============================================================
# 1. SETUP
# ============================================================

setwd("D:/Landscape_210")

# ---- Install required packages (run once, then comment out) ----
# install.packages("vegan")
# install.packages("geosphere")
# install.packages("vcfR")
# if (!require("devtools")) install.packages("devtools")
# devtools::install_github("thierrygosselin/radiator")
# devtools::install_github("thierrygosselin/assigner")
# install.packages("futile.matrix")
# install.packages("dartR")
# install.packages("PopGenome")
# install.packages("pegas")
# install.packages("corrplot")
# install.packages("robust")
# install_github("jdstorey/qvalue", force = TRUE)
# install.packages("qqman")
# install.packages("adegenet")
# install.packages("fit.models")
# BiocManager::install("qvalue")
# install.packages("remotes")
# remotes::install_github("koohyun-kwon/rdadapt")

# ---- Load libraries ----
library(fit.models)
library(adegenet)
library(vegan)
library(geosphere)
library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)
library(igraph)
library(philentropy)
library(radiator)
library(assigner)
library(ggplot2)
library(futile.matrix)
library(dartR)
library(PopGenome)
library(pegas)
library(corrplot)
library(robust)
library(qvalue)
library(qqman)
library(ggrepel)
library(ggpubr)
library(grid)
library(dplyr)


# ============================================================
# 2. IMPORT DATA
# ============================================================

# ---- Main dataset (climate + metadata) ----
# Option A: paste from clipboard (after Ctrl+C in Excel)
dataclim <- read.table(file = "clipboard", sep = "\t", header = TRUE)  # file: Landscape_genofile.xls

# Option B: read directly from the Excel file (preferred / reproducible)
library(readxl)
dataclim <- read_excel("D:/Landscape_210/landscape_genofile.xlsx")

# ---- Genetic data from VCF ----
genoLAND.VCF <- read.vcfR("EMCAP_533_SNPs_chr1-8_geolocalize_miss90_thinned10Kb.vcf.recode.vcf")
gl.genoLAND  <- vcfR2genind(genoLAND.VCF)          # convert to genind object
genotype     <- as.data.frame(gl.genoLAND)

# Genetic distance matrices
distgenEUCL <- dist(gl.genoLAND, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
distgenDISS <- diss.dist(gl.genoLAND, percent = FALSE, mat = FALSE)  # allelic differences between individuals

hist(distgenEUCL)

# Export a GenAlEx file (used to cross-check Mantel test results in GenAlEx)
genind2genalex(
  gl.genoLAND,
  filename  = "genealex_genoLand(1).xls",
  overwrite = FALSE,
  quiet     = FALSE,
  pop       = NULL,
  allstrata = TRUE,
  geo       = FALSE,
  geodf     = "xy",
  sep       = ",",
  sequence  = FALSE
)

# ---- Alternative genetic data source: allele counts from TASSEL5 ----
# Alleles are recoded as counts of the minor allele:
#   0 = homozygous major allele, 1 = heterozygote, 2 = homozygous alternative allele
# (loci-only file, no accession column: genoChickpeaLandscapebeagle.xlsx)
genotype <- read.table(file = "clipboard", sep = "\t", header = TRUE)
dist.geno <- dist(genotype, method = "euclidean")

# ---- Allele frequency dataset ----
AllFreq <- read.table(file = "clipboard", sep = "\t", header = TRUE)  # file: Sitefrequency


# ============================================================
# 3. DISTANCE MATRICES (BIOCLIM, GEOGRAPHIC)
# ============================================================

# ---- Bioclimatic (environmental) distance ----
PCbio       <- dataclim[, 19:37]
Env         <- scale(PCbio, center = TRUE, scale = TRUE)
dist.PCbio  <- dist(Env, method = "euclidean")

# ---- Geographic distance ----
geo      <- data.frame(dataclim$long, dataclim$lat)
dist.geo <- dist(geo, method = "euclidean")


# ============================================================
# 4. MANTEL TESTS
# ============================================================

# Shared plot theme used across the Mantel scatterplots below
mantel_theme <- theme(
  axis.text.x     = element_text(face = "bold", colour = "black", size = 18),
  axis.text.y     = element_text(face = "bold", size = 18, colour = "black"),
  axis.title      = element_text(face = "bold", size = 18, colour = "black"),
  panel.background = element_blank(),
  panel.border    = element_rect(fill = NA, colour = "black"),
  legend.title    = element_text(size = 12, face = "bold", colour = "black"),
  legend.text     = element_text(size = 10, face = "bold", colour = "black"),
  legend.position = "top",
  strip.background = element_rect(fill = "grey90", colour = "black"),
  strip.text      = element_text(size = 9, face = "bold")
)

# ---- 4.1 Genetic distance ~ geographic distance ----
geo_geno <- mantel(dist.geo, distgenEUCL, method = "spearman", permutations = 1000, na.rm = TRUE)
geo_geno
summary(lm(distgenEUCL ~ dist.geo))

graph <- mantel.correlog(distgenEUCL, dist.geo, XY = NULL, n.class = 0, break.pts = NULL,
                         cutoff = TRUE, r.type = "pearson", nperm = 999, mult = "holm", progressive = TRUE)
plot(graph)

xx <- as.vector(dist.geo)      # geographic distance
yy <- as.vector(distgenEUCL)   # genetic distance
manatelmatrix <- data.frame(xx, yy)

mm <- ggplot(manatelmatrix, aes(y = yy, x = xx)) +
  geom_point(size = 4, alpha = 0.75, colour = "black", shape = 21, fill = "grey") +
  geom_smooth(method = "lm", colour = "red", alpha = 0.2) +
  labs(y = "Euclidean genetic distance", x = "Euclidean geographic distance") +
  mantel_theme
mm

# ---- 4.2 Ecological distance ~ geographic distance ----
geo_eco <- mantel(dist.geo, dist.PCbio, method = "spearman", permutations = 1000, na.rm = TRUE)
geo_eco
summary(lm(dist.PCbio ~ dist.geo))

graph <- mantel.correlog(dist.PCbio, dist.geo, XY = NULL, n.class = 0, break.pts = NULL,
                         cutoff = TRUE, r.type = "pearson", nperm = 999, mult = "holm", progressive = TRUE)
plot(graph)

yy <- as.vector(dist.geo)      # geographic distance
zz <- as.vector(dist.PCbio)    # ecological distance
manatelmatrix <- data.frame(yy, zz)

mm <- ggplot(manatelmatrix, aes(y = yy, x = zz)) +
  geom_point(size = 4, alpha = 0.75, colour = "black", shape = 21, fill = "grey") +
  geom_smooth(method = "lm", colour = "red", alpha = 0.2) +
  labs(y = "Euclidean geographic distance", x = "Euclidean ecological distance") +
  mantel_theme
mm

# ---- 4.3 Genetic distance ~ ecological distance ----
geno_eco <- mantel(distgenEUCL, dist.PCbio, method = "spearman", permutations = 1000, na.rm = TRUE)
geno_eco
summary(lm(dist.PCbio ~ distgenEUCL))

graph <- mantel.correlog(distgenEUCL, dist.PCbio, XY = NULL, n.class = 0, break.pts = NULL,
                         cutoff = TRUE, r.type = "pearson", nperm = 999, mult = "holm", progressive = TRUE)
plot(graph)

xx <- as.vector(distgenEUCL)   # genetic distance
zz <- as.vector(dist.PCbio)    # ecological distance
manatelmatrix <- data.frame(zz, xx)

mm <- ggplot(manatelmatrix, aes(y = xx, x = zz)) +
  geom_point(size = 4, alpha = 0.75, colour = "black", shape = 21, fill = "grey") +
  geom_smooth(method = "lm", colour = "red", alpha = 0.2) +
  labs(y = "Euclidean genetic distance", x = "Euclidean ecological distance") +
  mantel_theme
mm

# ---- 4.4 Partial Mantel test: genetic ~ ecological | geographic ----
partial_mantel <- mantel.partial(distgenEUCL, dist.PCbio, dist.geo, method = "spearman",
                                 permutations = 1000, na.rm = TRUE)
partial_mantel
summary(lm(distgenEUCL ~ dist.PCbio | dist.geo))

xx <- as.vector(distgenEUCL)   # genetic distance
yy <- as.vector(dist.geo)      # geographic distance
zz <- as.vector(dist.PCbio)    # ecological distance
partial_mantel_matrix <- data.frame(xx, zz, yy)

mm <- ggplot(partial_mantel_matrix, aes(y = xx, x = zz)) +
  geom_point(size = 2.5, alpha = 0.75, colour = "black", shape = 21, aes(fill = yy)) +
  geom_smooth(method = "lm", colour = "red", alpha = 0.2) +
  scale_fill_continuous(high = "navy", low = "lightblue") +
  labs(y = "Euclidean genetic distance", x = "Euclidean ecological distance",
       fill = "Geographic distance") +
  mantel_theme
mm


# ============================================================
# 5. REDUNDANCY ANALYSIS (RDA)
# Reference: https://github.com/Capblancq/RDA-landscape-genomics
# ============================================================

# ---- Standardize climatic variables ----
ecobio <- dataclim[, 19:37]
Env    <- scale(ecobio, center = TRUE, scale = TRUE)
Env    <- as.data.frame(Env)

# ---- Neutral population structure (e.g. PCs from a structure/PCA analysis) ----
PopStruct <- dataclim[, 11:15]

# ---- Combine geographic, population structure, and environmental variables ----
Variables <- data.frame(dataclim$genovcf_code, geo, PopStruct, Env)
genotype  <- as.data.frame(gl.genoLAND)

# ---- Null and full models ----
RDA0     <- rda(genotype ~ 1, Variables)
RDAfull  <- rda(genotype ~ bio1 + bio2 + bio3 + bio4 + bio5 + bio6 + bio7 + bio8 + bio9 +
                  bio10 + bio11 + bio12 + bio13 + bio14 + bio15 + bio16 + bio17 + bio18 + bio19,
                Variables)

RDAtemp    <- rda(genotype ~ bio1 + bio2 + bio3 + bio4 + bio5 + bio6 + bio7 + bio8 + bio9 + bio10 + bio11, Variables)
RDAtempsel <- rda(genotype ~ bio8 + bio9, Variables)
RDAprec    <- rda(genotype ~ bio12 + bio13 + bio14 + bio15 + bio16 + bio17 + bio18 + bio19, Variables)
RDAprecsel <- rda(genotype ~ bio15 + bio18 + bio19, Variables)

# ---- Forward selection of variables ----
# Stopping criteria: p < 0.01 (1000 permutations) and adjusted R2 of the global model
mod <- ordiR2step(RDA0, RDAprec, Pin = 0.01, R2permutations = 1000, R2scope = TRUE)
mod$anova

# Check collinearity (VIF) among selected predictors
sqrt(vif.cca(RDAprec))
RDAtemp <- rda(genotype ~ bio7 + bio3 + bio8 + bio4 + bio6 + bio10 + bio11, Variables)
sqrt(vif.cca(RDAtempsel))
# Selected variables: bio3 & bio7 (temperature); bio17, bio18, bio19 (precipitation)


# ============================================================
# 6. VARIANCE PARTITIONING (PARTIAL RDA)
# Partitions variance among climate, population structure, and geography
# ============================================================

## Full model
pRDAfull <- rda(genotype ~ PC1 + PC2 + PC3 + dataclim.long + dataclim.lat +
                  bio8 + bio9 + bio15 + bio18 + bio19, Variables)
RsquareAdj(pRDAfull)
anova(pRDAfull)

## Pure climate model
pRDAclim <- rda(genotype ~ bio8 + bio9 + bio15 + bio18 + bio19 +
                  Condition(PC1 + PC2 + PC3 + dataclim.long + dataclim.lat), Variables)
RsquareAdj(pRDAclim)
anova.cca(pRDAclim)

## Pure neutral population structure model
pRDAstruct <- rda(genotype ~ PC1 + PC2 + PC3 +
                    Condition(dataclim.long + dataclim.lat + bio8 + bio9 + bio15 + bio18 + bio19), Variables)
RsquareAdj(pRDAstruct)
anova(pRDAstruct)

## Pure geography model
pRDAgeog <- rda(genotype ~ dataclim.long + dataclim.lat +
                  Condition(PC1 + PC2 + PC3 + bio8 + bio9 + bio15 + bio18 + bio19), Variables)
RsquareAdj(pRDAgeog)
anova(pRDAgeog)

## Geography + isolation-by-distance model
pRDAIBD <- rda(genotype ~ dataclim.long + dataclim.lat + bio8 + bio9 + bio15 + bio18, Variables)
RsquareAdj(pRDAIBD)
anova(pRDAgeog)  # NOTE: original script re-tests pRDAgeog here rather than pRDAIBD

# ---- Correlation matrix of bioclim + PC variables ----
corrplot(cor(dataclim[, c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9",
                          "bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19",
                          "PC1_prod","PC2_prod","PC1_flow","PC2_flow")]),
         type = "upper")


# ============================================================
# 7. RDA: GEOGRAPHY + ENVIRONMENT (COMBINED MODEL)
# ============================================================

colnames(Variables)[colnames(Variables) == "dataclim.lat"]  <- "Lat"
colnames(Variables)[colnames(Variables) == "dataclim.long"] <- "Long"

RDAgeo_env <- rda(genotype ~ Long + Lat + bio8 + bio9 + bio15 + bio18 + bio19, Variables)
summary(eigenvals(RDAgeo_env, model = "constrained"))

score <- scores(RDAgeo_env, display = "sites")
write.table(score, "Genotypevalue_RDAgeo_env")

# ---- Import site scores + genetic group assignments for plotting ----
data_RDAgeo_env <- read.csv("RDAgeo_env.csv", header = TRUE)

data_RDAgeo_env <- data_RDAgeo_env %>%
  mutate(
    GrK3 = recode(as.character(GrK3),
                  "1" = "Gr3.1", "2" = "Gr3.2", "3" = "Gr3.3", .default = GrK3),
    GrK6 = recode(as.character(GrK6),
                  "1" = "Gr6.1", "2" = "Gr6.2", "3" = "Gr6.3",
                  "4" = "Gr6.4", "5" = "Gr6.5", "6" = "Gr6.6", .default = GrK6)
  ) %>%
  mutate(
    GrK3 = factor(GrK3, levels = c("Gr3.1", "Gr3.2", "Gr3.3", "Admixed")),
    GrK6 = factor(GrK6, levels = c("Gr6.1", "Gr6.2", "Gr6.3", "Gr6.4", "Gr6.5", "Gr6.6", "Admixed"))
  )

TAB_gen <- data.frame(geno_names = row.names(score), score)

# Variable (biplot) loadings
TAB_var <- as.data.frame(scores(RDAgeo_env, choices = c(1, 2), display = "bp"))

# Common theme reused for the K3 / K6 RDA plots
common_theme <- theme_bw(base_size = 12, base_family = "Times") +
  theme(
    panel.background = element_blank(),
    legend.background = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank(),
    legend.text  = element_text(size = 11),
    legend.title = element_text(size = 12),
    axis.text    = element_text(size = 12),
    axis.title   = element_text(size = 13)
  )

# ---- 7.1 Plot: K3 genetic groups ----
loading_RDAgeo_env_K3 <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(.80), linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(.80), linewidth = 0.5) +
  geom_point(data = data_RDAgeo_env, aes(x = RDA1, y = RDA2, fill = GrK3),
             shape = 21, color = "black", size = 3.5, stroke = 0.6) +
  scale_fill_manual(values = c("blue", "darkorange", "chartreuse3", "darkgrey")) +
  geom_segment(data = TAB_var, aes(x = 0, y = 0, xend = RDA1 * 10, yend = RDA2 * 10),
               colour = "black", linewidth = 1, arrow = arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x = RDA1 * 10, y = RDA2 * 11, label = row.names(TAB_var)),
                   size = 4, family = "Times") +
  xlab("RDA 1 (58%)") + ylab("RDA 2 (19%)") +
  guides(fill = guide_legend(title = "Genetic group")) +
  common_theme

# ---- 7.2 Plot: K6 genetic groups ----
loading_RDAgeo_env_K6 <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(.80), linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(.80), linewidth = 0.5) +
  geom_point(data = data_RDAgeo_env, aes(x = RDA1, y = RDA2, fill = GrK6),
             shape = 21, color = "black", size = 3.5, stroke = 0.6) +
  scale_fill_manual(values = c("blue", "darkorange", "chartreuse3", "yellow", "brown", "purple", "darkgrey")) +
  geom_segment(data = TAB_var, aes(x = 0, y = 0, xend = RDA1 * 10, yend = RDA2 * 10),
               colour = "black", linewidth = 1, arrow = arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x = RDA1 * 10, y = RDA2 * 11, label = row.names(TAB_var)),
                   size = 4, family = "Times") +
  xlab("RDA 1 (58%)") + ylab("RDA 2 (19%)") +
  guides(fill = guide_legend(title = "Genetic group")) +
  common_theme

# ---- Combine K3 and K6 panels ----
final_plot <- ggarrange(
  loading_RDAgeo_env_K3, loading_RDAgeo_env_K6,
  nrow = 1, ncol = 2,
  labels = c("A", "B"),
  common.legend = FALSE,
  legend = "right",
  align = "hv"
)
final_plot

ggsave("RDA_geo_env_plot.tiff", plot = final_plot,
       width = 22, height = 8, units = "cm", dpi = 600, compression = "lzw")

anova(RDAgeo_env)


# ============================================================
# 8. GENOTYPE-ENVIRONMENT ASSOCIATION (GEA): IDENTIFYING LOCI UNDER SELECTION
# ============================================================

# Step 1: run RDA on genotype matrix with bioclim predictors,
# conditioning on the first 3 PCs to account for neutral population structure.
RDA_env  <- rda(genotype ~ bio8 + bio9 + bio15 + bio18 + bio19 +
                  Condition(gr1 + gr2 + gr3 + dataclim.long + dataclim.lat), Variables)
RDA_temp <- rda(genotype ~ bio9 + bio8 + Condition(PC1 + PC2 + PC3 + Long + Lat), Variables)
RDA_prec <- rda(genotype ~ bio15 + bio18 + bio19 + Condition(PC1 + PC2 + PC3 + Long + Lat), Variables)

plot(RDA_temp)
screeplot(RDA_prec, main = "Eigenvalues of constrained axes")
summary(eigenvals(RDA_temp, model = "constrained"))

# ---- rdadapt: computes Mahalanobis-distance-based p/q-values for RDA loadings ----
source("./src/rdadapt.R")
rdadapt <- function(rda, K) {
  zscores  <- rda$CCA$v[, 1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha  <- covRob(resscale, distance = TRUE, na.action = na.omit, estim = "pairwiseGK")$dist
  lambda   <- median(resmaha) / qchisq(0.5, df = K)
  reschi2test <- pchisq(resmaha / lambda, K, lower.tail = FALSE)
  qval <- qvalue(reschi2test)
  data.frame(p.values = reschi2test, q.values = qval$qvalues)
}

# ------------------------------------------------------------
# 8.1 PRECIPITATION
# ------------------------------------------------------------
rdadapt_env <- rdadapt(RDA_prec, 2)

# Bonferroni threshold
thres_env <- 0.05 / length(rdadapt_env$p.values)

# Loci below the Bonferroni threshold
top_outliers <- data.frame(
  Loci    = colnames(genotype)[which(rdadapt_env$p.values < thres_env)],
  p.value = rdadapt_env$p.values[which(rdadapt_env$p.values < thres_env)],
  contig  = unlist(lapply(strsplit(colnames(genotype)[which(rdadapt_env$p.values < thres_env)], split = "_"),
                          function(x) x[1]))
)
write.table(outliers, "Bonferroni_precipitation")

qvalue  <- data.frame(Loci = colnames(genotype), p.value = rdadapt_env$p.values, q.value = rdadapt_env$q.value)
outliers <- data.frame(
  Loci    = colnames(genotype)[which(rdadapt_env$q.values < 0.05)],
  p.value = rdadapt_env$p.values[which(rdadapt_env$q.values < 0.05)]
)

# Locus (species) scores for plotting
locus_scores <- scores(RDA_prec, choices = c(1:2), display = "species", scaling = "none")
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- "Not associated"
TAB_loci$type[TAB_loci$names %in% outliers$Loci]     <- "FDR"
TAB_loci$type[TAB_loci$names %in% top_outliers$Loci] <- "Bonferroni"
TAB_loci$type <- factor(TAB_loci$type, levels = c("Not associated", "FDR", "Bonferroni"))

TAB_var <- as.data.frame(scores(RDA_prec, choices = c(1, 2), display = "bp"))

loading_prec <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(.80), linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(.80), linewidth = 0.6) +
  geom_point(data = TAB_loci, aes(x = RDA1 * 40, y = RDA2 * 40, fill = type),
             shape = 21, color = "black", size = 3.5, stroke = 0.6) +
  scale_fill_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend = RDA1, yend = RDA2, x = 0, y = 0),
               colour = "black", linewidth = 1, linetype = 1, arrow = arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x = 1.1 * RDA1, y = 1.1 * RDA2, label = row.names(TAB_var)),
                   size = 4, family = "Times") +
  xlab("RDA 1: 48.8%") + ylab("RDA 2: 28.9%") +
  guides(color = guide_legend(title = "Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(),
        panel.grid = element_blank(), plot.background = element_blank(),
        legend.text = element_text(size = rel(.8)), strip.text = element_text(size = 11))
loading_prec

ggsave("RDA_prec.tiff", plot = loading_prec,
       width = 12, height = 8, units = "cm", dpi = 600, compression = "lzw")

# FDR results table
qvalue <- data.frame(Loci = colnames(genotype), p.value = rdadapt_env$p.values, q.value = rdadapt_env$q.value)
write.table(qvalue, "FDR_prec_3PCs")

loci_value <- data.frame(Loci = colnames(genotype), p.value = rdadapt_env$p.values)
write.table(loci_value, "GEA_RDA_prec.txt")

# ---- Manhattan plot for precipitation GEA results ----
Manhattan_prec <- read.table(file = "clipboard", sep = "\t", header = TRUE)  # p-values for precipitation

png(filename = "Manh_GEA_RDA_prec.png", width = 3000, height = 2000, res = 300)
manhattan(Manhattan_prec,
          col = c("darkblue", "gray60"),
          suggestiveline = -log10(0.00013416),
          genomewideline  = -log10(6.947918e-07),
          suggestivelinecol = "#6B4596FF",
          genomewidelinecol = "darkorange")
dev.off()

# Alternate version: plain Manhattan plot with threshold lines added manually
png(filename = "Manh_GEA_RDA_prec.png", width = 3000, height = 2000, res = 300)
manhattan(Manhattan_prec, col = c("darkblue", "gray60"), suggestiveline = FALSE, genomewideline = FALSE)
abline(h = -log10(0.00013416), col = "#6B4596FF", lwd = 2)
abline(h = -log10(6.947918e-07), col = "darkorange", lwd = 2)
dev.off()

# ---- RDA restricted to precipitation-enriched (significant) loci ----
geno_enrich <- genotype[which(rdadapt_env$q.values < 0.05)]
RDA_prec_enriched <- rda(geno_enrich ~ bio15 + bio18 + bio19 + Condition(PC1 + PC2 + PC3 + Long + Lat), Variables)
plot(RDA_prec_enriched)

locus_scores <- scores(RDA_prec_enriched, choices = c(1:2), display = "species", scaling = "none")
TAB_geno <- data.frame(genotype = row.names(TAB_gen), RDA_prec_enriched$CCA$u)
TAB_group_k3 <- data.frame(genotype = data_RDAgeo_env$genotype, group = data_RDAgeo_env$GrK3)
TAB_group_k6 <- data.frame(genotype = data_RDAgeo_env$genotype, group = data_RDAgeo_env$GrK6)
TAB_geno_group_k3 <- merge(TAB_geno, TAB_group_k3, by = "genotype")
TAB_geno_group_k6 <- merge(TAB_geno, TAB_group_k6, by = "genotype")
colnames(TAB_geno_group_k3) <- c("geno", "RDA1", "RDA2", "RDA3", "genetic_group")
colnames(TAB_geno_group_k6) <- c("geno", "RDA1", "RDA2", "RDA3", "genetic_group")

TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- factor("Adaptive loci", levels = c("Adaptive loci"))
TAB_var <- as.data.frame(scores(RDA_prec_enriched, choices = c(1, 2), display = "bp"))

loading_prec <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(.80), size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(.80), size = 0.6) +
  geom_point(data = TAB_geno, aes(x = RDA1, y = RDA2),
             shape = 21, fill = "darkblue", color = "black", size = 3.5, stroke = 0.6) +
  geom_segment(data = TAB_var, aes(xend = RDA1 / 5, yend = RDA2 / 5, x = 0, y = 0),
               colour = "black", size = 1, linetype = 1, arrow = arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x = RDA1 / 4.7, y = RDA2 / 4.7, label = row.names(TAB_var)),
                   size = 4, family = "Times") +
  xlab("RDA 1: 61.0 %") + ylab("RDA 2: 30.2 %") +
  guides(color = guide_legend(title = "Genetic group")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(),
        panel.grid = element_blank(), plot.background = element_blank(),
        legend.text = element_text(size = rel(.8)), strip.text = element_text(size = 11)) +
  ggtitle("RDA_prec") +
  theme(plot.title = element_text(color = "black", size = 14, face = "bold.italic"))
loading_prec

ggsave("enrichedRDA_prec.tiff", plot = loading_prec,
       width = 10, height = 8, units = "cm", dpi = 600, compression = "lzw")

# ---- Same plot, colored by genetic group (K3 / K6) ----
TAB_geno_group_k3 <- TAB_geno_group_k3 %>%
  mutate(genetic_group = factor(genetic_group, levels = c("Gr3.1", "Gr3.2", "Gr3.3", "Admixed")))

TAB_geno_group_k6 <- TAB_geno_group_k6 %>%
  mutate(genetic_group = factor(genetic_group,
                                levels = c("Gr6.1", "Gr6.2", "Gr6.3", "Gr6.4", "Gr6.5", "Gr6.6", "Admixed")))

loading_prec <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(.80), size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(.80), size = 0.6) +
  geom_point(data = TAB_geno_group_k3, aes(x = RDA1, y = RDA2, fill = genetic_group),
             shape = 21, color = "black", size = 3.5, stroke = 0.6) +
  scale_fill_manual(values = c("blue", "#F9A242FF", "limegreen", "darkgray")) +
  geom_segment(data = TAB_var, aes(xend = RDA1 / 5, yend = RDA2 / 5, x = 0, y = 0),
               colour = "black", size = 0.15, linetype = 1, arrow = arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x = RDA1 / 4.7, y = RDA2 / 4.7, label = row.names(TAB_var)),
                   size = 3.5, family = "Times") +
  xlab("RDA 1: 61.0 %") + ylab("RDA 2: 30.2 %") +
  guides(color = guide_legend(title = "Genetic group")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(),
        panel.grid = element_blank(), plot.background = element_blank(),
        legend.text = element_text(size = rel(.8)), strip.text = element_text(size = 11)) +
  ggtitle("RDA_prec") +
  theme(plot.title = element_text(color = "black", size = 14, face = "bold.italic"))
loading_prec

ggsave(loading_prec, filename = "RDA_prec_enriched(FDR).tiff", device = "tiff",
       limitsize = FALSE, dpi = 600, scale = 1.3)

plot(RDA_prec_enriched)
summary(eigenvals(RDA_prec_enriched, model = "constrained"))
write.table(RDA_prec_enriched$CCA$u, "Genotypevalue_RDA_prec_enriched_PC3")
anova(RDA_prec_enriched)

# ------------------------------------------------------------
# 8.2 TEMPERATURE
# ------------------------------------------------------------
rdadapt_env <- rdadapt(RDA_temp, 2)

thres_env <- 0.05 / length(rdadapt_env$p.values)

top_outliers <- data.frame(
  Loci    = colnames(genotype)[which(rdadapt_env$p.values < thres_env)],
  p.value = rdadapt_env$p.values[which(rdadapt_env$p.values < thres_env)],
  contig  = unlist(lapply(strsplit(colnames(genotype)[which(rdadapt_env$p.values < thres_env)], split = "_"),
                          function(x) x[1]))
)
write.table(outliers, "Bonferroni_temp")

qvalue <- data.frame(Loci = colnames(genotype), p.value = rdadapt_env$p.values, q.value = rdadapt_env$q.value)
outliers <- data.frame(
  Loci    = colnames(genotype)[which(rdadapt_env$q.values < 0.05)],
  p.value = rdadapt_env$p.values[which(rdadapt_env$q.values < 0.05)]
)

locus_scores <- scores(RDA_temp, choices = c(1:2), display = "species", scaling = "none")
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- "Not associated"
TAB_loci$type[TAB_loci$names %in% outliers$Loci]     <- "FDR"
TAB_loci$type[TAB_loci$names %in% top_outliers$Loci] <- "Bonferroni"
TAB_loci$type <- factor(TAB_loci$type, levels = c("Not associated", "FDR", "Bonferroni"))

TAB_var <- as.data.frame(scores(RDA_temp, choices = c(1, 2), display = "bp"))

loading_temp <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(.80), size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(.80), size = 0.6) +
  geom_point(data = TAB_loci, aes(x = RDA1 * 40, y = RDA2 * 40, fill = type),
             shape = 21, color = "black", size = 3.5, stroke = 0.6) +
  scale_fill_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend = RDA1, yend = RDA2, x = 0, y = 0),
               colour = "black", size = 1, linetype = 1, arrow = arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x = 1.1 * RDA1, y = 1.1 * RDA2, label = row.names(TAB_var)),
                   size = 4, family = "Times") +
  xlab("RDA 1: 69.7%") + ylab("RDA 2: 30.2%") +
  guides(color = guide_legend(title = "Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(),
        panel.grid = element_blank(), plot.background = element_blank(),
        legend.text = element_text(size = rel(.8)), strip.text = element_text(size = 11))
loading_temp

ggsave("RDA_temp.tiff", plot = loading_temp,
       width = 12, height = 8, units = "cm", dpi = 600, compression = "lzw")

qvalue <- data.frame(Loci = colnames(genotype), p.value = rdadapt_env$p.values, q.value = rdadapt_env$q.value)
write.table(qvalue, "FDR_temp_PC3")

# ---- Manhattan plot for temperature GEA results ----
Manhattan_temp <- read.table(file = "clipboard", sep = "\t", header = TRUE)  # p-values for temperature

manhattan(Manhattan_temp, col = c("darkred", "gray60"),
          suggestiveline = -log10(0.000129901), genomewideline = -log10(6.947918e-07))

# ---- RDA restricted to temperature-enriched (significant) loci ----
geno_enrich <- genotype[which(rdadapt_env$q.values < 0.05)]
RDA_temp_enriched <- rda(geno_enrich ~ bio8 + bio9 + Condition(PC1 + PC2 + PC3 + Long + Lat), Variables)

TAB_geno <- data.frame(RDA_temp_enriched$CCA$u)
locus_scores <- scores(RDA_temp_enriched, choices = c(1:2), display = "species", scaling = "none")
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- factor("Adaptive loci", levels = c("Adaptive loci"))
TAB_var <- as.data.frame(scores(RDA_temp_enriched, choices = c(1, 2), display = "bp"))

loading_temp <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(.80), size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(.80), size = 0.6) +
  geom_point(data = TAB_geno, aes(x = RDA1, y = RDA2),
             shape = 21, fill = "darkred", color = "black", size = 3.5, stroke = 0.6) +
  geom_segment(data = TAB_var, aes(xend = RDA1 / 5, yend = RDA2 / 5, x = 0, y = 0),
               colour = "black", size = 1, linetype = 1, arrow = arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x = RDA1 / 4.7, y = RDA2 / 4.7, label = row.names(TAB_var)),
                   size = 4, family = "Times") +
  xlab("RDA 1: 55.2 %") + ylab("RDA 2: 44.8 %") +
  guides(color = guide_legend(title = "Genetic group")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(),
        panel.grid = element_blank(), plot.background = element_blank(),
        legend.text = element_text(size = rel(.8)), strip.text = element_text(size = 11)) +
  ggtitle("RDA_temp") +
  theme(plot.title = element_text(color = "black", size = 14, face = "bold.italic"))
loading_temp

ggsave("enrichedRDA_temp.tiff", plot = loading_temp,
       width = 10, height = 8, units = "cm", dpi = 600, compression = "lzw")
# NOTE: original script saves `loading_prec` here (likely intended loading_temp)
ggsave(loading_prec, filename = "RDA_temp_enriched.tiff", device = "tiff",
       limitsize = FALSE, dpi = 600, scale = 1.3)

plot(RDA_temp_enriched)
summary(eigenvals(RDA_temp_enriched, model = "constrained"))
write.table(RDA_temp_enriched$CCA$u, "Genotypevalue_RDA_temp_enriched_PC3")
anova(RDA_temp_enriched)

# ---- Combined temperature + precipitation enriched RDA panel ----
ggarrange(loading_temp, loading_prec, nrow = 1, ncol = 2)


# ============================================================
# 9. TRAIT-RDA RELATIONSHIPS (FLOWERING / DROUGHT) & FINAL PLOTS
# ============================================================

df <- read.table(file = "clipboard", sep = "\t", header = TRUE)  # file: Landscape_genofile

df <- df %>%
  mutate(Kgroup_K3 = factor(Kgroup_K3, levels = c("Gr3.1", "Gr3.2", "Gr3.3", "Admixed")))

# Simple linear relationship between two log-transformed PCs
model <- lm(log_PC5 ~ log_PC3, data = df)
summary(model)
plot(model)
abline(model)
summary(model)$r.squared

boxplot(df$PC1_flow ~ df$adapt, data = df)
plot(df$RDA2_prec ~ df$PC1_flow)

# ---- PC1 flowering vs RDA1_temp, colored by K3 genetic group ----
B_K3 <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(.80), linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(.80), linewidth = 0.6) +
  geom_point(data = df, aes(x = PC1_flow, y = RDA1_temp, fill = Kgroup_K3),
             shape = 21, color = "black", size = 3.5, stroke = 0.6) +
  scale_fill_manual(values = c("blue", "darkorange", "chartreuse3", "darkgrey")) +
  xlab("PC1 flowering") + ylab("RDA1_temp") +
  guides(fill = guide_legend(title = "Genetic group")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(),
        panel.grid = element_blank(), plot.background = element_blank(),
        legend.text = element_text(size = rel(.8)), strip.text = element_text(size = 11))

# ---- PC1 flowering vs RDA2_prec, colored by K3 genetic group ----
C_K3 <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(.80), linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(.80), linewidth = 0.6) +
  geom_point(data = df, aes(x = PC1_flow, y = RDA2_prec, fill = Kgroup_K3),
             shape = 21, color = "black", size = 3.5, stroke = 0.6) +
  scale_fill_manual(values = c("blue", "darkorange", "chartreuse3", "darkgrey")) +
  xlab("PC1 flowering") + ylab("RDA2_prec") +
  guides(fill = guide_legend(title = "Genetic group")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(),
        panel.grid = element_blank(), plot.background = element_blank(),
        legend.text = element_text(size = rel(.8)), strip.text = element_text(size = 11))

ggarrange(B_K3, C_K3, nrow = 1, ncol = 2)

final_plot <- ggarrange(B_K3, C_K3, nrow = 1, ncol = 2,
                        common.legend = FALSE, legend = "right", align = "hv")
final_plot

ggsave("final_rda_floweirng_k3.tiff", plot = final_plot,
       width = 22, height = 7, units = "cm", dpi = 600, compression = "lzw")

# ---- Same two plots, colored by K6 genetic group ----
B_K6 <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(.80), linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(.80), linewidth = 0.6) +
  geom_point(data = df, aes(x = PC1_flow, y = RDA1_temp, fill = Kgroup_K6),
             shape = 21, color = "black", size = 2.5, stroke = 0.6) +
  scale_fill_manual(values = c("blue", "darkorange", "chartreuse3", "yellow", "brown", "purple", "darkgrey")) +
  xlab("PC1 flowering") + ylab("RDA1_temp") +
  guides(fill = guide_legend(title = "Genetic group")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(),
        panel.grid = element_blank(), plot.background = element_blank(),
        legend.text = element_text(size = rel(.8)), strip.text = element_text(size = 11))

C_K6 <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(.80), linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(.80), linewidth = 0.6) +
  geom_point(data = df, aes(x = PC1_flow, y = RDA2_prec, fill = Kgroup_K6),
             shape = 21, color = "black", size = 2.5, stroke = 0.6) +
  scale_fill_manual(values = c("blue", "darkorange", "chartreuse3", "yellow", "brown", "purple", "darkgrey")) +
  xlab("PC1 flowering") + ylab("RDA2_prec") +
  guides(fill = guide_legend(title = "Genetic group")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(),
        panel.grid = element_blank(), plot.background = element_blank(),
        legend.text = element_text(size = rel(.8)), strip.text = element_text(size = 11))

ggarrange(B_K6, C_K6, nrow = 1, ncol = 2)


# ============================================================
# 10. CHI-SQUARE TEST (e.g. early flowering ~ drought group)
# ============================================================

dt <- as.vector(df)
chisq <- chisq.test(dt$early, dt$drough, simulate.p.value = TRUE, B = 1000)
chisq
chisq$p.value