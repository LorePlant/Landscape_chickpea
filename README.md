# Landscape genomics highlights the adaptive evolution of chickpea across the Silks Roads 
In this README we are going to see the main landscape genomics analysis flow used in Rocchetti et al., 2025

- [**1. Data input**](https://github.com/NuCicer/AlphaSimR/blob/main/README.md#1change-in-the-additive-genetic-variance)
- [**2. Mantel Test**](https://github.com/NuCicer/AlphaSimR/blob/main/README.md#2-lenght-of-the-breeding-cycle)
- [**3. Redundancy Analysis RDA**](https://github.com/NuCicer/AlphaSimR/blob/main/README.md#3-trait-accuracy)

## Data input 
Environmental (bioclimatic data from WordClim), geographic and genetic data input and related distance calcolation (Euclidean distance)

```
#bioclim PCdata frame
PCbio = dataclim[,19:37]
Env <- scale(PCbio, center=TRUE, scale=TRUE)
dist.PCbio = dist(Env, method = "euclidean")

#geographic data
geo = data.frame(dataclim$long, dataclim$lat)
dist.geo = dist(geo, method = "euclidean")

#genetic  from VCF 
genoLAND.VCF <- read.vcfR("EMCAP_533_SNPs_chr1-8_geolocalize_miss90_thinned10Kb.vcf.recode.vcf")#import vcf file
gl.genoLAND <- vcfR2genind(genoLAND.VCF)#transfrom file in genind object
genotype<-as.data.frame(gl.genoLAND)
```
## **2. Mantel test**
We used Mantel test to asses the linear relationship between genetic distance with geographic distance and genetic distance with environmental distance. Ultimately we used partial Mantel test to asses the linear relationship between genetic distance matrix and environmental distance matrix considering the geopgraphic matrix as covariate.

> Mantel test genetic vs geography
```
geo_geno = mantel(dist.geo, distgenEUCL, method="spearman", permutations=1000,  na.rm = TRUE)
geo_geno
summary(lm(distgenEUCL~dist.geo))
graph = mantel.correlog(distgenEUCL, dist.geo, XY=NULL, n.class=0, break.pts=NULL, 
                        cutoff=TRUE, r.type="pearson", nperm=999, mult="holm", progressive=TRUE)

plot(graph)
#plotting Mantel test genetic and geographic distance
xx = as.vector(dist.geo) #convert distance matrix into a vector
yy= as.vector(distgenEUCL)
manatelmatrix = data.frame(xx,yy)
mm = ggplot(manatelmatrix, aes(y = yy, x = xx)) + 
  geom_point(size = 4, alpha = 0.75, colour = "black",shape = 21,fill = "grey") + 
  geom_smooth(method = "lm", colour = "red", alpha = 0.2)+
  labs(y = "Eucledian genetic distance", x = "Eucledian geographic distance")+
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 18), 
         axis.text.y = element_text(face = "bold", size = 18, colour = "black"), 
         axis.title= element_text(face = "bold", size = 18, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"), 
         legend.title = element_text(size =12, face = "bold", colour = "black"),
         legend.text = element_text(size = 10, face = "bold", colour = "black"), 
         legend.position = "top", strip.background = element_rect(fill = "grey90", colour = "black"),
         strip.text = element_text(size = 9, face = "bold"))

mm
```
> Mantel test genetic vs environment
```
 mantel test Genetic distance-ecological distance
geno_eco = mantel(distgenEUCL, dist.PCbio, method="spearman", permutations=1000,  na.rm = TRUE)
geno_eco
summary(lm(dist.PCbio~distgenEUCL))
graph = mantel.correlog(distgenEUCL, dist.PCbio, XY=NULL, n.class=0, break.pts=NULL, 
                        cutoff=TRUE, r.type="pearson", nperm=999, mult="holm", progressive=TRUE)

plot(graph)

xx = as.vector(distgenEUCL) #convert distance matrix into a vector
zz = as.vector(dist.PCbio)
manatelmatrix = data.frame(zz,xx)
mm = ggplot(manatelmatrix, aes(y = xx, x = zz))+
  geom_point(size = 4, alpha = 0.75, colour = "black",shape = 21,fill = "grey") + 
  geom_smooth(method = "lm", colour = "red", alpha = 0.2)+
  labs(y = "Euclidean genetic distance", x = "Euclidean ecological distance")+
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 18), 
         axis.text.y = element_text(face = "bold", size = 18, colour = "black"), 
         axis.title= element_text(face = "bold", size = 18, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"), 
         legend.title = element_text(size =12, face = "bold", colour = "black"),
         legend.text = element_text(size = 10, face = "bold", colour = "black"), 
         legend.position = "top", strip.background = element_rect(fill = "grey90", colour = "black"),
         strip.text = element_text(size = 9, face = "bold"))

mm
```
> Partial Mantel test genetic vs environment considering geography as covariate
```
#partial Mantel test 
partial_mantel = mantel.partial(distgenEUCL, dist.PCbio, dist.geo, method = "spearman", permutations = 1000,
                                na.rm = TRUE)
partial_mantel
summary(lm(distgenEUCL~dist.PCbio|dist.geo))
#plotting partial Mantel test
xx = as.vector(distgenEUCL) #convert distance matrix into a vector
yy= as.vector(dist.geo)
zz = as.vector(dist.PCbio)
partial_mantel_matrix = data.frame(xx,zz,yy)#new data frame with vectorize distance matrix

mm = ggplot(partial_mantel_matrix, aes(y = xx, x = zz)) + 
  geom_point(size = 2.5, alpha = 0.75, colour = "black",shape = 21, aes(fill = yy)) + 
  geom_smooth(method = "lm", colour = "red", alpha = 0.2) + 
  scale_fill_continuous(high = "navy", low = "lightblue")+
  labs(y = "Eucledian genetic distance", x = "Eucledian ecological distance", fill= "geographic distance")+
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 18), 
         axis.text.y = element_text(face = "bold", size = 18, colour = "black"), 
         axis.title= element_text(face = "bold", size = 18, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"), 
         legend.title = element_text(size =12, face = "bold", colour = "black"),
         legend.text = element_text(size = 10, face = "bold", colour = "black"), 
         legend.position = "top", strip.background = element_rect(fill = "grey90", colour = "black"),
         strip.text = element_text(size = 9, face = "bold"))


mm
```

## **3. Redundacy Analysis**

Redundancy analysis was used as landscape genomic tool to dissect the variance component given by environmental variable, geography and population structure on the total genetic variance present in the study population. RDA was ultimetly used for Genotype Environment Association (GEA) identifing SNP marker associated with the multivariate space of environmental variables.
We followed the approach described in Thibaut Capblancq & Brenna Forester (2021) https://github.com/Capblancq/RDA-landscape-genomics

To start out we scaled the environmental variable and create a uniaue dataset of bioclim environmental variable, latitude and longitude and principal componets for population structure.
```
#standardize climatic variable to ensure that the variable units are comparable
ecobio = dataclim[,19:37]
Env <- scale(ecobio, center=TRUE, scale=TRUE)
Env <- as.data.frame(Env)
## Neutral population structure table
PopStruct = dataclim[,11:15]
#combining geographic, Popstructure, environmental (scaled) variables
Variables <- data.frame(dataclim$genovcf_code, geo, PopStruct, Env)
#dataframe genotypic data
genotype<-as.data.frame(gl.genoLAND)
```

We run the RDA seperatly for temperature and precipitation bioclimatic variable. To reduce collinearity in the two groups we used Variance Inflation Factor (VIF) remouving variables with VIF >10

```
## Full model
RDAfull<- rda(genotype ~ bio1 + bio2 + bio3 + bio4 + bio5 + bio6 + bio7 + bio8 + bio9 + bio10 + bio11 + bio12 + bio13 + bio14 + bio15 + bio16 + bio17 + bio18 + bio19, Variables)
RDAtemp<-rda(genotype ~ bio1 + bio2 + bio3 + bio4 + bio5 + bio6 + bio7 + bio8 + bio9 + bio10 + bio11 , Variables )
sqrt(vif.cca(RDAtemp))
RDAtempsel<-rda(genotype ~  bio8 + bio9, Variables)
RDAprec<-rda(genotype ~ bio12 + bio13 + bio14 + bio15 + bio16 + bio17 + bio18 + bio19 , Variables )
sqrt(vif.cca(RDAprec))
RDAprecsel<-rda(genotype ~ bio15 + bio18 + bio19 , Variables )
```

RDA combine linear regression with pricipal component analysis. In Landscape genomics it represent a efficient tool to dissect the single effect of environment, spatial and demographic effect on totalt genetic diversity while considering the other as covariate. 

```
##full model
pRDAfull <- rda(genotype ~ PC1 + PC2 + PC3 + dataclim.long + dataclim.lat + bio8 + bio9 + bio15 + bio18 + bio19, Variables)
RsquareAdj(pRDAfull)
anova(pRDAfull)
## Pure climate model
pRDAclim <- rda(genotype ~ bio8 + bio9 + bio15 + bio18 + bio19 + Condition(PC1 + PC2 + PC3 + dataclim.long + dataclim.lat), Variables)
RsquareAdj(pRDAclim)
anova.cca(pRDAclim)
## Pure neutral population structure model  
pRDAstruct <- rda(genotype ~ PC1 + PC2 + PC3 + Condition(dataclim.long + dataclim.lat + bio8 + bio9 + bio15 + bio18 + bio19), Variables)
RsquareAdj(pRDAstruct)
anova(pRDAstruct)
##Pure geography model
pRDAgeog <- rda(genotype ~ dataclim.long + dataclim.lat + Condition(PC1 + PC2 + PC3 +bio8 + bio9 + bio15 + bio18 + bio19), Variables)
RsquareAdj(pRDAgeog)
anova(pRDAgeog)
```

To visualize the environment and geographic effect on differentianting genetic groups we use RDA considering only environment and geographic variable. Results are represented in the RDA biplot

```
#GEO-POPSTRUCTURE
pRDAIBD <- rda(genotype ~ dataclim.long + dataclim.lat + bio8 + bio9 + bio15 + bio18, Variables)
RsquareAdj(pRDAIBD)
anova(pRDAgeog)

#draw partial RDA geo+env
colnames(Variables)[colnames(Variables) == 'dataclim.lat'] <- 'Lat'
colnames(Variables)[colnames(Variables) == 'dataclim.long'] <- 'Long'
RDAgeo_env <- rda(genotype ~ Long + Lat + bio8 + bio9 + bio15 + bio18 + bio19, Variables)
summary(eigenvals(RDAgeo_env, model = "constrained"))
score<-scores(RDAgeo_env , display = "sites")
write.table(score, "Genotypevalue_RDAgeo_env")
data_RDAgeo_env <- read.table(file = "clipboard", 
                              sep = "\t", header=TRUE)
TAB_gen <- data.frame(geno_names = row.names(score), score)
install.packages("ggrepel")
library(ggrepel)
TAB_var <- as.data.frame(scores(RDAgeo_env, choices=c(1,2), display="bp"))
loading_RDAgeo_env_K3<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_point(data = data_RDAgeo_env, aes(x=RDA1, y=RDA2, color=GrK3), size = 4.5) +
  scale_color_manual(values = c("blue", "darkorange", "chartreuse3", "darkgrey")) + 
  geom_segment(data = TAB_var, aes(xend=RDA1*10, yend=RDA2*10, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x=RDA1*10, y=RDA2*11, label = row.names(TAB_var)), size = 4.5, family = "Times") +
  xlab("RDA 1: 58 %") + ylab("RDA 2: 19 %") +
  guides(color=guide_legend(title="Genetic group")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(1)), strip.text = element_text(size=13),axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))
loading_RDAgeo_env_K3

#same plot considering genetic groups at K=6
loading_RDAgeo_env_K6<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = data_RDAgeo_env, aes(x=RDA1, y=RDA2, color=GrK6), size = 4.5) +
  scale_color_manual(values = c("blue", "darkorange", "chartreuse3","yellow", "brown", "purple", "darkgrey")) + 
  geom_segment(data = TAB_var, aes(xend=RDA1*10, yend=RDA2*10, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x=RDA1*10, y=RDA2*11, label = row.names(TAB_var)), size = 4.5, family = "Times") +
  xlab("RDA 1: 58 %") + ylab("RDA 2: 19 %") +
  guides(color=guide_legend(title="Genetic group")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(1)), strip.text = element_text(size=13),axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))
loading_RDAgeo_env_K6

ggarrange(loading_RDAgeo_env_K3, loading_RDAgeo_env_K6, nrow=1, ncol=2)
```
We are going to use the RDA for Genotype Environment Association (GEA) to identify loci associated with multivariate environmental space. We run the analysis for the two groups of temperature and precipitation variables considering the geographic (latitude and longitude) and demographic (first three PCs from genomic data) as covariates on both of them. The procedure illustrated in Thibaut Capblancq & Brenna Forester (2021) is to calcolate the Mahalanobis distance of SNP markers from the origin in the RDA space. The Mahalanobis distance accounting for the different variance on the two RDA axis allow to correctely rank the SNP distqnces from the RDA origin. From the distribution of distances for all SNPs Pvalue and qvalus are determinated using the function _rdadapt_
```
#Genotype-Environment Associations: identifying loci under selection
##first step run a RDA model on the genotypic data matrix using all the normalized bioclim variable as explanotory variable and the first 3 PCs as conditioning variable to account for netrual pop structure 
RDA_env<- rda(genotype ~  bio8 + bio9 + bio15 + bio18 + bio19 + Condition(gr1 + gr2 + gr3 + dataclim.long + dataclim.lat), Variables )
RDA_temp<-rda(genotype ~ bio9 + bio8 + Condition(PC1 + PC2 + PC3 + dataclim.long + dataclim.lat), Variables )
RDA_prec<-rda(genotype ~  bio15 + bio18 + bio19 + Condition(PC1 + PC2 + PC3 + dataclim.long + dataclim.lat), Variables )
plot(RDA_temp)
screeplot(RDA_prec, main="Eigenvalues of constrained axes")
summary(eigenvals(RDA_temp, model = "constrained"))
#function radapt for FDR p-values and q-values
install.packages("remotes")
remotes::install_github("koohyun-kwon/rdadapt")
source("./src/rdadapt.R")
rdadapt<-function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}
```
We start the analysis for the precipitation variables applying two thresholds to identifying highly associated SNPs using a Bonferroni correction (Pvqlue = 0.05/numb. comparison) and less associated SNP using a False Discovery Rate (FDR) q values<0.05. The results are illustrated in the RDA biplot.

```
#precipitation
rdadapt_env<- rdadapt(RDA_prec, 2)
## P-values threshold after Bonferroni correction
thres_env <- 0.05/length(rdadapt_env$p.values)
## Identifying the loci that are below the p-value threshold
top_outliers <- data.frame(Loci = colnames(genotype)[which(rdadapt_env$p.values<thres_env)], p.value = rdadapt_env$p.values[which(rdadapt_env$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(genotype)[which(rdadapt_env$p.values<thres_env)], split = "_"), function(x) x[1])))
write.table(outliers, "Bonferroni_precipitation")
qvalue <- data.frame(Loci = colnames(genotype), p.value = rdadapt_env$p.values, q.value = rdadapt_env$q.value)
outliers <- data.frame(Loci = colnames(genotype)[which(rdadapt_env$q.values<0.05)], p.value = rdadapt_env$p.values[which(rdadapt_env$q.values<0.05)])
locus_scores <- scores(RDA_prec, choices=c(1:2), display="species", scaling="none")
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- "Not associated"
TAB_loci$type[TAB_loci$names%in%outliers$Loci] <- "FDR"
TAB_loci$type[TAB_loci$names%in%top_outliers$Loci] <- "Bonferroni"
TAB_loci$type <- factor(TAB_loci$type, levels = c("Not associated", "FDR", "Bonferroni"))
TAB_var <- as.data.frame(scores(RDA_prec, choices=c(1,2), display="bp"))
loading_prec<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*40, y=RDA2*40, colour = type), size = 2.5) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 1: 48.8%") + ylab("RDA 2: 28.9%") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
loading_prec
```
The same results can be illustrated by entering the genome position of each SNP to plot a Manhattan plot. We entered FDR and Bonferroni thresholds at the respective logP values -log10(0.00013416)and -log10(6.947918e-07). For major detail see the package _qqman_

```
#plotting Mhanattan plot using the library qqman
Manhattan_prec <- read.table(file = "clipboard", 
                             sep = "\t", header=TRUE) #import the p value result for precipitation


a<-manhattan(Manhattan_prec, col = c("darkblue", "gray60"),suggestiveline = -log10(0.00013416), genomewideline = -log10(6.947918e-07))
```
