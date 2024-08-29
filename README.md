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
#GEO-POPSTRUCTURE
pRDAIBD <- rda(genotype ~ dataclim.long + dataclim.lat + bio8 + bio9 + bio15 + bio18, Variables)
RsquareAdj(pRDAIBD)
anova(pRDAgeog)
```

