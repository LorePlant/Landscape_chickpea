# Landscape genomics highlights the adaptive evolution of chickpea across the Silks Roads 
In this README we are going to see the main landscape genomics analysis flow used in Rocchetti et al., 2025

- [**1. Data input**](https://github.com/NuCicer/AlphaSimR/blob/main/README.md#1change-in-the-additive-genetic-variance)
- [**2. Mantel Test**](https://github.com/NuCicer/AlphaSimR/blob/main/README.md#2-lenght-of-the-breeding-cycle)
- [**3. Redundancy Analysis RDA**](https://github.com/NuCicer/AlphaSimR/blob/main/README.md#3-trait-accuracy)

## Data input 
> Environmental (bioclimatic data from WordClim), geographic and genetic data input
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


