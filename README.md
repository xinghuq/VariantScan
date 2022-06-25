
# VariantScan: a machine leaning tool for variant association testing



# Introduction 

This package provides a set of tool for performing association tests to identify QTLs in genome-wide association studies (GWAS,MWAS,EWAS,PWAS). This packages integrates three methods, Linear Model, Local Polynomial Fitting (Nonlinear Model) 
and Generalized Additive Model (GAM) to perform association testing in genome wide scan studies. 
This packge also applies to case-control studies, where the ROC is used to access the model performance.

Welcome any [feedback](https://github.com/xinghuq/DA/issues) and [pull request](https://github.com/xinghuq/DA/pulls).  

# Installation

## Install the package from github:
```{R}
library(devtools)

install_github("xinghuq/VariantScan")

library("VariantScan")
```

## Testing the association between phenotypes and genotypes using genomic data

Get example file

```{R}
f <- system.file('extdata',package='VariantScan')
infile <- file.path(f, "sim1.csv")
## read genotype file
geno=read.csv(infile)

# traits
traitq=geno[,14]
genotype=geno[,-c(1:14)]

# get PCs as covariates

PCs=prcomp(genotype)
PCs$x[,1:2]

## do Vscan using local polynomial regression fitting without specifying covariates

loessW=VScan(x=genotype,y=(traitq),methods ="loess")

## do Vscan using local polynomial regression fitting using PCs as covariates

loessWcv=VScan(x=genotype,y=(traitq),U=PCs$x[,1:2],methods ="loess")

## try linear model

lmW=VScan(x=genotype,y=(traitq),methods ="lm")
lmWcv=VScan(x=genotype,y=(traitq),U=PCs$x[,1:2],methods ="lm")

``````

## Visualizing the association signatures 

Plot Manhattan plot

``````
## 
Loci<-rep("Neutral", 1000)
Loci[c(201,211,221,231,241,251,261,271,281,291)]<-"QT"
Selected_Loci<-Loci[-which(Loci=="Neutral")]

library(ggplot2)
## Manhattan plot

g1=ggplot() +
  geom_point(aes(x=which(Loci=="Neutral"), y=-log10(lmWcv$p_norm$p.value[-which(Loci!="Neutral")])), col = "gray83") +
  geom_point(aes(x=which(Loci!="Neutral"), y=-log10(lmWcv$p_norm$p.value[-which(Loci=="Neutral")]), colour = Selected_Loci)) +
  xlab("SNPs") + ylab("-log10(p-value)") +ylim(c(0,35))+theme_bw()

g1


g2=ggplot() +
  geom_point(aes(x=which(Loci=="Neutral"), y=-log10(loessWcv$p_norm$p.value[-which(Loci!="Neutral")])), col = "gray83") +
  geom_point(aes(x=which(Loci!="Neutral"), y=-log10(loessWcv$p_norm$p.value[-which(Loci=="Neutral")]), colour = Selected_Loci)) +xlab("SNPs") + ylab("-log10(p-value)") +ylim(c(0,35))+theme_bw()

g3


```

