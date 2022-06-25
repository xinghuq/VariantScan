\name{VScan}
\alias{VScan}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Variants (Biomarkers, e.g., genomic (genetic loci), transcriptomic (gene expression), epigenomic (methylations), proteomic(protein), metabolomic (metabolites) variants) Scanning and Association Tests Using Local Polynomial Fitting (Nonlinear Model) or Linear Model
%%  ~~function to do ... ~~
}
\description{Performing association tests for QTLs in genome-wide association studies (GWAS,MWAS,EWAS,PWAS) using nonlinear (Local Polynomial Fitting) or liner model.  When a nonlinear model loess smoother is selected, the intercept is used as a null model to test and calculate the R square statistic, then the R square is used as the weight for estimating the variant effect. In linear model, the beta value of the linear model (lm) is used as the weight for estimating the effect size of a variant. This function also applies to case-control studies, where the ROC is used to access the model performance.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VScan(x, y, U=NULL,methods = "loess",span = 0.65, family="gaussian", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Variant matrix, can be genomic (genetic loci), transcriptomic (gene expression), epigenomic (methylations), proteomic(protein), metabolomic (metabolites) variants.
%%     ~~Describe \code{x} here~~
}
  \item{y}{Traits
%%     ~~Describe \code{y} here~~
}
  \item{U}{Covariates,confounding factors.e.g.,  age, sex, PCs.
%%     ~~Describe \code{y} here~~
}
  \item{methods}{Model fit methods, whether "loess" or "lm". If method is "loess", the model will fit a local polynomial regression, if method is "lm", model will fit a linear model.
%%     ~~Describe \code{methods} here~~
}
  \item{span}{ The local polynomial regression parameter alpha which controls the degree of smoothing
%%     ~~Describe \code{span} here~~
}
  \item{family}{ The local polynomial regression parameter. if "gaussian", the fitting is done by least-squares, and if "symmetric" a re-descending M estimator is used with Tukey's biweight function. Can be abbreviated.
%%     ~~Describe \code{family} here~~
}
  \item{\dots}{ More arguments and parameters passing to loess and loess.control.
%%     ~~Describe \code{\dots} here~~
}
}
\details{ Fast association testing and variant scanning (GWAS,MWAS,EWAS,PWAS) using nonlinear (local polynomial fitting, loess) or liner model.When a nonlinear model loess smoother is selected, the intercept is used as a null model to test and calculate the R square statistic, then this R square is used as the weight for estimating the variant effect. In linear model, the beta value of the linear model (lm) is used as the weight for estimating the effect size of a variant.

Biomarkers can be biallelic loci, gene expression, methylations or protein expression. Biallelic markers have only two alleles, in GWAS, genotypes or the reference allele coypes are usually used to test the assocation between phenotypes and genotypes. Most of the cases, linear model is poweful enough to approximate the variant effect. However, if quantitive traits correspond to dynamic gene/protein expression data, and the associations may not be adequately approximated by simple linear model. 

Local polynomial regression (LOESS) build on "classical" methods, such as linear and nonlinear least squares regression. They address situations in which the classical procedures do not perform well or cannot be effectively applied without undue labor. LOESS combines much of the simplicity of linear least squares regression with the flexibility of nonlinear regression. It does this by fitting simple models to localized subsets of the data to build up a function that describes the deterministic part of the variation in the data, point by point.

When fitting a local polynomial regression, the model is fitted locally. For the fit at point x, the fit is made using points in a neighbourhood of x, weighted by their distance from x (with differences in "parametric" variables being ignored when computing the distance). The size of the neighbourhood is controlled by alpha (set by span or enp.target). For alpha < 1, the neighbourhood includes proportion alpha of the points, and these have tricubic weighting (proportional to (1 - (dist/maxdist)^3)^3). For alpha > 1, all points are used, with the "maximum distance"" assumed to be alpha^(1/p) times the actual maximum distance for p explanatory variables.

For the default family, fitting is by (weighted) least squares. For family="symmetric" a few iterations of an M-estimation procedure with Tukey's biweight are used. Be aware that as the initial value is the least-squares fit, this need not be a very resistant fit.

%%  ~~ If necessary, more details than the description above ~~
}
\value{Thus function gives the weights and p-values of various QTLs
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{W }{W, the weigts of the variants}
\item{p_norm}{p values assuming the normal distribution. The W is converted by arcsine transfromaton to normal distribution before estimating the p values }
\item{pvalue_chi}{p values assuming the Chisquare distribution with df= 1. The W is converted by arcsine transfromaton to normal distribution before estimating the p values }
%% ...
}
\references{
Opsomer, J. D., Ruppert, D. (1997). Fitting a bivariate additive model by local polynomial regression. The Annals of Statistics, 25(1), 186-211.


W. S. Cleveland, E. Grosse and W. M. Shyu (1992) Local regression models. Chapter 8 of Statistical Models in S eds J.M. Chambers and T.J. Hastie, Wadsworth & Brooks/Cole.
%% ~put references to the literature/web site here ~
}
\author{qinxinghu@gmail.com
%%  ~~who you are~~
}

\examples{
# require time to run

f <- system.file('extdata',package='VariantScan')
infile <- file.path(f, "sim1.csv")
geno=read.csv(infile)

traitq=geno[,14]
genotype=geno[,-c(1:14)]
#PCs=prcomp(genotype)
#str(PCs$x[,1:2])

loessW=VScan(x=genotype,y=(traitq),methods ="loess")
#loessWcv=VScan(x=genotype,y=(traitq),U=PCs$x[,1:2],methods ="loess")

lmW=VScan(x=genotype,y=(traitq),methods ="lm")
#lmWcv=VScan(x=genotype,y=(traitq),U=PCs$x[,1:2],methods ="lm")


Loci<-rep("Neutral", 1000)
Loci[c(201,211,221,231,241,251,261,271,281,291)]<-"QT"
Selected_Loci<-Loci[-which(Loci=="Neutral")]

library(ggplot2)
## Manhattan plot
g1=ggplot() +
  geom_point(aes(x=which(Loci=="Neutral"), 
  y=-log10(loessW$p_norm[-which(Loci!="Neutral")])), col = "gray83") +
  geom_point(aes(x=which(Loci!="Neutral"),
  y=-log10(loessW$p_norm[-which(Loci=="Neutral")]), colour = Selected_Loci)) +
  xlab("SNPs") + ylab("-log10(p-value)") +ylim(c(0,35))+theme_bw()

g1


g2=ggplot() +
  geom_point(aes(x=which(Loci=="Neutral"), 
  y=-log10(lmW$p_norm[-which(Loci!="Neutral")])), col = "gray83") +
  geom_point(aes(x=which(Loci!="Neutral"), 
  y=-log10(lmW$p_norm[-which(Loci=="Neutral")]), colour = Selected_Loci)) +
  xlab("SNPs") + ylab("-log10(p-value)") +ylim(c(0,35))+theme_bw()

g2

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Association testing using local polynomial fitting}% use one of  RShowDoc("KEYWORDS")

