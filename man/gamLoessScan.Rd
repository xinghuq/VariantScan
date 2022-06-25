\name{gamLoessScan}
\alias{gamLoessScan}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Variants (Biomarkers, e.g., genomic (genetic loci), transcriptomic (gene expression), epigenomic (methylations), proteomic(protein), metabolomic (metabolites) variants) Scanning and Association Tests Using Generalized Additive Model with  local polynomial regression (LOESS).
%%  ~~function to do ... ~~
}
\description{Fitting a Generalized Additive Mixed Models (GAMM) with local polynomial regression in association testing.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gamLoessScan(genotype, traits, U, cv_method = "adaptive_cv", 
model_metric = "RMSE", n_hyperparameter_search = 10, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{genotype}{Varants/genotypes matrix coding in reference allele (0,1,2) or variant count
%%     ~~Describe \code{genotype} here~~
}
  \item{traits}{Traits
%%     ~~Describe \code{traits} here~~
}
  \item{U}{ Covariates or confounding factors
%%     ~~Describe \code{U} here~~
}
  \item{cv_method}{Cross-validation
%%     ~~Describe \code{cv_method} here~~
}
  \item{model_metric}{Model performance metrics, based on which the optimal model is determined.
%%     ~~Describe \code{model_metric} here~~
}
  \item{n_hyperparameter_search}{Number of hyperparameters for tuning
%%     ~~Describe \code{n_hyperparameter_search} here~~
}
  \item{\dots}{other arguments passing to generalized additive mixed models (gam)
%%     ~~Describe \code{\dots} here~~
}
 
}
\details{Fits the specified generalized additive mixed model (GAMM) with LOESS smoothness.
%%  ~~ If necessary, more details than the description above ~~
}
\value{The weights of variants as well as their p-values
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Wood S.N. (2006b) Generalized Additive Models: An Introduction with R. Chapman and Hall/CRC Press.

Wang, Y. (1998) Mixed effects smoothing spline analysis of variance. J.R. Statist. Soc. B 60, 159-174.

Lin, X and Zhang, D. (1999) Inference in generalized additive mixed models by using smoothing splines. JRSSB. 55(2):381-400.

%% ~put references to the literature/web site here ~
}

\examples{
# not run
\donttest{f <- system.file('extdata',package='VariantScan')
infile <- file.path(f, "sim1.csv")
geno=read.csv(infile)

traitq=geno[,14]
genotype=geno[,-c(1:14)]
PCs=prcomp(genotype)
test=gamLoessScan(genotype =genotype,traits =(traitq),U=PCs$x[,1:2],n_hyperparameter_search=5)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{gamLoessScan}% use one of  RShowDoc("KEYWORDS")