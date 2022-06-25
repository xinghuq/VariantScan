\name{genmat}
\alias{genmat}
\alias{genmat.bed}
\alias{genmat.vcf}
\alias{genmat.gds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get genotype matrix from genomic data
%%  ~~function to do ... ~~
}
\description{Get genotype matrix from genomic data
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
genmat(genfile, sample.id = NULL, snp.id = NULL, 
snpfirstdim = NA, .snpread = NA, with.id = FALSE, verbose = TRUE, ...)


## S3 method for class 'genmat.bed'
genmat.bed(genfile, sample.id = NULL, snp.id = NULL, 
snpfirstdim = NA, .snpread = NA, with.id = FALSE, verbose = TRUE, ...)



## S3 method for class 'genmat.vcf'
genmat.vcf(genfile, sample.id = NULL, snp.id = NULL, 
snpfirstdim = NA, .snpread = NA, with.id = FALSE, verbose = TRUE, ...)



## S3 method for class 'genmat.gds'
genmat.gds(genfile, sample.id = NULL, snp.id = NULL, 
snpfirstdim = NA, .snpread = NA, with.id = FALSE, verbose = TRUE, ...)


}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genfile}{Genetic datasets containg sample ID and SNP ID, format includes bed (plink), vcf, or GDS file.
%%     ~~Describe \code{genfile} here~~
}
  \item{sample.id}{Sample ID
%%     ~~Describe \code{sample.id} here~~
}
  \item{snp.id}{SNP ID
%%     ~~Describe \code{snp.id} here~~
}
  \item{snpfirstdim}{whether genotypes are stored in the individual-major mode (TRUE), (i.e, list all SNPs for the first individual, and then list all SNPs for the second individual, etc) or (FALSE) for snp-major mode; if NA, determine automatically
%%     ~~Describe \code{snpfirstdim} here~~
}
  \item{.snpread}{internal use
%%     ~~Describe \code{.snpread} here~~
}
  \item{with.id}{whether return "sample.id" and "snp.id".
%%     ~~Describe \code{with.id} here~~
}
  \item{verbose}{whether printing information
%%     ~~Describe \code{verbose} here~~
}
  \item{\dots}{more arguments
%%     ~~Describe \code{\dots} here~~
}
}
\details{Effectively get genoptype matrix from various genotype formats, including bed, vcf, or gds.
%%  ~~ If necessary, more details than the description above ~~
}
\value{The function returns an integer matrix with values 0, 1, 2 or NA representing the number of reference allele when with.id=FALSE; or list(genotype, sample.id, snp.id) when with.id=TRUE. The orders of sample and SNP IDs in the genotype matrix are actually consistent with sample.id and snp.id in the GDS file, which may not be as the same as the arguments sampel.id and snp.id specified by users.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Zheng, X., & Weir, B. S. (2016). Eigenanalysis of SNP data with an identity by descent interpretation. Theoretical population biology, 107, 65-76.

%% ~put references to the literature/web site here ~
}

\examples{

\donttest{inp=SNPRelate::snpgdsExampleFileName()

genomat1=genmat.gds(inp)}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Genotype matrix}% use one of  RShowDoc("KEYWORDS")

