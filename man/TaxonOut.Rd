\name{TaxonOut}
\alias{TaxonOut}
\alias{TAXONOUT}
\title{single Taxon out input file for leave-one-out analysis in BEAST}
\description{
    This function produces a single input file for leave-one-out analyses using BEAST version 1 software.
	In this analysis, the date of the chosen taxon is estimated using the remaining taxa dates.
	The function "ListTaxa" is intended to help identifying the order (parameter "takeOut")
	of the taxon desired for the leave-out-analysis. 
}
\usage{
TaxonOut(name, lBound = 0, hBound = 1.0E100, takeOut, writeTrees = T)
}
\arguments{
    \item{name}{
        Name of the input file should be a .xml file generated using BEAUTi. Quote the name ("example").
		Do not included .xml.
    }
	\item{lBound}{
		The uniform prior lower bound for the age of the missing taxa (default = 0)
	}
	\item{hBound}{
		The uniform prior higher bound for the age of the missing taxa (default = 1.0E100)
	}
	\item{takeOut}{
        Taxon order for the take-one-out analysis.
    }
	\item{writeTrees}{
		This argument has no function in the current version (default = T).
	}
}
\details{
    The function works only with a .xml file generated with BEAUti
}
\value{
    The function returns a single file to perform a leave-one-out analyses using BEAST software
	for the specific taxon. 
}
\references{
	Rieux, A. and Khatchikian, C.E., 2017. TipDatingBeast: An R package to assist the implementation of 
	phylogenetic tip-dating tests using BEAST. Molecular ecology resources, 17(4), pp.608-613.
	Drummond AJ, Suchard MA, Xie D & Rambaut A (2012) Bayesian phylogenetics with BEAUti and the BEAST 1.7
	Molecular Biology And Evolution 29: 1969-1973.
}
\examples{
\dontrun{
    # using the example files "Flu_BEAST_1.8.xml" found in example folder.
	# example using the 5th taxon ("CHICKEN_HONGKONG_915_1997")
	TaxonOut("Flu_BEAST_1.8", takeOut = 5)
	# produce a single input files without the corresponding taxon ("CHICKEN_HONGKONG_915_1997")
}
}
\keyword{BEAST Software}
\keyword{phylogenetics}