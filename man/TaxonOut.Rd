\name{TaxonOut}
\alias{TaxonOut}
\alias{TAXONOUT}
\title{single Taxon out input file for leave-one-out analysis in BEAST}
\description{
    This function produces a single input file for leave-one-out analyses using BEAST version 1 software.
	In this analysis, the date of the choosen taxon is estimated using the remaining taxa dates.
	The function "ListTaxa" is intended to help identifying the order (parameter "takeOut")
	of the taxon desired for the leave-out-analysis. 
}
\usage{
TaxonOut(name, takeOut)
}
\arguments{
    \item{name}{
        Name of the input file should be a .xml file generated using BEAUTi. Quote the name ("example").
		Do not included .xml.
    }
	\item{takeOut}{
        Taxon order for the take-one-out analysis.
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
	Rieux A & Khatchikian, C. Unpublished.
	Drummond AJ, Suchard MA, Xie D & Rambaut A (2012) Bayesian phylogenetics with BEAUti and the BEAST 1.7
	Molecular Biology And Evolution 29: 1969-1973.
}
\examples{
\dontrun{
    # using the example files "Flu.xml" found in example folder.
	# example using the 5th taxon ("CHICKEN_HONGKONG_915_1997")
	TaxonOut("Flu", 5)
	# produce a single input files without the corresponding taxon ("CHICKEN_HONGKONG_915_1997")
}
}
\keyword{BEAST Software}
\keyword{phylogenetics}