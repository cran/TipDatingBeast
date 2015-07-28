\name{TaxaOut}
\alias{TaxaOut}
\alias{TAXAOUT}
\title{generates input files for leave-one-out analyses in BEAST}
\description{
    This function produces input files to perform leave-one-out analyses using BEAST version 1 software.
	As many files as taxa present in the input file are produced; each one leaving each
	consecutive taxon out for analysis.
}
\usage{
TaxaOut(name)
}
\arguments{
    \item{name}{
        Name of the input file should be a .xml file generated using BEAUTi.
		Place the name between quotes ("example").
		Does not include .xml.
    }
}
\details{
    The function works only with a .xml file generated with BEAUti version 1.
}
\value{
    The function returns as many files as taxa are present in the input file; each one leaving each
	consecutive taxon out for analysis.
}
\references{
	Rieux A & Khatchikian, C. Unpublished.
	Drummond AJ, Suchard MA, Xie D & Rambaut A (2012) Bayesian phylogenetics with BEAUti and the BEAST 1.7
	Molecular Biology And Evolution 29: 1969-1973.
}
\examples{
\dontrun{
    # using the example files "Flu.xml" found in example folder.
	TaxaOut("Flu")
	# produce 21 input files, each one without the corresponding taxon
}
}
\keyword{BEAST Software}
\keyword{phylogenetics}