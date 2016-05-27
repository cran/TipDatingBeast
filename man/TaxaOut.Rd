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
TaxaOut(name, lBound = 0, hBound = 1.0E100, writeTrees = T)
}
\arguments{
    \item{name}{
        The name of the original XML-formatted input file on which to apply the LOOCV procedure
		(the .xml extension should be excluded). This xml should be set up so that earlier dates have
		lower numerical values (i.e., set direction = "forwards" when setting up date in BEAUti).
		Place the name between quotes ("example").
    }
    \item{lBound}{
		The uniform prior lower bound for the age of the missing taxa (default = 0)
	}
	\item{hBound}{
		The uniform prior higher bound for the age of the missing taxa (default = 1.0E100)
	}
	\item{writeTrees}{
		Set to False (F) if you do not want the trees to be written when running the date-randomized
		datasets in BEAST. To proceed the LOOCV, only the log files are needed (default = T).
		}
}
\details{
    The function works only with a .xml file generated with BEAUti.
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
    # using the example files "Flu_BEAST_1.8.xml" found in example folder.
	TaxaOut("Flu_BEAST_1.8")
	# produce 21 input files, each one without the corresponding taxon
}
}
\keyword{BEAST Software}
\keyword{phylogenetics}