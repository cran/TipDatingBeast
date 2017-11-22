\name{ListTaxa}
\alias{ListTaxa}
\alias{LISTTAXA}
\title{List taxa names present in the BEAST input file}
\description{
    This function list the taxa names present in an input xml file for the BEAST version 1 software. 
	The function is intended to help using the function "TaxonOut" as it allows
    the identification of the order of each taxa present in the input file.
}
\usage{
ListTaxa(name)
}
\arguments{
    \item{name}{
        Name of the input file should be a .xml file generated using BEAUTi. Quote the name ("example").
		Do not included .xml. 
    }
}
\details{
    The function works only with a .xml file generated with BEAUti
}
\value{
    The function returns the names and order of the taxa present in an input xml file for the BEAST software.
}
\references{
	Rieux, A. and Khatchikian, C.E., 2017. TipDatingBeast: An R package to assist the implementation of 
	phylogenetic tip-dating tests using BEAST. Molecular ecology resources, 17(4), pp.608-613.
	Drummond AJ, Suchard MA, Xie D & Rambaut A (2012) Bayesian phylogenetics with BEAUti and the BEAST 1.7.
	Molecular Biology And Evolution 29: 1969-1973.
}
\examples{
\dontrun{
    # using the example files "Flu-BEAST-1.8.xml" found in example folder.
	# example file can be found in the example folder
	ListTaxa("Flu-BEAST-1.8")
	# list all 21 taxa in the file in the console 
}
}
\keyword{BEAST Software}
\keyword{phylogenetics}
