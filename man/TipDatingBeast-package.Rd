\name{TipDatingBeast-package}
\alias{TipDatingBeast-package}
\alias{TipDatingBeast}
\docType{package}
\title{Tip dating of phylogenetic trees with BEAST}
\description{
    Assist performing tip-dating of phylogenetic trees with BEAST. Main functions include
	randomization of dates among tips and producing new input files with such randomization
	and generating input files for leave-one-out analyses using BEAST software.
	}
\details{
\tabular{ll}{
Package: \tab TipDatingBeast\cr
Type: \tab Package\cr
Version: \tab 1.6\cr
Date: \tab 2015-07-28\cr
License: \tab GPL (>= 2)\cr
}
The RandomDates function randomize tip dates (BEAST v1).
The RandomDates2 function randomize tip dates (BEAST v2).
The RandomCluster function ramdomize tip dates between group of tips (BEAST v1).
The RandomCluster2 function ramdomize tip dates between group of tips (BEAST v2).
The TaxaOut function generates input files for leave-one-analysis.
The TaxonOut function generates a single input file by removing the date of a particular taxon.
The ListTaxa function displays the names and order of the taxa in the xml file.
}
\author{
	Rieux A & Khatchikian, C.
Maintainer: C. Khatchikian <ckhatchikian@gmail.com>
}
\references{
	Rieux A & Khatchikian, C. Unpublished.
	Drummond AJ, Suchard MA, Xie D & Rambaut A (2012) Bayesian phylogenetics with BEAUti and the BEAST 1.7.
	Molecular Biology And Evolution 29: 1969-1973.
}
\keyword{BEAST Software}
\keyword{phylogenetics}
