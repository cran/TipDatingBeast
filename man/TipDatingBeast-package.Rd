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
Version: \tab 1.0.8\cr
Date: \tab 2018-01-29\cr
License: \tab GPL (>= 2)\cr
}
The RandomDates function randomize tip dates.
The RandomCluster function ramdomize tip dates between group of tips.
The TaxaOut function generates input files for leave-one-analysis.
The TaxonOut function generates a single input file by removing the date of a particular taxon.
The ListTaxa function displays the names and order of the taxa in the xml file.
The PlotDRT function plots BEAST output of ramdomize tip dates analysis.
The PlotLOOCV function plots BEAST output of leave-one-analysis.
}
\author{
	Rieux A & Khatchikian, C.
Maintainer: C. Khatchikian <ckhatchikian@gmail.com>
}
\references{
	Rieux, A. and Khatchikian, C.E., 2017. TipDatingBeast: An R package to assist the implementation of 
	phylogenetic tip-dating tests using BEAST. Molecular ecology resources, 17(4), pp.608-613.
	Drummond AJ, Suchard MA, Xie D & Rambaut A (2012) Bayesian phylogenetics with BEAUti and the BEAST 1.7.
	Molecular Biology And Evolution 29: 1969-1973.
}
\keyword{BEAST Software}
\keyword{phylogenetics}
