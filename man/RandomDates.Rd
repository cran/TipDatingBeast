\name{RandomDates}
\alias{RandomDates}
\alias{RANDOMDATES}
\title{randomize dates among tips in the BEAST input file}
\description{
    This function read an XML file (BEAST version 1), randomize dates among all tips and produce
	a new XML input file with such randomization.
	The process is repeated up to the number of replicates (default is 20). 
    }
\usage{
RandomDates(name, reps = 20, writeTrees = T)
}
\arguments{
    \item{name}{
        The name of the original XML-formatted input file on which to apply the date-randomization procedure.
		Quote the name ("example"). The .xml extension should not be included. 
    }
    \item{reps}{
        The number of repetions required by the user. There will be as many date-randomized datasets produced
		as the value of reps (default = 20).
	}
	\item{writeTrees}{
		This argument has no function in the current version (default = T).
		}
}
\details{
    The function works only with a .xml file generated with BEAUti
}
\value{
    The function returns one or many files (the number is set by the "reps" argument; default is 20)
	In each new file, the date values are randomized among all tips.
}
\references{
	Rieux, A. and Khatchikian, C.E., 2017. TipDatingBeast: An R package to assist the implementation of 
	phylogenetic tip-dating tests using BEAST. Molecular ecology resources, 17(4), pp.608-613.
	Drummond AJ, Suchard MA, Xie D & Rambaut A (2012) Bayesian phylogenetics with BEAUti and the BEAST 1.7
	Molecular Biology And Evolution 29: 1969-1973.
}
\examples{
\dontrun{
    # using the example files "Flu-BEAST-1.8.xml" found in example folder.
	# example file can be found in the example folder
	RandomDates("Flu-BEAST-1.8", reps = 20)
	# produce 20 replicate input files (.xml) in working directory
}
}
\keyword{BEAST Software}
\keyword{phylogenetics}

