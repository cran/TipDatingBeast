\name{RandomDates2}
\alias{RandomDates2}
\alias{RANDOMDATES2}
\title{randomize dates among tips in the BEAST input file}
\description{
    This function read an XML file (BEAST version 2), randomize dates among all tips and produce
	a new XML input file with such randomization.
	The process is repeated up to the number of replicates (default is 20). 
	}
\usage{
RandomDates2(name, reps = 20)
}
\arguments{
    \item{name}{
        Name of the input file should be a .xml file generated using BEAUTi. Quote the name ("example").
		Do not included .xml. 
    }
    \item{reps}{
        Number of replicates to be produced (default 20).
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
	Rieux A & Khatchikian, C. Unpublished.
	Drummond AJ, Suchard MA, Xie D & Rambaut A (2012) Bayesian phylogenetics with BEAUti and the BEAST 1.7
	Molecular Biology And Evolution 29: 1969-1973.
}
\examples{
\dontrun{
    # using the example files "Flu2.xml" found in example folder.
	# example file can be found in the example folder
	RandomDates("Flu2", 10)
	# produce 10 replicate input files (.xml) in working directory
}
}
\keyword{BEAST Software}
\keyword{phylogenetics}

