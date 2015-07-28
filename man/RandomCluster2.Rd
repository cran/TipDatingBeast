\name{RandomCluster2}
\alias{RandomCluster2}
\alias{RANDOMCLUSTER2}
\title{randomize dates among tips in the BEAST input file}
\description{
    This function is an extension of the RandomDates function.
    It reads an XML file (BEAST version 2), randomize dates between clusters of tips
	(instead of among all tips as done by RandomDates) and produce a new input file with such randomization.
    Clusters of tips, in general based on dates, can either be provided by the user or automatically
	estimated using the package "mclust".
    The process is repeated up to the number of replicates (default is 20).
}
\usage{
RandomCluster2(name, reps = 20, loadCluster = T)
}
\arguments{
    \item{name}{
        Name of the input file should be a .xml file generated using BEAUTi. Quote the name ("example").
		Does not include .xml. 
    }
    \item{reps}{
        Number of replicates (and files) to be produced (default 20).
    }
    \item{loadCluster}{
        F or T (default T). If T, clusters are loaded from a cluster structure file.
		The file containing the cluster structure needs to follow the example provided.
		Any tip assigned to cluster "0" will not be included in any randomization.
		Tip dates will only be randomized between (and not within) clusters.
		The cluster file should be named "clusters.NAME.csv" where NAME is the XML file name.
		If F, clusters are calculated using the package "mclust" procedure and an output cvs file containing
		the cluster structure is produced.
		}
}
\details{
    The function works only with a .xml file generated with BEAUti
}
\value{
    The function returns one or many files (the number is set by the "reps" argument; default is 20)
	In each new file, the date values are randomized among tips.
}
\references{
	Rieux A & Khatchikian, C. Unpublished. \cr
	Drummond AJ, Suchard MA, Xie D & Rambaut A (2012) Bayesian phylogenetics with BEAUti and the BEAST 1.7.
	Molecular Biology And Evolution 29: 1969-1973.
	Fraley C & Raftery AE (2002) Model-based clustering, discriminant analysis, and density
	estimation. Journal of the American Statistical Association 97: 611-631.
}
\examples{
\dontrun{
    # using the example files "Flu2.xml" and "clusters.Flu2.csv" found in example folder.
	RandomDates2("Flu2", 10, T)
	# produce 10 replicate input files (.xml) in working directory
}
}
\keyword{BEAST Software}
\keyword{phylogenetics}

