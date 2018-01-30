\name{RandomCluster}
\alias{RandomCluster}
\alias{RANDOMCLUSTER}
\title{randomize dates among tips in the BEAST input file}
\description{
	This function is similar to "RandomDates" excepts that in "RandomCluster", samples are grouped into clusters
	and the shuffling procedure randomizes dates among the clusters but not within 
	(see manual for more details on this procedure). There are two distinct ways to group the samples into clusters.
	The first one is through the upload of a csv.file containing the names of the samples (as given in the XML) 
	and a cluster number. Any positive integer (i.e., positive number) can be used to identified cluster;
	if a "0" is given to any sample, it would be excluded from the procedure.
	The file containing the classification should be labeled: clusters."name".cvs.
	An example for such a file in the case of the Influenza dataset can be found distributed with the package. 

	In a second approach, a model-based clustering classification is automatically performed using the mclust library.
	In this case, the option loadCluster should be set to FALSE (loadCluster = F).
	If this option is chosen, the new classification is written in a csv file that is labeled: clusters."name".cvs.
	}
\usage{
RandomCluster(name, reps = 20, loadCluster = T, writeTrees = T)
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
    \item{loadCluster}{
        F or T (default T). If T, clusters are loaded from a cluster structure file.
		The file containing the cluster structure needs to follow the example provided.
		Any tip assigned to cluster "0" will not be included in any randomization.
		Tip dates will only be randomized between (and not within) clusters.
		The cluster file should be named "clusters.NAME.csv" where NAME is the XML file name.
		If F, clusters are calculated using the package "mclust" procedure and an output cvs file
		containing the cluster structure is produced.
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
	In each new file, the date values are randomized among tips.
}
\references{
	Rieux, A. and Khatchikian, C.E., 2017. TipDatingBeast: An R package to assist the implementation of 
	phylogenetic tip-dating tests using BEAST. Molecular ecology resources, 17(4), pp.608-613.
	Drummond AJ, Suchard MA, Xie D & Rambaut A (2012) Bayesian phylogenetics with BEAUti and the BEAST 1.7.
	Molecular Biology And Evolution 29: 1969-1973.
	Fraley C & Raftery AE (2002) Model-based clustering, discriminant analysis, and density
	estimation. Journal of the American Statistical Association 97: 611-631.
}
\examples{
\dontrun{
    # using the example files "Flu_BEAST_1.8.xml" and "clusters.Flu.csv" found in example folder
	RandomCluster("Flu_BEAST_1.8", reps = 20, loadCluster = T)
	# produce 20 replicate input files (.xml) in working directory
}
}
\keyword{BEAST Software}
\keyword{phylogenetics}

