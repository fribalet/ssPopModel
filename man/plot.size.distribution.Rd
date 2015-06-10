\name{plot.size.distribution}
\alias{plot.size.distribution}
\title{Plot the size distribution of your population of interest in 3D}
\description{The function display in 3D the time-binned size distribution for your population of interest created by the size.distribution function}

\usage{
plot.size.distribution(distribution, mode="log")

}
\arguments{
\item{distribution}{A size distribution, either the frequency or counts of the size distribution}
\\item{mode}{Whether to display data in logarithmic ("log") or linear scale ("lin"). Default value is set to "log".}
}
\value{
An interactive 3D plot}

\examples{

### Select the opp data
# distribution <- system.file("extdata", package="ssPopModel")
# Vhist <- distribution[[1]]
# plot.size.distribution(Vhist)
}