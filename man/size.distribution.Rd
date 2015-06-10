\name{size.distribution}
\alias{size.distribution}
\title{Generate the size distribution of your population of interest based on SeaFlow data}
\description{The function constructed a time-binned size distribution for your population of interest using light scattering measurements from SeaFlow database. The measured size distribution is divided equally on a logarithmic scale of base 2}

\usage{
size.distribution(popcycle.location, popname, n.breaks=57, time.interval = 60)
}

\arguments{
\item{popcycle.location}{location of the SeaFlow database that contains the light scatter values of your population of interest}
\item{popname}{the name of your population of interest. The name must be the sane than in the SeaFlow database}
\item{param}{the name of the parameter that contains the light scattering values. Default value is `fsc_small`}
\item{n.breaks}{the number of size classes for the size distribution. Default value is 57}
\item{time.interval}{the time interval of the size distribution, in minutes. Default value is 60.}
}

\value{
A size distribution made of two dimensions, for frequency and counts of the size distribution, respectively.
}

\examples{
### Select the opp data
# popcycle.location <- system.file("extdata", package="ssPopModel")

### Below is an example to create hourly-binned size distribution of 57 size classes for Prochlorococcus  
# dist <- size.distribution(popcycle.location, popname = 'prochloro', param="fsc_small", n.breaks=57, time.interval = 60)

}