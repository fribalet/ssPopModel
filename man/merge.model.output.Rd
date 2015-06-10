\name{merge.model.output}
\alias{merge.model.output}
\title{Merge the different simulations of the model}
\description{The function merge all data output and calculate the mean and standard deviation of the division rates, model parameters and projected size distribution.}

\usage{
merge.model.output(output.files, plot.raw=TRUE)
}

\arguments{
\item{output.files}{a list files that contains the model output}
\item{plot.raw}{the function plots either each model simulation individually (plot.raw=TRUE), or the calculated average (+ sd) of the simulations (plot.raw=FALSE). Default value is set to TRUE}
}

\value{
The function return a list of 5 objects:
1- Estimates. It contains the average (+ sd) of the estimates of division rates and model parameters
2- Vproj. It contains the mean frequencies of the projected size distribution
3- Nproj. It contains the mean count of the projected size distribution
4- Vproj.sd. It contains the  standard deviation for the frequencies of the projected size distribution
5- Vproj.sd. It contains the standard deviation for the count of the projected size distribution
}

\examples{
# output.files <- system.file("extdata", package="ssPopModel")
# binned <- merge.model.output(output.files, plot.raw=TRUE)
}