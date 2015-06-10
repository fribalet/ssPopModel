\name{plot.parameters}
\alias{plot.parameters}
\title{Plot the Gamma and Delta function based on the 4 model parameters}
\description{The function display the Delta function (i.e., how the probability of cell division increases as a power of cell volume) and the Gamma function (i.e., how the probability to grow depends on the incident radiations) based on the 4 model parameters.}

\usage{
plot.parameters(merged.model.output)

}
\arguments{
\item{merged.model.output}{An object from the merge.model.output function that contains the mean of the model parameters from the different simulations.}
}
\value{
Two panels that include the Delta and Gamma function.}

\examples{

### Select the opp data
# merged.data <- system.file("extdata", package="ssPopModel")
# plot.parameters(merged.data)
}