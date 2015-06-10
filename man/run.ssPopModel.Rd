\name{run.ssPopModel}
\alias{run.ssPopModel}
\title{Size-structured matrix population model}
\description{The function generate estimates cell division rates based on the dynamics of the size distribution over a 24-hr period using a size-structured matrix population model. The model  uses 4 parameters, b, dmax, E* and gmax. To estimate these parameters, the model is fitted to observations of the number of cells in each size class during a 24-hr period. The projected size distribution at time t+1 is calculated from the observed distribution at time t by matrix multiplication (using dt = 10), where the fraction of cells in each size class is based on the formulas for cell growth, cell division and stasis. The parameters are chosen by minimizing the weighted sum of the squared deviations between the observed and projected size distributions, assuming that the parameters remain constant within a day but may change from day to day.  The result of the optimization is a best-fitting projection matrix and an estimated division rate for each time point of the 24-hr period. 

For more details on the matrix population model, see Ribalet et al. Light-driven synchrony of Prochlorococcus cell growth and mortality in the subtropical Pacific gyre. Proceedings of the National Academy of Sciences USA, in press
}

\usage{
run.ssPopModel(distribution, Par, time.delay=0,  dt=10)
}

\arguments{
\item{distribution}{size distribution object generated from the size.distribution function. }
\item{Par}{The Photosynthetic Active Radiation data. The dataframe must include a column named time that contains data with the following time format: Y-m-d H:M:S}
\item{time.delay}{the time delay (in hour) of the 24-hr simulation with respect to the first data point of the size distribution. Value can range from 0 to 24.}
\item{dt}{The projection matrix gives the fraction of cells in size class j at time t that become cells in size class i at time t+dt. Since the model allows a cell to undergo only one transition from one size to another, dt must be low. Default value is set dt =10, similar to the original model}
}

\value{
A projection matrix, estimates of division rates and the values of the 4 parameters for each 24-hr period of the time series.
}

\examples{
# Par <- read.csv('/Volumes/seaflow/SCOPE_1/Par_SCOPE_1.csv')

# distribution <- load('/Volumes/seaflow/SCOPE_1/prochloro_dist_Ncat57_SCOPE_1')

# time.delay <- 0 # time.delay (in hour) set the start of the time series with respect to t0

# model <- run.ssPopModel(distribution, Par, time.delay) 
}