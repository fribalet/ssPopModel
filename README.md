ssPopModel
==========

ssPopModel is a R package that uses size-structured matrix population model to estimate hourly division rates of microbial populations from SeaFlow data. These estimates are independent of variations in cell abundance and can be used to study physiologically-driven changes in population dynamics. The model is described in:

Ribalet et al. Light-driven synchrony of <i>Prochlorococcus</i> cell growth and mortality in the subtropical Pacific gyre. Proceedings of the National Academy of Sciences USA, in press

Here is a quick tutorial to get you started (a more detailed tutorial is in progress).

1. Create a time series of hourly-binned size distribution of <i>Prochlorococcus</i>
 ```r
  library(ssPopModel)
  popcycle.location <- "/Volumes/seaflow/SCOPE_1" # location of the SeaFlow database
  popname <- 'prochloro'
 
  distribution <- size.distribution(popcycle.location, popname, param="fsc_small", n.breaks=57, time.interval = 60)
 ```

 You can visualize the size distribution using the function `plot.size.distribution`
 ```r
 dist <- distribution[[1]] # choose 1 or 2 if you want to see the frequency or count for the size distribution, respectively
 mode <- "log" # choose "log" or "lin" if you want to plot in logarithmic or linear scale, respectively

 plot.size.distribution(dist, mode)
 ```

2. Run the model
 ```r
 # The model need the Photosynthetic Active Radiation in order to estimate the growth rate
 Par <- read.csv('/Volumes/seaflow/SCOPE_1/Par_SCOPE_1.csv')
 distribution <- load('/Volumes/seaflow/SCOPE_1/prochloro_dist_Ncat57_SCOPE_1')
 time.delay <- 0 # time.delay (in hour) set the start of the time series with respect to t0

 model <- run.ssPopModel(distribution, Par, time.delay) 
 ```

3. Merge data output
 ```r
 all.simulations <- list.files('path/to/model/output', full.names=T)
 output <- merge.model.output(all.simulations)
 ```
 
 You can visualize the `gamma` and `delta` function using the function `plot.parameters`
 ```r
 plot.paramters(output)
 ```
 
 You can visualize the projection of the size distribution using the function `plot.size.distribution`
 ```r
 dist <- output$Vproj # choose Vproj or Noj if you want to see the frequency or count for the size distribution, respectively
 mode <- "log" # choose "log" or "lin" if you want to plot in logarithmic or linear scale, respectively

 plot.size.distribution(dist, mode)
 ```

