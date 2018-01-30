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
  save(distribution, "path/to/size/distribution")
 ```

2. You can visualize the size distribution using the function `plot.size.distribution`
 ```r
 load("path/to/size/distribution"))
 dist <- distribution[[1]] # choose 1 or 2 if you want to see the frequency or count for the size distribution, respectively
 mode <- "log" # choose "log" or "lin" if you want to plot in logarithmic or linear scale, respectively

 plot.size.distribution(dist, mode)
 ```

3. Run the model
 ```r
 # The model needs the Photosynthetic Active Radiation in order to estimate the growth rate. Here is an example:
 Par.path <- system.file("extdata", "Par.csv", package="ssPopModel")
 Par <- read.csv(Par.path)

 path.distribution <- system.file("extdata", "prochloro_distribution", package="ssPopModel")


 time.delay <- 0 # time.delay (in hour) set the start of the time series with respect to t0

 model <- run.ssPopModel(path.distribution, Par, time.delay) 
 ```

4. Merge data output
 ```r
 all.simulations <- list.files('path/to/model/output', full.names=T)
 output <- merge.model.output(all.simulations)
 ```
 
5. You can visualize the `gamma` and `delta` function using the function `plot.parameters`
 ```r
 plot.paramters(output)
 ```
 
6. You can visualize the projection of the size distribution using the function `plot.size.distribution`
 ```r
 dist <- output$Vproj # choose Vproj or Noj if you want to see the frequency or count for the size distribution, respectively
 mode <- "log" # choose "log" or "lin" if you want to plot in logarithmic or linear scale, respectively

 plot.size.distribution(dist, mode)
 ```

