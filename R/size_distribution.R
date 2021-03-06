#' Create particle size distribution
#'
#' @param db SQLite3 database file path.
#' @param vct.dir VCT file directory..
#' @param quantile Filtering function.
#' @param param Parameter from VCT to create the size distribution.
#'   Can be either diameter (diam_lwr, diam_mid or diam_upr)
#'   or carbon quotas (Qc_lwr, Qc_mid or Qc_upr)
#' @param breaks Breaks must be a vector of values defining the breaks for the size distribution.
#' @return Size distribution
#' @examples
#' \dontrun{
#'  
#' breaks <- 'something here'
#' distribution <- create_PSD(db, vct.dir, quantile=50, channel="Qc_mid", breaks)
#' }
#' @export 
create_PSD <- function(db, vct.dir, quantile=50, param = 'Qc_mid', breaks){

  QUANT <- as.numeric(quantile)
  PARAM <- as.character(param)

  require(popcycle)
  require(tidyverse)
  require(zoo)

  # Get list of files, with list of outliers
  vct.table.all <- subset(get.vct.table(db), quantile == QUANT)
  outliers <- get.outlier.table(db)
  vct.table <- merge(outliers, vct.table.all, by='file', all.y=T, all.x=F)

  # Remove outliers from list of files
  vct.table <- subset(vct.table, flag==0)

  # Get Time from metadata
  sfl <- get.sfl.table(db)

  # Get volume analyzed per sample
    # instrument serial Number
    inst <- get.inst(db)
    # volume of virtual core for each sample
  opp <- subset(get.opp.table(db), quantile == QUANT)


  # Return list of unique files
  vct.list <- unique(vct.table$file)

  # Get list of files, with list of outliers
  vct.table.all <- subset(get.vct.table(db), quantile == QUANT)
  outliers <- get.outlier.table(db)
  vct.table <- merge(outliers, vct.table.all, by='file', all.y=T, all.x=F)

  # Remove outliers from list of files
  vct.table <- subset(vct.table, flag==0)

  # Get Time from metadata
  sfl <- get.sfl.table(db)

  # Get volume analyzed per sample
    # instrument serial Number
    inst <- get.inst(db)
    # volume of virtual core for each sample
  opp <- subset(get.opp.table(db), quantile == QUANT)


  # Return list of unique files
  vct.list <- unique(vct.table$file)




  ### create PSD for each timepoint 
  i <- 1
  distribution <- tibble()
  for(file.name in vct.list){

    #file.name <- vct.list[2]
    message(round(100*i/length(vct.list)), "% completed \r", appendLF=FALSE)

    #retrieve time
    time <- sfl[sfl$file == file.name, 'date']

    # retrieve volume of stream
    # flow rate (mL min-1)
    fr <- popcycle::flowrate(sfl[sfl$file == file.name, 'stream_pressure'], inst=inst)$flow_rate
    # convert to microL min-1
    fr <- fr * 1000
    # acquisition time (min)
    acq.time <- sfl[sfl$file == file.name, 'file_duration']/60
    # retrieve opp/evt
    opp.evt <- opp[opp$file == file.name,'opp_evt_ratio']
    # volume in microL
    volume <- round(fr * acq.time * opp.evt, 6)

 
    # retrieve PSD data
    vct <- get.vct.by.file(vct.dir, file.name, quantile=QUANT)

    # remove Beads from data
    dat <- vct[vct$pop != 'beads', c(PARAM,'pop')]

    # Get particle concentration in each bin for each population
    PSD <-NULL
    phyto <- unique(dat$pop)
    for(p in phyto){
      # get particle count
      psd <- table(cut(dat[which(dat$pop == p),PARAM], breaks))
      PSD <- rbind(PSD, psd)
    }
  
    # add time and population name
    PSD <- as_tibble(PSD)
    PSD <- PSD %>% add_column(time=as.POSIXct(time,format = "%FT%T", tz = "GMT"), pop=as.character(phyto), volume, .before=1)

    # bind data together
    distribution <- bind_rows(distribution, PSD)

    i <- i + 1
    flush.console()
  }
 
return(distribution)

}


#' Manipulate the size distribution created by FCSplankton::create_PSD(). 
#' Calculate the sum of particles in each size class over specific temporal resolution; transform the header
#'
#' @param distribution Particle size disitribution created by FCSplankton::create_PSD().
#'  i.e., a tibble of size distribution over time. First column must be time (POSIXt class object);
#'  Second column must name of the population; other columns represent the different size classes. 
#'  Size classes can represent either diameter or carbon quota (assuming spherical particles).
#' @param time.step Time step over which to sum the number of particles in each size class. Default 1 hour, must be higher than 3 minutes
#' @param diam.to.Qc Convert diameter to carbon quotas as described in
#'  Menden-Deuer, S. and Lessard, E. J. Carbon to volume relationships for dinoflagellates, diatoms, and other protist plankton.
#'  Limnol. Oceanogr. 45, 569–579 (2000).
#' @param Qc.to.diam Convert carbon quota into diameter (reciprocal of diam.to.Qc)
#' @param interval.to.geomean Transform size class intervals to geometric mean values 
#' (i.e. convert breaks (min, max] to geometric mean defined as sqrt(mean*max). 
#' @return Size distribution 
#' @name transform_PSD
#' @examples
#' \dontrun{
#' distribution <- transform_PSD(distribution, time.step="1 hour")
#' }
#' @export
transform_PSD <- function(distribution, time.step="1 hour", 
                                        diam.to.Qc=FALSE, 
                                        Qc.to.diam=FALSE, 
                                        interval.to.geomean=FALSE){
  
  # Check that 'time' is a POSIXt class object 
  if(! lubridate::is.POSIXt(distribution$time)){
  print("Time is not recognized as POSIXt class")
  stop
  }

  # Check that 'pop' column is there 
  if(!any(names(distribution)=='pop')){
    print("column 'pop' is missing")
  stop
  }

  # Menden-Deuer, S. & Lessard conversion factors
  d <- 0.261; e <- 0.860
  # convert size interval (factors) into data.frame
  breaks <- strsplit(sub("\\]","",sub("\\(","",colnames(distribution)[-c(1:3)])),",")

  if(Qc.to.diam){
    #convert Qc into diam using the Menden-Deuer conversion
    b <- lapply(breaks, function(x) round(2*(3/(4*pi)*(as.numeric(x)/d)^(1/e))^(1/3),6))
    colnames(distribution)[-c(1:3)] <- sub("\\)","\\]", sub("c","",as.character(b)))
  }

  if(diam.to.Qc){
    # convert diam into Qc using the Menden-Deuer conversion
    b <- lapply(breaks, function(x) round(d*(4/3*pi*(0.5*as.numeric(x))^3)^e,6))
    colnames(distribution)[-c(1:3)] <- sub("\\)","\\]", sub("c","",as.character(b)))
  }

  if(interval.to.geomean){
    # transform size class intervals to mean values (i.e. convert breaks (min, max] to geom mean). 
    midval <- unlist(list(lapply(breaks, function(x) sqrt(mean(as.numeric(x))*max(as.numeric(x))))))
    colnames(distribution)[-c(1:3)] <- round(midval,6)
  }

  # Calculate the mean in each size class over new time interval
  dist.bin <- distribution %>%
        group_by(time = cut(time, breaks=time.step), pop) %>%
        summarise_all(list(sum))

  # time was converted to factor, and need to be convereted back to POSIXt
  dist.bin$time <- as.POSIXct(dist.bin$time, tz='GMT')

  return(dist.bin)

}


#' Plot size distribution
#'
#' @param distribution Data frame of size distribution over time, (x time; y size classes). 
#' First column must be time (POSIXt class object), second column must name of the population;
#' other columns represent the different size classes.
#' Size classes can represent either diameter or carbon quota (assuming spherical particles).
#' @param lwd Line width for the lines
#' @param z.type "lin" for linear scaling of z values, "log" for logarithmic scaling
#' @return Plot carbon biomass in each size class
#' @name plot_PSD
#' @examples
#' \dontrun{
#' plot_PSD(distribution)
#' }
#' @export
plot_PSD <- function(distribution, lwd=4, z.type='log'){

    require(plotly)
    group.colors <- c(unknown="grey", 
                      prochloro=viridis::viridis(4)[1],
                      synecho=viridis::viridis(4)[2],
                      picoeuk=viridis::viridis(4)[3], 
                      croco=viridis::viridis(4)[4])

    # remove volume
    distribution <- distribution[,-c(3)]

    # convert time as factor to be compatible with plotting
    distribution$time <- as.factor(distribution$time)

    # format dat to be compatible with scatter3d
    d <- reshape2::melt(distribution)

    # order data by time
    d <- d[order(d$time),]

    plotly::plot_ly() %>%
          plotly::add_trace(data=d, x= ~ time, y = ~ variable, z = ~ value, 
                            type='scatter3d', mode='lines', line=list(width=lwd), 
                            color=~pop, colors=group.colors) %>%
          plotly::layout(scene = list(xaxis = list(autorange = "reversed"),
                              yaxis = list(title="size classes"),
                              zaxis = list(title="", type= z.type)))
                              #"Carbon (mg L<sup>-1</sup>"
                              #"Abundance (cells µL<sup>-1</sup>"
}
