#' Create particle size distribution
#'
#' @param db SQLite3 database file path.
#' @param vct.dir VCT file directory..
#' @param quantile Filtering function.
#' @param param Parameter from VCT to create the size distribution.
#'   Can be either diameter (diam_lwr, diam_mid or diam_upr)
#'   or carbon quotas (Qc_lwr, Qc_mid or Qc_upr)
#' @param breaks Breaks must be a sequence defining the breaks for the size distribution.
#' @return Size distribution
#' @examples
#' \dontrun{
#' distribution <- size.distribution(db, vct.dir, quantile=50, popname="synecho", channel="Qc_mid", delta=0.125, m=60)
#' }
#' @export size.distribution

size.distribution <- function(db, vct.dir, quantile=50,
                      param ='Qc_mid',
                      breaks){

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

  ##################################
  ### create PDF for each sample ###
  ##################################
  i <- 1
  distribution <- NULL
  for(file.name in vct.list){

    #file.name <- vct.list[2]
    message(round(100*i/length(vct.list)), "% completed \r", appendLF=FALSE)

    #retrieve time
    time <- sfl[sfl$file == file.name, 'date']

    # retrieve volume of virtual correct
      # flow rate (mL min-1)
      fr <- flowrate(sfl[sfl$file == file.name, 'stream_pressure'], inst=inst)$flow_rate
      # convert to µL min-1
      fr <- fr * 1000
      # acquisition time (min)
      acq.time <- sfl[sfl$file == file.name, 'file_duration']/60
      # opp/evt
      opp.evt.ratio <- opp[opp$file == file.name,'opp_evt_ratio']
    volume <- opp.evt.ratio * fr * acq.time

    # retrieve data
    vct <- get.vct.by.file(vct.dir, file.name, quantile=QUANT)

    # remove Beads from data
    dat <- vct[vct$pop != 'beads', c(PARAM,'pop')]

    # Get particle concentration in each bin for each population
    PSD <-NULL
    phyto <- unique(dat$pop)
    for(p in phyto){
      # get particle count
      psd <- table(cut(dat[which(dat$pop == p),PARAM], breaks))
      # Get particle concentration in each bin (10^-3 cells µL-1 or 10^3 cells L-1)
      psd <- round(1000 * psd / volume, 3)
      PSD <- rbind(PSD, psd)
    }

    # add time and population name
    PSD <- cbind(time, pop=as.character(phyto), PSD)

    # bind data together
    distribution <- data.frame(rbind(distribution, PSD),check.names=F)

    i <- i + 1
    flush.console()
  }
  #convert data frame to tibble, with correct classes (tibble wrongly assumed the class of each column, arghh!!!!)
  distribution <- as_tibble(distribution)
  distribution$time <- as.POSIXct(distribution$time,format = "%FT%T", tz = "GMT")
  distribution$pop <- as.character(distribution$pop)
  distribution[,-c(1,2)] <- mutate_all(distribution[,-c(1,2)], function(x) as.numeric(as.character(x)))

  return(distribution)

}



#' Bin size distribution by time
#'
#' @param distribution Size distribution created by size.distribution(). Time format must be compatible with POSIXt class
#'   Can be either diameter (diam_lwr, diam_mid or diam_upr)
#'   or carbon quotas (Qc_lwr, Qc_mid or Qc_upr)
#' @param time.step Time resolution (must be higher than 3 minutes). Default is 1 hour
#' @param diam.to.Qc Convert diameters into carbon quotas as described in
#' Menden-Deuer, S. & Lessard, E. J. Carbon to volume relationships for dinoflagellates, diatoms, and other protist plankton.
#' Limnol. Oceanogr. 45, 569–579 (2000).
#' @param Qc.to.diam Convert carbon quotas into diameters (reciprocal of Menden-Deuer, S. & Lessard, E. J. 2000)
#' @param abundance.to.biomass Calcualte total carbon biomass in each size class (abundance x Qc). Warning: If size class values represent diameters, make sure to set diam.to.Qc = TRUE.
#' @return Size distribution with temporal resolution defined by time.step
#' @examples
#' \dontrun{
#' distribution <- bin.distribution.by.time(distribution, time.step="1 hour")
#' }
#' @export bin.distribution.by.time

transform.size.distribution <- function(distribution, time.step="1 hour", diam.to.Qc=T, Qc.to.diam=F, abundance.to.biomass=F){

  if(! lubridate::is.POSIXt(distribution$time)){
  print("Time is not recognized as POSIXt class")
  stop
  }

  # Menden-Deuer, S. & Lessard conversion factors
  d <- 0.261; e <- 0.860
  # convert size interval (factors) into data.frame
  breaks <- strsplit(sub("\\]","",sub("\\(","",colnames(distribution)[-c(1,2)])),",")

  if(Qc.to.diam){
    #convert Qc into diam using the Menden-Deuer conversion
    b <- lapply(breaks, function(x) round(2*(3/(4*pi)*(as.numeric(x)/d)^(1/e))^(1/3),5))
    colnames(distribution)[-c(1,2)] <- sub("\\)","\\]", sub("c","",as.character(b)))
  }

  if(diam.to.Qc){
    # convert diam into Qc using the Menden-Deuer conversion
    b <- lapply(breaks, function(x) round(d*(4/3*pi*(0.5*as.numeric(x))^3)^e,5))
    colnames(distribution)[-c(1,2)] <- sub("\\)","\\]", sub("c","",as.character(b)))
    breaks <- strsplit(sub("\\]","",sub("\\(","",colnames(distribution)[-c(1,2)])),",")
  }

  if(abundance.to.biomass){
    # multiply abundance by carbon quotas to get carbon biomass in each size class
    midval <- unlist(list(lapply(breaks, function(x) mean(as.numeric(x)))))
    distribution[-c(1,2)] <- sweep(distribution[-c(1,2)], MARGIN=2, midval, `*`)
  }

  # Calculate the mean in each size class over new time interval
  dist.bin <- distribution %>%
        group_by(time = cut(time, breaks=time.step), pop) %>%
        summarise_all(list(mean))

  # time was converted to factor, and need to be convereted back to POSIXt
  dist.bin$time <- as.POSIXct(dist.bin$time, tz='GMT')

  return(dist.bin)

}



#' Plot size distribution
#'
#' @param distribution Size distribution created by size.distribution().
#'   Can be either diameter (diam_lwr, diam_mid or diam_upr)
#'   or carbon quotas (Qc_lwr, Qc_mid or Qc_upr)
#'   Must be based on carbon quotas (Qc_lwr, Qc_mid or Qc_upr) to be meaningful
#' @param lwd line width for the lines
#' @return Plot carbon biomass in each size class
#' @examples
#' \dontrun{
#' plot.size.distribution(distribution)
#' }
#' @export plot.biomass.distribution

plot.size.distribution <- function(distribution, lwd=4){

    require(plotly)
    group.colors <- c(unknown="grey", prochloro=viridis::viridis(4)[1],synecho=viridis::viridis(4)[2],picoeuk=viridis::viridis(4)[3], croco=viridis::viridis(4)[4])

    # convert time as factor to be compatible with plotting
    distribution$time <- as.factor(distribution$time)

    # format dat to be compatible with scatter3d
    d <- reshape2::melt(distribution)

    # order data by time
    d <- d[order(d$time),]

    plotly::plot_ly() %>%
          plotly::add_trace(data=d, x= ~ time, y = ~ variable, z = ~ value, type='scatter3d', mode='lines', line=list(width=lwd), color=~pop, colors=group.colors) %>%
          plotly::layout(scene = list(xaxis = list(autorange = "reversed"),
                              yaxis = list(title="size classes"),
                              zaxis = list(title="")))
                              #"Carbon (mg L<sup>-1</sup>"
                              #"Abundance (cells µL<sup>-1</sup>"

}
