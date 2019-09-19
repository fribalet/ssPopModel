# cruise <- c <- "SCOPE_1"
# path <- "/Volumes/OPPdata/SeaFlow-OPP/latest/"
# opp.dir <- paste0(path,cruise, "/",cruise,"_opp")
# vct.dir <- paste0(path,cruise, "/", cruise, "_vct")
# db <- paste0(path,cruise, "/",cruise,".db")
#
# min <- 0.1
# breaks <- round(min*2^(((1:m)-1)*0.125),5)
# breaks <- round(2^seq(log2(min), log2(max), length.out=200), 5)
#
# distribution <- size.distribution(db, vct.dir, quantile=2.5, popname='synecho', param="diam_mid", breaks=NULL)









#' Create particle size distribution
#'
#' @param db SQLite3 database file path.
#' @param vct.dir VCT file directory..
#' @param quantile Filtering function.
#' @param popname Population name. If NULL, all particles (except beads) will be used.
#' @param param Parameter from VCT to create the size distribution.
#'   Can be either diameter (diam_lwr, diam_mid or diam_upr)
#'   or carbon quotas (Qc_lwr, Qc_mid or Qc_upr)
#' @param delta Delta is a constant that define te resolution of the size disitribution.
#'   For Matrix model application, Delta must be chosen so that 1/delta is an integer
#' @param m M is the numnber of size classes. Must be chossen so it covers the size range
#' @param breaks Breaks must be a sequence defining the breaks for the size distribution (overwrite delta and m).
#' @return Size distribution
#' @examples
#' \dontrun{
#' distribution <- size.distribution(db, vct.dir, quantile=2.5, popname="synecho", channel="Qc_mid", delta=0.125, m=60)
#' }
#' @export size.distribution
size.distribution <- function(db, vct.dir, quantile=c(2.5, 50,97.5),
                      popname='prochloro',
                      param ='Qc_mid',
                      delta = 0.125, m= 60,
                      breaks = NULL){

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

  # Select on files that contains the population of interest
  if(!is.null(popname)) vct.table <- vct.table[vct.table$pop == popname,]


  if(is.null(breaks) & any(PARAM == c("Qc_lwr","Qc_mid","Qc_upr"))){
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  # set PDF template
    if(!is.wholenumber(1/delta)) print("WARNING: For Matrix model application, Delta must be chosen so that 1/delta is an integer")
    min <- min(vct.table[,paste0(PARAM,"_1q")]) / 10 # 1/10 the value of 25% percentile
    breaks <- round(min*2^(((1:m)-1)*delta),5)
  }

  if(is.null(breaks) & any(PARAM == c("diam_lwr","diam_mid","diam_upr"))){
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  # set PDF template
    if(!is.wholenumber(1/delta)) print("WARNING: For Matrix model application, Delta must be chosen so that 1/delta is an integer")
    min <- 3/4*pi*(0.5*min(vct.table[,paste0(PARAM,"_1q")])^3) / 10 # 1/10 the value of 25% percentile converted into sphere volume
    breaks <- min*2^(((1:m)-1)*delta)
    breaks <- round(2 * (4*breaks/(pi*3))^(1/3),5)
  }

  if(max(breaks) < max(vct.table[,paste0(PARAM,"_3q")])) print("WARNING: PDF template doesn't cover the the range of values")



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

  # Select on files that contains the population of interest
  if(!is.null(popname)) vct.table <- vct.table[vct.table$pop == popname,]

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
  dist <- list()
  for(file.name in vct.list){

    #file.name <- vct.list[2]
    message(round(100*i/length(vct.list)), "% completed \r", appendLF=FALSE)

    #retrieve time
    time <- sfl[sfl$file == file.name, 'date']

    # retrieve volume of virtual correct
      # flow rate (mL min-1)
      fr <- flowrate(sfl[sfl$file == file.name, 'stream_pressure'], inst=inst)$flow_rate
      # (µL min-1)
      fr <- fr * 1000
      # acquisition time (min)
      acq.time <- sfl[sfl$file == file.name, 'file_duration']/60
      # opp/evt
      opp.evt.ratio <- opp[opp$file == file.name,'opp_evt_ratio']
    volume <- opp.evt.ratio * fr * acq.time

    #retrieve data
    vct <- get.vct.by.file(vct.dir, file.name, quantile=QUANT)

    if (!is.null(popname)) {
      dat <- vct[vct$pop == popname, PARAM]
    }else{
      dat <- vct[vct$pop != 'beads', PARAM]
    }

    # Get particle count in each bin
    d <- table(cut(dat, breaks))

    # Get particle concentration in each bin (10^-3 cells µL-1 or 10^3 cells L-1)
    PSD <- round(1000 * d / volume, 3)

    PSD <- tibble(t(PSD))
    PSD <- add_column(PSD, .before=1, time=time)

    dist[[i]] <- PSD

    i <- i + 1
    flush.console()
  }

  distribution <- data.frame(matrix(unlist(dist), nrow=length(dist), byrow=T))
  colnames(distribution) <- c("time",names(d))

  return(distribution)

}




#####################################
### Bin size distribution by time ###
#####################################

bin.size.distribution.by.time(distribution, time.resolution="1 hour")

distribution$time <- as.POSIXct(distribution$time, format = "%FT%T", tz="GMT")
distribution %>%
        group_by(time = cut(time, breaks=time.resolution)) %>%
        summarise_all(LAT = mean(LAT, na.rm=T), LON = mean(LON, na.rm=T))











##############################
### PLOT size distribution ###
##############################
plot.size.distribution <- function(distribution, log=TRUE, smooth=0){

    require(plotly)
    require(oce)

    # smooth distribution
    param <- data.frame(oce::matrixSmooth(as.matrix(distribution[,-c(1)]), pass=smooth))
    colnames(param) <- colnames(distribution)[-1]

    abundance <- t(as.matrix(param))
    time <- as.POSIXct(distribution$time, format = "%FT%T", tz="GMT")
    size <- as.factor(colnames(distribution)[-1])

    p <- plotly::plot_ly(data=param, x= ~ time, y = ~ size, z = ~ abundance) %>%
                    add_surface() %>%
                    layout(scene = list(xaxis = list(autorange = "reversed")))

    if(log) p <- p %>% plotly::layout(scene = list(yaxis = list(type = "log",title=paste(varname))))

    return(p)
}
