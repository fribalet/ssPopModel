################################
### CREATE size distribution ###
################################
size.distribution <- function(db, vct.dir,binwidth=c(0.05, 0.002), log=c(TRUE, FALSE),
                      quantile=c(2.5, 50,97.5),
                      popname=c(NULL, 'prochloro','synecho','picoeuk','croco'),
                      channel=c('diam_mid','Qc_mid'),
                      lim=c(min, max)){

  QUANT <- as.numeric(quantile)
  CHANNEL <- as.character(channel)

  require(popcycle)
  require(zoo)

  # Get list of files, with list of outliers
  vct.table.all <- subset(get.vct.table(db), quantile == QUANT)
  outliers <- get.outlier.table(db)
  vct.table <- merge(outliers, vct.table.all, by='file', all.y=T, all.x=F)

  # Remove outliers from list of files
  vct.table <- subset(vct.table, flag==0)

  # Select on files that contains the population of interest
  if(!is.null(popname)) vct.table <- vct.table[vct.table$pop == popname,]

  # set bining of PDF
    # range of mean value
    if(is.null(lim)) lim <- range(vct.table[,CHANNEL])
    # 1/10 the smallest mean cell value
    min <- lim[1] / 10
    # 10 x the largest mean cell value
    max <- 10 * lim[2]

  breaks <- round(seq(min, max, by=binwidth),5)
  if(log) breaks <- round(2^seq(log2(min), log2(max), by=binwidth),5)

  # Get Time from metadata
  sfl <- get.sfl.table(db)

  # Get volume analyzed per sample
    # instrument serial Number
    inst <- get.inst(db)
    # volume of virtual core for each sample
  opp <- subset(get.opp.table(db), quantile == QUANT)


  # Return list of unique files
  vct.list <- unique(vct.table$file)

  i <- 2
  distribution <- data.frame(round(zoo::rollmean(breaks, 2),5))
  names(distribution) <- channel


  ##################################
  ### create PDF for each sample ###
  ##################################
  for(file.name in vct.list){

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
      dat <- vct[vct$pop == popname, CHANNEL]
    }else{
      dat <- vct[vct$pop != 'beads', CHANNEL]
    }

    # Get particle count in each bin
    dist <- data.frame(table(cut(dat, breaks)))

    # Get particle concentration in each bin (10^-3 cells µL-1 or 10^3 cells L-1)
    PDF <- round(1000 * dist$Freq / volume, 3)

    # add new data to the table
    distribution <- data.frame(cbind(distribution, PDF),check.names = FALSE)
    colnames(distribution)[i] <- time

    i <- i + 1
    flush.console()
  }

  # remove rows outside actual size range (zeros), but include the two extremes
  dim <- range(which(apply(distribution[,-c(1)], 1, sum) != 0))
  distribution <- distribution[(dim[1]-1):(dim[2]+1),]

  return(distribution)

}

##############################
### PLOT size distribution ###
##############################
plot.size.distribution <- function(distribution, log=TRUE, smooth=0){

    require(plotly)
    require(oce)

    # smooth distribution
    param <- data.frame(oce::matrixSmooth(as.matrix(distribution[,-c(1)]), pass=smooth))
    param <- data.frame(cbind(distribution[,1], param))
    colnames(param) <- colnames(distribution)
    varname <- colnames(distribution)[1]

    abundance <- as.matrix(param)[,-c(1)]
    time <- as.POSIXct(colnames(param), format = "%FT%T", tz="GMT")[-1]
    names(param)[1] <- 'tbd'

    p <- plotly::plot_ly(data=param, x= ~ time, y = ~ tbd, z = ~ abundance) %>%
                    add_surface() %>%
                    layout(scene = list(xaxis = list(autorange = "reversed"), yaxis= list(title=paste(varname))))

    if(log) p <- p %>% plotly::layout(scene = list(yaxis = list(type = "log",title=paste(varname))))

    return(p)
}
