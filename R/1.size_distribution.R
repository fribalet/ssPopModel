################################
### CREATE size distribution ###
################################
size.distribution <- function(db, vct.dir,binwidth=0.0125, quantile=c(2.5, 50,97.5),
                              popname=c('prochloro','synecho','picoeuk',NULL),
                              channel=c('diam_lwr', 'diam_mid', 'diam_upr')){

  QUANT <- as.numeric(quantile)

  require(popcycle)
  require(zoo)
  # create template of size distribution
  min.volume <- 4/3 * pi * (0.4/2)^3 # smallest Pro
  max.volume <- 4/3 * pi * (6/2)^3 # largest Picoeuks
  breaks <- round(2^seq(log2(min.volume), log2(max.volume), by=binwidth),5)

  # Get Time from metadata
  sfl <- get.sfl.table(db)

  # Get info to determine volume analyzed per sample
  inst <- get.inst(db) # instrument serial Number
  opp <- subset(get.opp.table(db), quantile == QUANT) # volume of virtual core for each sample

  # Get list of files, with list of outliers
  vct.table.all <- subset(get.vct.table(db), quantile == QUANT)
  outliers <- get.outlier.table(db)
  vct.table <- merge(outliers, vct.table.all, by='file', all.y=T, all.x=F)

  # Remove outliers from list of files
  vct.table <- subset(vct.table, flag==0)

  # Select on files that contains the population of interest
  if(!is.null(popname)) vct.table <- vct.table[vct.table$pop == popname,]

  # Return list of unique files
  vct.list <- unique(vct.table$file)

  i <- 1
  distribution <- NULL
  for(file.name in vct.list){
    message(round(100*i/length(vct.list)), "% completed \r", appendLF=FALSE)

    #retrieve time
    time <- as.POSIXct(sfl[sfl$file == file.name, 'date'], format = "%FT%T", tz = "GMT")

    #retrieve volume of virtual correct
    fr <- 1000*flowrate(sfl[sfl$file == file.name, 'stream_pressure'], inst=inst)$flow_rate # flow rate (µL min-1)
    acq.time <- sfl[sfl$file == file.name, 'file_duration']/60 # acquisition time (min)
    opp.evt.ratio <- opp[opp$file == file.name,'opp_evt_ratio'] # opp/evt
    vol <- opp.evt.ratio * fr * acq.time

    #retrieve data
    vct <- get.vct.by.file(vct.dir, file.name, quantile=QUANT)

    if (!is.null(popname)) {
      size <- vct[vct$pop == popname, channel]
    }else{
      size <- vct[vct$pop != 'beads', channel]
    }

    # calcualte Volume from diameter, assuming spherical shape
    volume <- 4/3 * pi * (size/2)^3

    # Get particle count in each bin
    dist <- data.frame(table(cut(size, breaks)))

    # Get particle concentration in each bin (cells µL-1)
    PDF <- round(dist$Freq / vol, 5)

    # add new data to the table
    distribution <- data.frame(cbind(distribution, PDF),check.names = FALSE)
    colnames(distribution)[i] <- as.character(time)
    rownames(distribution) <- round(rollmean(breaks, 2),5)

    i <- i + 1
    flush.console()
  }

  return(distribution)

}


##############################
### PLOT size distribution ###
##############################
plot.size.distribution <- function(distribution, mode = c('log', 'lin'), ...){

    require(rgl)
    jet.colors <- viridis::viridis
    mode <- as.character(mode[1])

    param <- distribution
    percentile <- cut(unlist(param), 100)

    # in linear scale
    if(mode =='lin'){
        plot3d(rep(as.numeric(row.names(param)), dim(param)[2]),
                rep(as.POSIXct(colnames(param), tz="GMT"), each=dim(param)[1]) ,
                unlist(param),
                col=jet.colors(100)[percentile])#, xlab="size class", ylab="time", zlab="Frequency",...)
     }


    # in log scale
    if(mode =='log'){
        plot3d(log2(rep(as.numeric(row.names(param)), dim(param)[2])),
                rep(as.POSIXct(colnames(param)), each=dim(param)[1]) ,
                unlist(param),
                col=jet.colors(100)[percentile],  xlab="size class", ylab="time", zlab="Frequency",...)
    }
}
