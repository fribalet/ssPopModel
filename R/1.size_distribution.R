size.distribution <- function(db, opp.dir, vct.dir, popname, volume.width=0.07, time.interval = 60){

    require(popcycle)


      ############################################
      ### Get the range of 'param' for 'phyto' ###
      ############################################

        # Get the time range
        stat <- get.stat.table(db, flag=TRUE)
        try(stat <- subset(stat, flag ==0), silent=T)
        stat$time <- as.POSIXct(stat$time,format="%FT%T",tz='GMT')

          if(is.null(popname)){phyto.stat <- stat
              }else{phyto.stat <- subset(stat, pop == popname)}

        time.range <- range(phyto.stat$time)
        time <- seq(time.range[1],time.range[2] , by=time.interval*60) # cut the time series according to time interval


        # Get the beads data
        # print(paste("obtaining the median ", param, "of beads for normalization"))
         m.beads <- median(subset(stat, pop =='beads' & time > time.range[1] & time < time.range[2])[,"fsc_small_mean"])

        # Get Volbins
             # Get the volume range for 'phyto'
             #print(paste("obtaining the range in", param, "for", popname))
             param.phyto <- get.vct.stats.by.date(db, time.range[1], time.range[2])

               if(is.null(popname)){param.phyto <- subset(param.phyto, pop!='beads')
                   }else{param.phyto <- subset(param.phyto, pop==popname)}

             param.range <- c(quantile(param.phyto[,"fsc_small_min"],0.01), quantile(param.phyto[,"fsc_small_max"],0.99))
             norm.param.range <- param.range / m.beads
             volume.range <- round(1.918*(norm.param.range^0.524),4)
             # if(inst == 740) biomass.range <- round(4.753*(norm.param.range^1.235),4)
             # if(inst == 751) biomass.range <- round(5.401*(norm.param.range^1.622),4)
             biomass.range <- volume.range * 0.220 # Booth 1988 (Burbage & Binder found Qc <- 0.5 * norm.param.range ^(1/1.74))

             volbins.cut <- 2^seq(log2(volume.range[1]), log2(volume.range[2]), by=volume.width)
             biobins.cut <- 2^seq(log2(biomass.range[1]), log2(biomass.range[2]), length.out=length(volbins.cut))




        ##################################
        ### Generate SIZE distribution ###
        ##################################
        print(paste("generating", popname, "size distribution binned in", length(volbins.cut), "size classes, in ",time.interval, "minutes time interval"))

         i <- 0
        Vhist <- Bhist <- Ntot  <- Time <- NULL
        for( t in time){
             message(round(100*i/length(time)), "% completed \r", appendLF=FALSE)

            tryCatch({
            #get the opp for phyto
            #t <- time[133]
            t <- as.POSIXct(t, origin="1970-01-01", tz='GMT')
            stat.subset <- subset(stat, pop== popname & flag==0 & time >=t & time < t+60*time.interval)
            opp <- try(get.opp.by.file(opp.dir, stat.subset$file, quantile=50, vct.dir=vct.dir, pop=popname, channel='fsc_small'))
            if(class(opp) == "try-error" | nrow(opp) < 10){
                next
                }

            # convert normalized FSC to Volume and Biomass
                norm.fsc <- opp[,"fsc_small"]/m.beads
                volume <- round(1.918*(norm.fsc^0.524),4)
                # if(inst == 740) biomass <- round(4.753*(norm.fsc^1.235),4)
                # if(inst == 751) biomass <- round(5.401*(norm.fsc^1.622),4)
                biomass <- volume * 0.220

            # create the frequency distribution of Volume and Biomass
              dens <- hist(volume, breaks=volbins.cut, plot=F)
              freq.dist <-  dens$density*diff(dens$breaks) # convert density to frequency
            Vhist <- data.frame(cbind(Vhist, freq.dist))

              dens2 <- hist(biomass, breaks=biobins.cut, plot=F)
              freq.dist2 <-  dens2$density*diff(dens2$breaks) # convert density to frequency
            Bhist <- data.frame(cbind(Bhist, freq.dist2))

                n <- mean(stat.subset$abundance)
            Ntot <- data.frame(cbind(Ntot, n))

            Time <- c(Time, t)

            } , error = function(e) {print(paste("Encountered error at ", t))})
            i <-  i + 1
            flush.console()
        }


        #######################
        ### SAVE the output ###
        #######################
        volbins <- dens$mids
        biobins <- dens2$mids

        colnames(Vhist) <- colnames(Bhist) <- colnames(Ntot) <- as.character(Time)
        rownames(Vhist) <- volbins
        rownames(Bhist) <- biobins
        rownames(Ntot) <- "Ntot"

        distribution <- list()
            distribution[[1]] <- Vhist
            distribution[[2]] <- Bhist
            distribution[[3]] <- Ntot

         print("done")
        return(distribution)

        }







##############################
### PLOT size distribution ###
##############################
plot.size.distribution <- function(freq.distribution, mode = c('log', 'lin'), ...){

    require(rgl)
    jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

    mode <- as.character(mode[1])

    param <- freq.distribution
    percentile <- cut(unlist(param), 100)

    # in linear scale
    if(mode =='lin'){
        plot3d(rep(as.numeric(row.names(param)), dim(param)[2]),
                rep(as.numeric(colnames(param)), each=dim(param)[1]) ,
                unlist(param),
                col=jet.colors(100)[percentile], xlab="size class", ylab="time", zlab="Frequency",...)
     }


    # in log scale
    if(mode =='log'){
        plot3d(log2(rep(as.numeric(row.names(param)), dim(param)[2])),
                rep(as.numeric(colnames(param)), each=dim(param)[1]) ,
                unlist(param),
                col=jet.colors(100)[percentile],  xlab="size class", ylab="time", zlab="Frequency",...)
    }
}
