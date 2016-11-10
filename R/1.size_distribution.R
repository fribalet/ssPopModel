# popcycle.location <- "/Volumes/seaflow/SCOPE_1"
# popname <- "prochloro"
# param <- "fsc_small"
# n.breaks <- 57
# time.interval <- 60 #minutes


############################################
### Get the range of 'param' for 'phyto' ###
############################################

size.distribution <- function(db, opp.dir, vct.dir, popname, param="fsc_small", n.breaks=57, time.interval = 60){

    require(popcycle)

        # Get the time range
        #print("getting the stat table from the database")
        stat <- get.stat.table(db, flag=TRUE)
        try(stat <- subset(stat, flag ==0), silent=T)
        stat$time <- as.POSIXct(stat$time,format="%FT%T",tz='GMT')

          if(is.null(popname)){phyto.stat <- stat
              }else{phyto.stat <- subset(stat, pop == popname)}

        time.range <- range(phyto.stat$time)
        time <- seq(time.range[1],time.range[2] , by=60*time.interval) # cut the time series according to time interval

        # Get the range of 'param' for 'phyto'
        #print(paste("obtaining the range in", param, "for", popname))
        param.phyto <- get.vct.stats.by.date(db, time.range[1], time.range[2])

          if(is.null(popname)){param.phyto <- subset(param.phyto, pop!='beads')
              }else{param.phyto <- subset(param.phyto, pop==popname)}

        param.range <- c(mean(param.phyto[,paste0(param, "_min")]), mean(param.phyto[,paste0(param, "_max")]))


        # Get the beads data
        # print(paste("obtaining the median ", param, "of beads for normalization"))
         m.beads <- median(subset(stat, pop =='beads' & time > time.range[1] & time < time.range[2])[,paste0(param, "_mean")])
        # Plot the light scattering of beads over time
        # plot.time(stat, popname='beads',param='fsc_small', ylim=c(1, 10^3.5))
        # abline(h=m.beads, col=2, lwd=3)

        # # BInned data according to 'time'interval'
        # beads <- subset(df, pop=='beads')
        # time.binned <- cut(beads$time, time, labels=F)
        # param.beads.binned <- as.vector(tapply(beads$fsc_small, time.binned, median))
        # time.beads.binned <- as.POSIXct(as.vector(tapply(beads$time, time.binned, mean)), origin="1970-01-01", tz='GMT')
        # points(time.beads.binned , param.beads.binned, type='o', col=2, pch=16)

        # # Smooth the data
        # spar <- 0.45 # smooothing parameter, the higher the more smoothing is applied.
        # smooth <- smooth.spline(time.beads.binned, param.beads.binned,spar=spar)
        # smooth.param.beads.binned <- spline(as.POSIXct(smooth$x,origin="1970-01-01",tz="GMT"), smooth$y, xout=as.POSIXct(smooth$x,origin="1970-01-01",tz="GMT"))
        # lines(smooth.param.beads.binned ,col=3,lwd=3) # visualize the smooth data and re-adjust the spar parameter if necessary.


        ###############################################################
        ### Generate SIZE distribution, binned into 'n.breaks' for each time interval ###
        ###############################################################
        print(paste("generating size distribution made of", n.breaks, "size class, binned by",time.interval, "minutes time interval"))

         i <- 0
        Vhist <- Ndist  <- Time <- NULL
        for( t in time){

            # t <- time[1]
             message(round(100*i/length(time)), "% completed \r", appendLF=FALSE)

            tryCatch({
            #get the opp for phyto
            t <- as.POSIXct(t, origin="1970-01-01", tz='GMT')
            opp.list <- subset(stat, pop== popname & flag==0 & time >=t & time < t+60*time.interval)$file

            pop <- try(get.opp.by.file(opp.dir, opp.list, vct.dir=vct.dir, pop=popname, channel=param))

            if(class(pop) == "try-error" | nrow(pop) < 10){
                next
                }

            # get Beads signal

            # get opp/evt ratio (used to calculate Ndist)
          ## opp.evt.ratio <- median(get.opp.evt.ratio.by.date(t, t+60*time.interval)$ratio)

            # create the size distribution of normalized forward scatter, using a Gaussian filter
            dens <- density(log2(pop[,param]/m.beads), n=n.breaks,from=log2(param.range[1]/m.beads) , to=log2(param.range[2]/m.beads), kernel='gaussian')
            freq.dist <-  dens$y*diff(dens$x)[1] # convert density to frequency
            freq.dist <- freq.dist/sum(freq.dist) # normailize the frequency to 1
                Vhist <- data.frame(cbind(Vhist, freq.dist))
             size.dist <- round(freq.dist * nrow(pop))
                Ndist <- data.frame(cbind(Ndist, size.dist))
                Time <- c(Time, t)

            } , error = function(e) {print(paste("Encountered error at ", t))})
            i <-  i + 1
            flush.console()
        }


        ######################################################################
        ### CONVERT normalized forward SCATTER by 1 micron beads to VOLUME ###
        ######################################################################
        print(paste("converting", param, "into volume"))
        norm.fsc <- 2^dens$x
        if(popname == "synecho" | popname == "prochloro"){volbins <- round(10^(0.524*log10(norm.fsc) + 0.283),3)
          }else{volbins <- round(10^(0.404*log10(norm.fsc)^2 + 1.802*log10(norm.fsc) + 1.174),3)}



        #######################
        ### SAVE the output ###
        #######################

        colnames(Vhist) <- colnames(Ndist) <- as.character(Time)
        rownames(Vhist) <- rownames(Ndist) <- volbins

        distribution <- list()
            distribution[[1]] <- Vhist
            distribution[[2]] <- Ndist

         print("done")
        return(distribution)

        }







##############################
### PLOT size distribution ###
##############################
plot.size.distribution <- function(distribution, mode = c('log', 'lin'), ...){

    require(rgl)
    jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

    mode <- as.character(mode[1])

    param <- distribution
    percentile <- cut(unlist(param), 100)

    # in linear scale
    if(mode =='lin'){
        plot3d(rep(as.numeric(row.names(param)), dim(param)[2]),
                rep(as.numeric(colnames(param)), each=dim(param)[1]) ,
                unlist(param),
                col=jet.colors(100)[percentile], xlab="size class", ylab="time", zlab="Frequency", ...)
     }


    # in log scale
    if(mode =='log'){
        plot3d(log2(rep(as.numeric(row.names(param)), dim(param)[2])),
                rep(as.numeric(colnames(param)), each=dim(param)[1]) ,
                unlist(param),
                col=jet.colors(100)[percentile],  xlab="size class", ylab="time", zlab="Frequency", ...)
    }
}
