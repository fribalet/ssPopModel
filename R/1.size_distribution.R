# popcycle.location <- "/Volumes/seaflow/SCOPE_1"
# popname <- "prochloro"
# param <- "fsc_small"
# n.breaks <- 57
# time.interval <- 60 #minutes


#################################
### Get the range of 'param' for 'phyto' ###
#################################
size.distribution <- function(popcycle.location, popname, param="fsc_small", n.breaks=57, time.interval = 60){

    require(popcycle)

        # Get the time range
        print("getting the stat table from the database")
        set.project.location(popcycle.location)
        stat <- get.stat.table()
        phyto.stat <- subset(stat, pop == popname)
        phyto.stat$time <- as.POSIXct(phyto.stat$time,format="%FT%T",tz='GMT')
        time.range <- range(phyto.stat$time)
        time <- seq(time.range[1],time.range[2] , by=60*time.interval) # cut the time series according to time interval

        # Get the range of 'param' for 'phyto' 
        print(paste("obtaining the range in", param, "for", popname, 'be patient, this can take several minutes depending of the amount of particles'))
        param.phyto <- get.opp.by.date(time.range[1], time.range[2], pop=popname, channel=param)
        param.range <- range(param.phyto[,param])


        #########################
        ### SMOOTH bead signal  ###
        #########################

        # Get the beads data
         print(paste("obtaining the median ", param, "of beads for normalization"))
         m.beads <- median(subset(stat, pop =='beads' & time > time.range[1] & time < time.range[2])[,param])
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

             message(round(100*i/length(time)), "% completed \r", appendLF=FALSE)
          
            tryCatch({
            #get the opp for phyto
            t <- as.POSIXct(t, origin="1970-01-01", tz='GMT')
            pop <- try(get.opp.by.date(t, t+60*time.interval, pop=popname, channel=param))

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


        #################################################################
        ### CONVERT normalized forward SCATTER by 1 micron beads to VOLUME ###
        #################################################################       
        print(paste("converting", param, "into volume"))
        norm.fsc <- 2^dens$x
           if(popname == "synecho" | popname == "pico" | popname == "prochloro"){
                volbins <- round(10^(0.524*log10(norm.fsc) + 0.283),3)
                # Size$volume <- 10^(0.5*log10(Size$stages/Size$fsc_beads))# MIE THEORY
                }else{
              #Size$volume <- 10^(0.75*log10(Size$stages/Size$fsc_beads)) # MIE THEORY
              volbins <- round(10^(1.2384*log10(norm.fsc) + 1.003),3)
            }
          


        ####################
        ### SAVE the output ###
        ####################

        colnames(Vhist) <- colnames(Ndist) <- as.character(Time)
        rownames(Vhist) <- rownames(Ndist) <- volbins
           
        distribution <- list()
            distribution[[1]] <- Vhist
            distribution[[2]] <- Ndist

         print("done")   
        return(distribution)

        }







#########################
### SHOW size distribution ###
#########################
plot.size.distribution <- function(distribution, mode = c('log', 'lin')){

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
                col=jet.colors(100)[percentile], type='l', lwd=3, xlab="size class", ylab="time", zlab="Frequency")
     }


    # in log scale
    if(mode =='log'){
        plot3d(log2(rep(as.numeric(row.names(param)), dim(param)[2])), 
                rep(as.numeric(colnames(param)), each=dim(param)[1]) , 
                unlist(param), 
                col=jet.colors(100)[percentile], type='l', lwd=3, xlab="size class", ylab="time", zlab="Frequency")
    }
}

