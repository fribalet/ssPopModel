#' Create particle size distribution
#'
#' @param db SQLite3 database file path.
#' @param vct.dir VCT file directory..
#' @param breaks Breaks must be a vector of values defining the breaks for the size distribution.
#' @param quantile OPP Filtering quantile.
#' @return Size distribution
#' @examples
#' \dontrun{
#'  
#' breaks <- # for Qc
#' min <- 0.003 # minimum quotas (3 fgC / cell)
#' delta <- 1/18 # to define the width of each bin
#' m <- 300 # number of bins
#' breaks <- round(min*2^(((1:m)-1)*delta),4) # log2 space bin
#' print(breaks)
#'
#' distribution <- create_PSD(db, vct.dir, breaks, quantile = 50)
#' }
#' @export 
create_PSD <- function(db, vct.dir, breaks, quantile = 50){

  QUANT <- as.numeric(quantile)

  require(popcycle)
  require(tidyverse)

  # Get list of files, with list of outliers
  vct.list <- list.files(vct.dir, "parquet", full.names=T)
  
  ### 1. create PSD for each timepoint 
  i <- 1
  distribution <- tibble()
  for(file.name in vct.list){

    #file.name <- vct.list[2]
    message(round(100*i/length(vct.list)), "% completed \r", appendLF=FALSE)

    ## retrieve PSD data
    # Select data for QUANT
    if(QUANT == 2.5){
      vct <- arrow::read_parquet(file.name, col_select=c("date", c(!contains("diam") & contains('q2.5')))) %>% filter(q2.5)
      vct <- rename(vct, Qc_lwr = Qc_lwr_q2.5, Qc_mid = Qc_mid_q2.5, Qc_upr = Qc_upr_q2.5, pop = pop_q2.5)    
    }
    if(QUANT == 50){ 
      vct <- arrow::read_parquet(file.name, col_select=c("date", c(!contains("diam") & contains('q50')))) %>% filter(q50)
      vct <- rename(vct, Qc_lwr = Qc_lwr_q50, Qc_mid = Qc_mid_q50, Qc_upr = Qc_upr_q50, pop = pop_q50)    
    }
    if(QUANT == 97.5){
      vct <- arrow::read_parquet(file.name, col_select=c("date", c(!contains("diam") & contains('q97.5')))) %>% filter(q97.5)
      vct <- rename(vct, Qc_lwr = Qc_lwr_q97.5, Qc_mid = Qc_mid_q97.5, Qc_upr = Qc_upr_q97.5, pop = pop_q97.5)    
    }

    ## Get particle count in each bin 
    # group by population and by breaks
    # for each refractive index

    psd_lwr <-  vct %>% 
            group_by(date, pop, breaks=cut(Qc_lwr, breaks), .drop=F) %>% 
            count(breaks) %>%
            pivot_wider(names_from = breaks, values_from = n) 
    psd_lwr <- psd_lwr %>% add_column(n='lwr', .after="pop")

    psd_mid <-  vct %>% 
            group_by(date, pop=pop, breaks=cut(Qc_mid, breaks), .drop=F) %>% 
            count(breaks) %>%
            pivot_wider(names_from = breaks, values_from = n) 
    psd_mid <- psd_mid %>% add_column(n='mid', .after="pop")

    psd_upr <- vct %>% 
            filter(pop != "beads") %>%
            group_by(date, pop=pop, breaks=cut(Qc_upr, breaks), .drop=F) %>% 
            count(breaks) %>%
            pivot_wider(names_from = breaks, values_from = n) 
    psd_upr <- psd_upr %>% add_column(n='upr', .after="pop")

    # add data of each refractive index
    psd <- bind_rows(psd_lwr, psd_mid, psd_upr)

    # bind data together
    distribution <- bind_rows(distribution, psd)

    i <- i + 1
    flush.console()
  }

  ### 2. Retrieve metadata
  ## Retrieve SFL table
  sfl <- get.sfl.table(db)
  # format time
  sfl$time <- as.POSIXct(sfl$date, format="%FT%T", tz="UTC")
  # retrieve flow rate (mL min-1) of detectable volume
  fr <- popcycle::flowrate(sfl$stream_pressure, inst= get.inst(db))$flow_rate
  # convert to microL min-1
  fr <- fr * 1000
  # acquisition time (min)
  acq.time <- sfl$file_duration/60
  # volume in microL
  sfl$volume <- round(fr * acq.time , 0)
    
  ## Retrive Outlier table
  outliers <- get.outlier.table(db)
  # merge with sfl
  sfl.all <- merge(sfl, outliers, by="file")

  ## Retrive OPP table
  # retrieve opp/evt
  opp <- as_tibble(get.opp.table(db))
  opp <- opp %>% filter(quantile == QUANT)
    
  ## merge all metadata
  meta <- as_tibble(merge(sfl.all, opp, by="file")[c("time", "lat","lon","volume","opp_evt_ratio","flag")])
    
  ## merge metadata with PSD
  PSD <- as_tibble(merge(distribution, meta, by.x="date",by.y="time",all.x=T)) %>%
          relocate(contains("]"), .after=flag) # reorder column

  ### 3. calculate cell abundance in each bin (cells / microliter)
  ## calculate the volume of virtual core for each population, 
  # volume of virtual core is calculated based on a median value of opp_evt ratio for the entire cruise
  # except for small particles (i.e prochloro and synecho) where it is calcualted based on the opp_evt_ratio at that time
  id <- which(PSD$pop == "prochloro" | PSD$pop == "synecho")
  PSD[id, "volume"] <- PSD[id, "volume"] * PSD[id,"opp_evt_ratio"]
  PSD[-c(id), "volume"] <- PSD[-c(id), "volume"] * median(PSD[["opp_evt_ratio"]][-c(id)])
  
  ## calculate cell abundance in each bin (cells / microliter)
  clmn <- grep("]", names(PSD))
  PSD[,clmn] <- PSD[,clmn] / PSD[["volume"]]

return(PSD)
}


#' Manipulate the size distribution created by FCSplankton::create_PSD(). 
#' Calculate the sum of particles in each size class over specific temporal resolution; transform the header
#'
#' @param distribution Particle size disitribution created by FCSplankton::create_PSD().
#'  i.e., a tibble of size distribution over time. First column must be time (POSIXt class object);
#'  Second column must name of the population; other columns represent the different size classes. 
#'  Size classes can represent either diameter or carbon quota (assuming spherical particles).
#' @param time.step Time step over which to sum the number of particles in each size class. Default 1 hour, must be higher than 3 minutes
#' @param Qc.to.diam Convert carbon quotas to diameter as described in
#'  Menden-Deuer, S. and Lessard, E. J. Carbon to volume relationships for dinoflagellates, diatoms, and other protist plankton.
#'  Limnol. Oceanogr. 45, 569–579 (2000).
#' @param abundance.to.biomass Calculate carbon biomass in each population (i.e. cell abundance x Qc)
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
                                        Qc.to.diam=FALSE, 
                                        interval.to.geomean=FALSE,
                                        abundance.to.biomass=FALSE){
  
  # Check that 'time' is a POSIXt class object 
  if(! lubridate::is.POSIXt(distribution$date)){
  print("Date is not recognized as POSIXt class")
  stop
  }

  # Check that 'pop' column is there 
  if(!any(names(distribution)=='pop')){
    print("column 'pop' is missing")
  stop
  }

   # Calculate the mean in each size class over new time interval
  if(!is.null(time.step)){
    distribution <- distribution %>%
                      group_by(date = cut(date, breaks=time.step), pop) %>%
                      summarise(across(lat:lon, mean), across(volume, sum), across(contains("]"), sum))
  }                    



  # Menden-Deuer, S. & Lessard conversion factors
  d <- 0.261; e <- 0.860
  # convert size interval (factors) into data.frame
  breaks <- strsplit(sub("\\]","",sub("\\(","",colnames(distribution)[clmn])),",")
  # select column that have PSD data
  clmn <- grep("]", names(distribution))

  if(Qc.to.diam){
    #convert Qc into diam using the Menden-Deuer conversion
    b <- lapply(breaks, function(x) round(2*(3/(4*pi)*(as.numeric(x)/d)^(1/e))^(1/3),6))
    colnames(distribution)[clmn] <- sub("\\)","\\]", sub("c","",as.character(b)))
  }

  if(interval.to.geomean){
    # transform size class intervals to mean values (i.e. convert breaks (min, max] to geom mean). 
    if(Qc.to.diam){
      midval <- unlist(list(lapply(b, function(x) sqrt(mean(as.numeric(x))*max(as.numeric(x))))))
    }else{
      midval <- unlist(list(lapply(breaks, function(x) sqrt(mean(as.numeric(x))*max(as.numeric(x))))))
      }
    colnames(distribution)[clmn] <- round(midval,4)
  }
  
  if(abundance.to.biomass){
    # calculate biomass in each bin (pgC / L)
    midval <- unlist(list(lapply(breaks, function(x) sqrt(mean(as.numeric(x))*max(as.numeric(x))))))
    distribution[,clmn] <-  t(diag(midval) %*%  t(as.matrix(distribution[,clmn])))
  }


  # time converted to factor needs to be converted back to POSIXt
  distribution$date <- as.POSIXct(distribution$date, tz='GMT')

  return(distribution)

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

    # select column that have PSD data
    clmn1 <- grep("]", names(distribution))
    clmn2 <- grep("pop", names(distribution))

    distribution <- distribution[,c(clmn1, clmn2, clmn3)]

    # format data to be compatible with scatter3d
    d <- reshape2::melt(distribution, id.vars=c("date","lat","lon","pop"))

    # order data by time
    d <- d[order(d$time),]

    plotly::plot_ly() %>%
          plotly::add_trace(data=d, x= ~ date, y = ~ variable, z = ~ value, 
                            type='scatter3d', mode='lines', line=list(width=lwd), 
                            color=~pop, colors=group.colors) %>%
          plotly::layout(scene = list(xaxis = list(autorange = "reversed"),
                              yaxis = list(title="size classes"),
                              zaxis = list(title="", type= z.type)))
                              #"Carbon (mg L<sup>-1</sup>"
                              #"Abundance (cells µL<sup>-1</sup>"
}
