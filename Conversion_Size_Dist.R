# [ribalet@bloom Cell_Division]
# for i in $(seq 6 1 8); do echo "Rscript Conversion_Size_Dist.R $i crypto Rhodomonas_Feb2014" | qsub -lwalltime=1:00:00,nodes=1:ppn=1 -N crypto_conv$i -d.; done


args <- commandArgs(TRUE)
cat <- as.numeric(args[1])
phyto <- as.character(args[2])
cruise <- as.character(args[3])

home <- '~/Cell_Division/'



#library(rgl)
library(zoo)

# home <- "/Volumes/ribalet/Cell_division/" 
# cruise <- "Thompson_10"
# phyto <- "synecho"



jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",	"#FF7F00", "red", "#7F0000"))


	#######################	
	## SIZE DISTRIBUTION ##
	#######################

	list <- list.files(paste(home,cruise,"/",sep=""),pattern=paste("HD.size.class_",cruise,"_",phyto,sep=""))
	Size <- NULL

	for(l in list){
		print(l)
		s <- read.csv(paste(home,cruise,"/",l,sep=""))
		Size <- rbind(Size, s)
	}

	Size$time <- as.POSIXct(Size$time, tz="GMT")
	Size$num.time <- as.numeric(Size$time)
	
	
	
	# Size$corrected_stages <-  10^((Size$stages/2^16)*3.5)
	# Size$corrected_fsc_beads <- 10^((Size$fsc_beads/2^16)*3.5)
	# if(cruise =="MBARI_1"){
				# Size$corrected_stages <-  10^(((Size$stages+5000)/2^16)*3.5) 
				# Size$corrected_fsc_beads <- 10^((median(Size$fsc_beads)/2^16)*3.5)
				# }
		
	############################### OLD CONVERSION ########################################	
	# Size$volume <-  21.853267*((Size$corrected_stages/Size$corrected_fsc_beads)^1.834432)
	#######################################################################################	
	
	
	
	############################### NEW CONVERSION #################################################################################################
	if(phyto == "synecho" | phyto == "pico" | phyto == "prochloro"){
		Size$volume <- 10^(0.5237*log10(Size$stages/Size$fsc_beads) + 0.2828)
		}
	
	if(phyto == "ultra"){
		Size$volume <- 10^(0.0.5461* log10(Size$stages/Size$fsc_beads)^2 + 1.8712*log10(Size$stages/Size$fsc_beads) + 1.0772)

	}
	
	
	if(phyto == "crypto" | phyto == "nano"){
		Size$volume <- 10^(2.384*log10(Size$stages/Size$fsc_beads) + 1.003)

	}
	
	################################################################################################################################################
	
		#volume.range <- range(Size[which(Size[,"size.dist"] > 10), "volume"]); print(volume.range)
		#volume.range <- range(Size[which(Size[,"freq.dist"] > 10^-2), "volume"]); print(volume.range)
		mean.volume <- median(Size[which(Size[,"freq.dist"] == max(Size[,"freq.dist"])), "volume"]); print(mean.volume)
		mean.diameter <- 2*((mean.volume *3)/(pi*4))^(1/3) ; print(mean.diameter)

		# percentile <- cut(Size[,"freq.dist"], 100); plot3d(x=log(Size$volume), y=Size$num.time, z=Size$freq.dist, col=jet.colors(100)[percentile], type='l', lwd=2)
	
		volume.range <- c(mean.volume/20, mean.volume*20); print(volume.range)
		diameter.range <- 2*((volume.range *3)/(pi*4))^(1/3) ; print(diameter.range)
		
	Size.phyto <- subset(Size, volume > volume.range[1] & volume < volume.range[2])


	# percentile <- cut(Size.phyto[,"freq.dist"], 100); plot3d(x=log(Size.phyto$volume), y=Size.phyto$num.time, z=Size.phyto$freq.dist, col=jet.colors(100)[percentile], type='l', lwd=2)

	n.day <- round(diff(range(Size.phyto$time))); print(paste("Number of days in the dataset:",n.day))
	start <- min(Size.phyto$time)

	##############################	
	## CELL VOLUME DISTRIBUTION ##
	##############################
	
# cat <- 8
	
	###############################
	m <- 2^cat # number of Size class
	###############################
	
	## where to cut Size class
		diff.volume <- log(max(Size.phyto$volume)/min(Size.phyto$volume), base=2)/(m+2)
		volbins.cut.ext <- min(Size.phyto$volume) * 2^((1:(m+3) -1)*diff.volume)
		volbins.cut <- volbins.cut.ext[-c(1, m+3)]
	
	## Size class
		diff.volume <- log(max(Size.phyto$volume)/min(Size.phyto$volume), base=2)/(m-1)
		volbins <- min(Size.phyto$volume) * 2^((1:(m) -1)*diff.volume)


	##############################
	## RUN Size.model.functions ##
	##############################

	resol <-  60 # number of minutes per interval
	hours <- 25
	breaks <- hours*60/resol


		### SELECT Size DISTRIBUTION for DAY i


		### rebuild Size distribution according to volbins
		HD <- cut(Size.phyto$volume, volbins.cut)
		HD.volume <- as.vector(rep(volbins, length(unique(Size.phyto$time))))
		HD.time <- rep(unique(Size.phyto$time), each=m)
		HD.hist <- tapply(Size.phyto$freq.dist, list(HD,Size.phyto$time), mean)
			HD.hist <- as.vector(apply(HD.hist, 2, function(x) na.approx(x, na.rm=F)))
		HD.size <- tapply(Size.phyto$size.dist, list(HD,Size.phyto$time), mean)
			HD.size <- as.vector(apply(HD.size, 2, function(x) na.approx(x, na.rm=F)))

	    # para <- HD.hist; percentile <- cut(para, 100); plot3d(log(HD.volume), HD.time, HD.hist, col=jet.colors(100)[percentile], type='l', lwd=2, xlab="size class", ylab="time", zlab="Frequency")

		Size.volume <- data.frame(cbind(HD.volume,HD.time,HD.hist,HD.size))
		
		### binned the data by 1-h interval
		
		h.time <- as.numeric(seq(min(Size$time), max(Size$time), 60*60))
		h <- cut(Size.volume$HD.time, breaks=h.time, include.lowest = T)
		time <- as.vector(tapply(Size.volume$HD.time, h, mean))
		Vhists <- t(tapply(Size.volume$HD.hist, list(h,Size.volume$HD.volume), mean))
		N_dist <- t(tapply(Size.volume$HD.size, list(h,Size.volume$HD.volume), mean))
		
	        ### NA interpolation
	        Vhists <- try(t(apply(Vhists, 1, function(x) na.approx(x, na.rm=F))))
	        N_dist <- try(t(apply(N_dist, 1, function(x) na.approx(x, na.rm=F))))
	     	
	     	id <- findInterval(h.time, na.approx(time, na.rm=F))

	     	colnames(Vhists) <- colnames(N_dist) <- h.time[id]
	    
	    # para <- Vhists; percentile <- cut(unlist(para), 100); plot3d(log(rep(as.numeric(row.names(para)), dim(para)[2])), rep(as.numeric(colnames(para)), each=dim(para)[1]) , Vhists , col=jet.colors(100)[percentile], type='l', lwd=6, xlab="size class", ylab="time", zlab="Frequency")
	    
	  	    
	    distribution <- list()
	    distribution[[1]] <- Vhists
	    distribution[[2]] <- N_dist
		
	    save(distribution, file=paste(home,cruise,"/", phyto,"_dist_Ncat",m,"_",cruise,sep=""))
	    print(paste("saving ", home,cruise,"/", phyto,"_dist_Ncat",m,"_",cruise,sep=""))