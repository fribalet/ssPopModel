
run.ssPopModel <- function(path.distribution, Par, time.delay=0, dt=10){

	jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",	"#FF7F00", "red", "#7F0000"))

	require(zoo)


		##################
		### PARAMTERS ###
		##################
		load(path.distribution)
	            t <- as.numeric(time.delay)
		resol <- as.numeric(dt)

		Vhists <- distribution[[1]]
		N_dist <- distribution[[2]]
		volbins <- as.numeric(row.names(Vhists))

		#################################
		### Get the range of 'para' for 'phyto' ###
		#################################

		## Define the time series
		# t <- 1
		time.range <- range(as.numeric(colnames(Vhists)))
		time.interval <- median(diff(as.numeric(colnames(Vhists))))
		days <- as.POSIXct(seq(time.range[1]+t*60*60,time.range[2]+t*60*60 , by=time.interval*24),origin="1970-01-01", tz="GMT") # cut the time series according to time interval

		print(paste("Number of days in the dataset:",length(days)))


		###################
		## LOAD PAR DATA ##
		###################

		Par$time <- as.POSIXct(Par$time, tz="GMT")
		Par[which(Par$par < 0),'par'] <- 0 # remove negative Par values
		id <- which(is.na(Par$par))
			if(length(id) > 0) Par <- Par[-id,] # remove NA in Par

		##########################
		## RUN size.model.functions ##
		##########################

	model <- array(NA, dim=c(4,1))
	for(i in 1:length(days)){

			# i <- 1

			print(paste("24-h growth projection starting at day", i,":",days[i]))

			hours <- seq(days[i], days[i]+60*60*24, by=time.interval)
		#plot(Par$time, Par$par, type='o'); points(c(start, end),c(0,0), col='red',pch=16, cex=2)

		### SELECT SIZE DISTRIBUTION for DAY i
			id <- match(hours, as.numeric(colnames(Vhists)), nomatch=0)
			V.hists <- Vhists[,id]
			N.dist <- N_dist[,id]

			print(paste("the time series has",dim(V.hists)[2] , "/ 25 data points"))
			if(is.null(dim(V.hists))){
				print(paste("Not enough data point, skipping to the next 24-h period"))
				next
				}
			if(dim(V.hists)[2]  < 16){
				print(paste("Not enough data point, skipping to the next 24-h period"))
				next
				}




		    # para <- V.hists; percentile <- cut(unlist(para), 100); plot3d(log2(rep(as.numeric(row.names(para)), dim(para)[2])), rep(as.numeric(colnames(para)), each=dim(para)[1]) , unlist(para), col=jet.colors(100)[percentile], type='l', lwd=3, xlab="size class", ylab="time", zlab="Frequency")


	### SELECT PAR corresponding to V.hists

		light <- subset(Par, time >= hours[1] & time <= hours[25])
		pEdata <- smooth.spline(light[,"time"], light[,"par"])
		Edata <- as.matrix(cbind(pEdata$x, pEdata$y))

	### RUN size.class.model_functions
		proj <- try(.determine.opt.para(V.hists=V.hists,N.dist=N.dist,Edata=Edata, resol=resol))

		if(class(proj) !='try-error'){
		model <- matrix(cbind(as.array(model), as.array(proj)), nrow=4,ncol=ncol(model)+1)
	  save(model, file=paste(home,"/",cruise,"/",phyto,"_model_growth_",cruise,"_Ncat",m,"_t",t, sep=""))

	  }else{print("error during optimization")}
	}

	return(model)
}
