# [ribalet@bloom Cell_Division]
# for i in $(seq 0 1 24); do echo "Rscript Model_HD_Division_Rate.R $i synecho Thompson_4" | qsub -lwalltime=24:00:00,nodes=1:ppn=1 -N synGR$i -d.; done

#library(rgl)
library(DEoptim)
library(zoo)

home <- "/Volumes/ribalet/Cell_division/"; cruise <- "DeepDOM"
code <- '~/Documents/DATA/Codes/ssPopModel/'
t <- 0
phyto <- 'prochloro'

source(paste(code,'functions_model.R',sep=""), chdir = TRUE)

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",	"#FF7F00", "red", "#7F0000"))




	##################
	### PARAMTERS ###
	##################
	args <- commandArgs(TRUE)
	t <- as.numeric(args[1])
	phyto <- as.character(args[2])
	cruise <- as.character(args[3])


	m <- 57 # 2^6 # number of size class
	time.interval <- 60 #  number of minutes
	resol <-  10 # number of minutes per interval


	###########################	
	## LOAD SIZE DISTRIBUTION ##
	###########################

	# t <- 0	
	# phyto <- "prochloro"
   	
    	print(paste("time delay:", t))
	print(paste("phytoplankton population:", phyto))
	
	load(paste(home,"/",cruise,"/", phyto,"_dist_Ncat",m,"_",cruise,sep=""))
	Vhists <- distribution[[1]]
	Vhists <- sweep(Vhists, 2, colSums(Vhists), '/') # Normalize each column of VHists to 1
	N_dist <- distribution[[2]]

	volbins <- as.numeric(row.names(Vhists))
		

	#################################
	### Get the range of 'para' for 'phyto' ###
	#################################

	## Define the time series
	# t <- 1 
	time.range <- range(as.numeric(colnames(Vhists)))
	days <- as.POSIXct(seq(time.range[1]+t*60*60,time.range[2]+t*60*60 , by=time.interval*60*24),origin="1970-01-01", tz="GMT") # cut the time series according to time interval

	print(paste("Number of days in the dataset:",length(days)))

	# para <- Vhists; percentile <- cut(unlist(para), 100); plot3d(log(rep(as.numeric(row.names(para)), dim(para)[2])), rep(as.numeric(colnames(para)), each=dim(para)[1]) , as.matrix(para), col=jet.colors(100)[percentile], type='l', lwd=3, xlab="size class", ylab="time", zlab="Frequency")
	
	
	###################
	## LOAD PAR DATA ##
	###################
	
	Par.path <- paste(home,"/", cruise,"/Par_",cruise,sep="")
	Par <- read.csv(Par.path, sep=",")
	Par$time <- as.POSIXct(Par$time, tz="GMT")
	Par[which(Par$par < 0),'par'] <- 0 # remove negative Par values
	id <- which(is.na(Par$par))
		if(length(id) > 0) Par <- Par[-id,] # remove NA in Par

	##########################
	## RUN size.model.functions ##
	##########################

model <- array(NA, dim=c(4,1))
for(i in 8:length(days)){
		
		# i <- 14

		print(paste("24-h growth projection starting at day", i,":",days[i]))
		
		hours <- seq(days[i], days[i]+60*60*24, by=time.interval*60)
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
		if(dim(V.hists)[2]  < 12){
			print(paste("Not enough data point, skipping to the next 24-h period"))
			next
			}




	    # para <- V.hists; percentile <- cut(unlist(para), 100); plot3d(log2(rep(as.numeric(row.names(para)), dim(para)[2])), rep(as.numeric(colnames(para)), each=dim(para)[1]) , unlist(para), col=jet.colors(100)[percentile], type='l', lwd=3, xlab="size class", ylab="time", zlab="Frequency")


	### SELECT PAR corresponding to V.hists
		
		light <- subset(Par, time >= hours[1] & time <= hours[25])
		pEdata <- smooth.spline(light[,"time"], light[,"par"])
		Edata <- as.matrix(cbind(pEdata$x, pEdata$y))

	### RUN size.class.model_functions
		proj <- try(determine.opt.para(V.hists=V.hists,N.dist=N.dist,Edata=Edata,volbins=volbins))
		
		para <- proj$Vproj; percentile <- cut(unlist(para), 100); plot3d(log(rep(volbins, dim(para)[2])), rep(as.numeric(colnames(para)), each=nrow(para)), z=as.matrix(para), col=jet.colors(100)[percentile], type='l', lwd=3, xlab="size class", ylab="time", zlab="Frequency")
		
		if(class(proj) !='try-error'){
		model <- matrix(cbind(as.array(model), as.array(proj)), nrow=4,ncol=ncol(model)+1)
	   # save(model, file=paste(home,"/",cruise,"/",phyto,"_modelHD_growth_",cruise,"_Ncat",m,"_t",t, sep=""))

	  }else{print("error during optimization")}
}
