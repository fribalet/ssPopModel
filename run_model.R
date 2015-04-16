# [ribalet@bloom Cell_Division]
# for i in $(seq 0 1 24); do echo "Rscript Model_HD_Division_Rate.R $i synecho Thompson_4" | qsub -lwalltime=24:00:00,nodes=1:ppn=1 -N synGR$i -d.; done

#  library(rgl)
library(DEoptim)
library(zoo)

#home <- "/Volumes/ribalet/Cell_division/"; folder <- NULL; cruise <- "Crypto_TimeCourse_June2013"

home <- '~/Cell_Division/'; folder <- NULL

source(paste(home,'functions_modelHD.R',sep=""), chdir = TRUE)

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",	"#FF7F00", "red", "#7F0000"))


args <- commandArgs(TRUE)
t <- as.numeric(args[1])
phyto <- as.character(args[2])
cruise <- as.character(args[3])



###############################
m <- 2^6 # number of size class
###############################




###############
## REFERENCE ##
###############
# library(R.matlab); mat <- readMat("/Users/francois/Documents/DATA/SeaFlow/Cell_Division/Matlab/day733320data.mat"); res <- readMat("/Users/francois/Documents/DATA/SeaFlow/Cell_Division/Matlab/results.mat")
# volbins <- mat$volbins[1,]
# Edata <- mat$Edata
# V.hists <- mat$Vhists
# N.dist <- mat$N.dist
# Vproj <- res$Vproj
# para <- V.hists; percentile <- cut(para, 100); plot3d(log(rep(volbins , breaks)), rep(1:ncol(para), each=nrow(para)), z=matrix(para), col=jet.colors(100)[percentile], type='l', lwd=6)


	##############
	## PAR DATA ##
	##############
	
	Par.path <- paste(home, folder,cruise,"/Par_",cruise,sep="")
	Par <- read.csv(Par.path, sep=",")
	Par$time <- as.POSIXct(Par$time, tz="GMT")
	Par$num.time <- as.numeric(Par$time)


	#######################	
	## SIZE DISTRIBUTION ##
	#######################

	# t <- 0	
	# phyto <- "crypto"
   
    print(paste("time delay:", t))
	print(paste("phytoplankton population:",phyto))
	
	load(paste(home,folder,cruise,"/", phyto,"_dist_Ncat",m,"_",cruise,sep=""))
	Vhists <- distribution[[1]]
	Vhists <- sweep(Vhists, 2, colSums(Vhists), '/') # Normalize each column of VHists to 1
	N_dist <- distribution[[2]]

	volbins <- as.numeric(row.names(Vhists))
			sizebins <- 2*(volbins*3/(pi*4))^(1/3)# to check the actual diameter
	volbins <- volbins/max(volbins) # to make sure values are never > 1, for compatibility issue with the Delta function

	time.numc <- as.numeric(colnames(Vhists))	
	time <- as.POSIXct(time.numc, origin="1970-01-01" ,tz="GMT")	
	n.day <- round(diff(range(time))); print(paste("Number of days in the dataset:",n.day))

	# para <- Vhists; percentile <- cut(unlist(para), 100); plot3d(log(rep(as.numeric(row.names(para)), dim(para)[2])), rep(as.numeric(colnames(para)), each=dim(para)[1]) , Vhists , col=jet.colors(100)[percentile], type='l', lwd=6, xlab="size class", ylab="time", zlab="Frequency")
	
	

	##############################
	## RUN size.model.functions ##
	##############################

	resol <-  60 # number of minutes per interval
	breaks <- 25*60/resol

	model <- array(NA, dim=c(4,1))

# t <- 1

	for(i in seq(1,length(time)-24, 24)){
		print(paste("starting hour:",i+t))
		#i <- 96
		start <- time[i+t]
		end <- time[(i+t)+24]
		
		if(is.na(end)){
			print("cycle is less than 24h")
			next
			}
		print(paste("calculating growth projection from ",start , "to",end))
	
	
	#plot(Par$time, Par$par, type='o'); points(c(start, end),c(0,0), col='red',pch=16, cex=2)

		### SELECT SIZE DISTRIBUTION for DAY i
		V.hists <- Vhists[,c(i:(i+24)+t)]
		N.dist <- N_dist[,c(i:(i+24)+t)]

	    # para <- V.hists; percentile <- cut(unlist(para), 100); plot3d(log(rep(as.numeric(row.names(para))), dim(para)[2]), rep(as.numeric(colnames(para)), each=dim(para)[1]), para , col=jet.colors(100)[percentile], type='l', lwd=6, xlab="size class", ylab="time", zlab="Frequency")


		### SELECT PAR corresponding to each sample
		light <- subset(Par, num.time > start & num.time < end)
		h <- cut(light$num.time, breaks=breaks)
		h.par <- tapply(light$par, h, mean)
		t.Edata <- matrix(cbind(time[c(i:(i+24)+t)], h.par), ncol=2)
        
	        ### NA interpolation
	        Edata <- apply(t.Edata, 2, function(x) na.approx(x, na.rm=F))

		
		### RUN size.class.model_functions
		proj <- try(determine.opt.para(V.hists=V.hists,N.dist=N.dist,Edata=Edata,volbins=volbins))
		
		#para <- proj$Vproj; percentile <- cut(unlist(para), 100); plot3d(log(rep(volbins, 24)), rep(1:ncol(para), each=nrow(para)), z=matrix(para), col=jet.colors(100)[percentile], type='l', lwd=6, xlab="size class", ylab="time", zlab="Frequency")
		
		if(class(proj) !='try-error'){
		model <- matrix(cbind(as.array(model), as.array(proj)), nrow=4,ncol=ncol(model)+1)
	    save(model, file=paste(home,folder,cruise,"/",phyto,"_modelHD_growth_",cruise,"_Ncat",m,"_t",t, sep=""))

	  }else{print("error during optimization")}
}
