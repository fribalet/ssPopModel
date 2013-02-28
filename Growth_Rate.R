# [ribalet@bloom Cell_Division]
# for i in $(seq 0 1 24); do echo "Rscript Growth_Rate.R $i synecho MBARI_2" 0 | qsub -lwalltime=06:00:00,nodes=1:ppn=1 -N synGR$i -d.; done
# for i in $(seq 0 1 24); do echo "Rscript Growth_Rate.R $i pico MBARI_2" 0 | qsub -lwalltime=06:00:00,nodes=1:ppn=1 -N picoGR$i -d.; done
# for i in $(seq 0 1 24); do echo "Rscript Growth_Rate.R $i ultra MBARI_2" 0 | qsub -lwalltime=06:00:00,nodes=1:ppn=1 -N ultraGR$i -d.; done
# for i in $(seq 0 1 24); do echo "Rscript Growth_Rate.R $i nano MBARI_2" 0 | qsub -lwalltime=06:00:00,nodes=1:ppn=1 -N nanoGR$i -d.; done

# for i in $(seq 0 1 24); do echo "Rscript Growth_Rate.R $i prochloro Med4_TimeCourse_July2012" | qsub -lwalltime=06:00:00,nodes=1:ppn=1 -N proGR$i -d.; done
# for i in $(seq 0 1 24); do echo "Rscript Growth_Rate.R $i ultra Taps_TimeCourse_Dec2012" | qsub -lwalltime=12:00:00,nodes=1:ppn=1 -N tapsGR$i -d.; done



#  library(rgl)
library(DEoptim)
library(zoo)

#home <- "/Volumes/ribalet/Cell_division/"; folder <- NULL; cruise <- "MBARI_1"; 

home <- '~/Cell_Division/'; folder <- NULL

source(paste(home,'size.class.model_functions.R',sep=""), chdir = TRUE)

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",	"#FF7F00", "red", "#7F0000"))


args <- commandArgs(TRUE)
t <- as.numeric(args[1])
phyto <- as.character(args[2])
cruise <- as.character(args[3])
norm <- as.numeric(args[4])







###############
## REFERENCE ##
###############
# library(R.matlab); mat <- readMat("/Users/francois/Documents/DATA/SeaFlow/Cell_Division/Matlab/day733320data.mat"); res <- readMat("/Users/francois/Documents/DATA/SeaFlow/Cell_Division/Matlab/results.mat")
# volbins <- mat$volbins[1,]
# Edata <- mat$Edata
# Vhists <- mat$Vhists
# N_dist <- mat$N.dist
# Vproj <- res$Vproj
# para <- Vproj; percentile <- cut(para, 100); plot3d(rep(1:nrow(para), breaks), rep(1:ncol(para), each=nrow(para)), z=matrix(para), col=jet.colors(100)[percentile], type='l', lwd=6)


	##############
	## PAR DATA ##
	##############
	
	Par.path <- paste(home, folder,cruise,"/Par_",cruise,sep="")
	Par <- read.csv(Par.path, sep=",")
	#unix <- paste(Par$date, Par$time)
	#Par$as.time <- as.POSIXct(strptime(unix, "%d-%m-%Y %H:%M:%S", tz="GMT"))
	#Par$as.time <- as.POSIXct(strptime(unix, "%Y-%m-%d %H:%M:%S", tz="GMT"))
	Par$as.time <- as.POSIXct(Par$time,tz="GMT")
	Par$num.time <- as.numeric(Par$as.time)

# night <- Par[Par[,"par"] < 2,] ## select Dusk and Dawn time for Day i
# id <- which(diff(night$UNIXtime)> 60*60); id <- c(0, id)
# dawn <- night[id,"UNIXtime"]
# n.day <- length(id); print(paste("number of day:",n.day))

	#######################	
	## SIZE DISTRIBUTION ##
	#######################

	# t <- 4	
	# phyto <- "ultra"
    print(paste("time delay:", t))
	print(paste("phytoplankton population:",phyto))
	if(norm == 1) Size <- read.csv(paste(home,folder,cruise,"/norm.size.class_",cruise,"_",phyto,".csv", sep=""))
	if(norm == 0) Size <- read.csv(paste(home,folder,cruise,"/size.class_",cruise,"_",phyto,".csv", sep=""))
	Size$time <- as.POSIXct(Size$time, tz="GMT")
	Size$num.time <- as.numeric(Size$time)
	Size[Size[,"size.dist"] == 0,"freq.dist"] <- 0

	n.day <- round(diff(range(Size$time))); print(paste("Number of days in the dataset:",n.day))
	start <- min(Size$time)

		#percentile <- cut(Size[,"freq.dist"], 100); plot3d(x=Size$stages, y=Size$num.time, z=Size$freq.dist, col=jet.colors(100)[percentile], type='l', lwd=2)






	##############################
	## RUN size.model.functions ##
	##############################

resol <-  60 # number of minutes per interval
hours <- 24
breaks <- 1 + hours*60/resol

	model_growth <- array(NA, dim=c(4,1))
	
	for(i in 1:(n.day)){

		#i <- 1
		print(paste("calculating growth projection starting ", start+t*3600))
	
	#plot(Par$as.time, Par$par, type='o'); points(c(start+t*3600, start + 3600*hours + 3600*t),c(0,0), col='red',pch=16)

		### SELECT SIZE DISTRIBUTION for DAY i
		#size <- subset(Size, UNIXtime < dawn[i+1]+t & UNIXtime > dawn[i]+t)
		size <- subset(Size, time > start + 3600*t & time < start + 3600*hours + 3600*t)

		#percentile <- cut(size[,"freq.dist"], 100); plot3d(x=size$stages, y=size$num.time, z=size$freq.dist, col=jet.colors(100)[percentile], type='l', lwd=2)

		h <- cut(size$num.time, breaks=breaks)
		h.time.size <- tapply(size$num.time, h, mean)
		h.hist <- t(tapply(size$freq.dist, list(h,size$stages), mean))
		h.size <- t(tapply(size$size.dist, list(h,size$stages), mean))
	    Vhists <- matrix(h.hist, ncol=breaks)
       	N_dist <- round(matrix(h.size, ncol=breaks))
		volbins <-  10^((unique(size$stages)/2^16)*3.5)/10 # 2.0 to fit with Sosik size distribution
		sizebins <- 2*(volbins*3/(pi*4))^(1/3)# to check the actual diameter

	        ### NA interpolation
	        Vhists <- try(t(apply(Vhists, 1, na.approx)))
	        N_dist <- try(t(apply(N_dist, 1, na.approx)))
	    #para <- Vhists; percentile <- cut(unlist(para), 100); plot3d(rep(1:nrow(para), 24), rep(1:ncol(para), each=nrow(para)), z=matrix(para), col=jet.colors(100)[percentile], type='l', lwd=6, xlab="size class", ylab="time", zlab="Frequency")


	    colnames(Vhists) <- colnames(N_dist) <- round(na.approx(h.time.size))
        rownames(Vhists) <- rownames(N_dist) <- round(volbins,5)

	
		### SELECT PAR corresponding to each sample
		light <- subset(Par, num.time < max(size$num.time, na.rm=T) & num.time > min(size$num.time, na.rm=T))
		h <- cut(light$num.time, breaks=breaks)
		h.time.par <- tapply(light$num.time, h, mean)
		h.par <- tapply(light$par, h, mean)
		Edata <- matrix(cbind(h.time.par, h.par), ncol=2)
        
	        ### NA interpolation
	        Edata <- try(t(apply(Edata, 1, na.approx)))

		
		start <- start + 3600*hours


			if(abs(min(light$as.time)- min(size$time)) > resol){
				print(abs(min(light$as.time)- min(size$time)))
				next
			}
			
			if(abs(max(light$as.time)- max(size$time)) > resol){
				print(abs(max(light$as.time)- max(size$time)))
				next
			}
		

		### RUN size.class.model_functions
		proj <- try(determine.opt.para(Vhists=Vhists,N_dist=N_dist,Edata=Edata,volbins=volbins))
		
		#para <- proj$Vproj; percentile <- cut(unlist(para), 100); plot3d(rep(1:nrow(para), 24), rep(1:ncol(para), each=nrow(para)), z=matrix(para), col=jet.colors(100)[percentile], type='l', lwd=6, xlab="size class", ylab="time", zlab="Frequency")
		
		
		if(class(proj) !='try-error'){
		model_growth <- matrix(cbind(as.array(model_growth), as.array(proj)), nrow=4,ncol=i+1)
		model <- model_growth[,-1]
	    	if(norm == 1) save(model, file=paste(home,folder,cruise,"/model_norm_growth_",cruise,"_",t,"_", phyto, sep=""))
	    	if(norm == 0) save(model, file=paste(home,folder,cruise,"/model_growth_",cruise,"_",t,"_", phyto, sep=""))

	  }else{print("error during optimization")}
}
