library(rgl)
library(zoo)

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",	"#FF7F00", "red", "#7F0000"))

	home <- "~/Documents/DATA/SeaFlow/Thompson/"
	home <- "/Volumes/ribalet/Cell_Division/"
	cruise <- "Rhodomonas_Feb2014" # "Med4_TimeCourse"
	phyto <- "crypto"
	
############################
## OLD SIZE DISTRIBUTION ###
############################	
	
	# Size <- read.csv(paste(home,cruise,"/norm.size.class_",cruise,"_",phyto,".csv",sep=""))
	# Size$time <- as.POSIXct(Size$time, tz="GMT")
	# Size$num.time <- as.numeric(Size$time)
	# Size[Size[,"size.dist"] == 0,"freq.dist"] <- 0
	# n.day <- round(diff(range(Size$time))); print(paste("Number of days in the dataset:",n.day))
	# start <- min(Size$time)

#############################
## FULL SIZE DISTRIBUTION ###
#############################	

list <- list.files(paste(home,cruise,"/",sep=""),pattern=paste("HD.size.class_",cruise,"_",phyto,sep=""))
Size <- NULL

	for(l in list){
		print(l)
		s <- read.csv(paste(home,cruise,"/",l,sep=""))
		Size <- rbind(Size, s)
	}

	Size$time <- as.POSIXct(Size$time, tz="GMT")
	Size$num.time <- as.numeric(Size$time)
	Size[Size[,"size.dist"] == 0,"freq.dist"] <- 0

percentile <- cut(Size[,"freq.dist"], 100); plot3d(x=log10(Size$stages/Size$fsc_beads), y=Size$time, z=Size$freq.dist, col=jet.colors(100)[percentile], type='p', lwd=1, xlab="Size distribution", ylab="Frequency", zlab="Time",axes=F)

percentile <- cut(Size[,"freq.dist"], 100); plot3d(x=Size$stages/Size$fsc_beads, y=Size$time, z=Size$freq.dist, col=jet.colors(100)[percentile], type='p', lwd=1, xlab="Size distribution", ylab="Frequency", zlab="Time",axes=F)



################################
### BINNED SIZE DISTRIBUTION ###
################################
m <- 64

print(paste("phytoplankton population:",phyto))
	
	load(paste(home,cruise,"/", phyto,"_dist_Ncat",m,"_",cruise,sep=""))

	Vhists <- distribution[[1]]
	Vhists <- sweep(Vhists, 2, colSums(Vhists), '/') # Normalize each column of VHists to 1
	N_dist <- distribution[[2]]

	volbins <- as.numeric(row.names(Vhists))
			sizebins <- 2*(volbins*3/(pi*4))^(1/3)# to check the actual diameter

	time.numc <- as.numeric(colnames(Vhists))	
	time <- as.POSIXct(time.numc, origin="1970-01-01" ,tz="GMT")	

para <- Vhists
percentile <- cut(unlist(para), 100); plot3d((rep(as.numeric(row.names(para)), dim(para)[2])), rep(as.numeric(colnames(para)), each=dim(para)[1]) , Vhists , col=jet.colors(100)[percentile], type='l', lwd=6, xlab="size class", ylab="time", zlab="Frequency")
percentile <- cut(unlist(para), 100); plot3d(log10(rep(as.numeric(row.names(para)), dim(para)[2])), rep(as.numeric(colnames(para)), each=dim(para)[1]) , Vhists , col=jet.colors(100)[percentile], type='l', lwd=6, xlab="size class", ylab="time", zlab="Frequency")

head(distribution)
