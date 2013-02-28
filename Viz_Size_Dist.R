library(rgl)
library(zoo)

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",	"#FF7F00", "red", "#7F0000"))

	home <- "~/Documents/DATA/SeaFlow/Thompson/"
	home <- "/Volumes/ribalet/Cell_Division/"
	cruise <- "MBARI_1" # "Med4_TimeCourse"
	phyto <- "nano"
	Size <- read.csv(paste(home,cruise,"/norm.size.class_",cruise,"_",phyto,".csv",sep=""))
	Size$time <- as.POSIXct(Size$time, tz="GMT")
	Size$num.time <- as.numeric(Size$time)
	Size[Size[,"size.dist"] == 0,"freq.dist"] <- 0

	n.day <- round(diff(range(Size$time))); print(paste("Number of days in the dataset:",n.day))
	start <- min(Size$time)
	size <- Size
	
	#Size <- subset(size, time > min(size$time) + 3600*24*1.9 & time < min(size$time) + 3600*24*6.3)
		percentile <- cut(Size[,"freq.dist"], 100); plot3d(x=Size$stages, y=Size$time, z=Size$freq.dist, col=jet.colors(100)[percentile], type='l', lwd=1, xlab="Size distribution", ylab="Frequency", zlab="Time",axes=F)


		
resol <-  60 # number of minutes per interval
hours <- 24
breaks <- 1 + hours*60/resol

t<- 8		
		
		
		size <- subset(Size, time > start + 3600*t & time < start + 3600*hours + 3600*t)
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
	    para <- Vhists; percentile <- cut(unlist(para), 100); plot3d(rep(1:nrow(para), 24), rep(1:ncol(para), each=nrow(para)), z=matrix(para), col=jet.colors(100)[percentile], type='l', lwd=6, xlab="size class", ylab="time", zlab="Frequency")




scatterplot3d(Size$stages, y= as.POSIXct(Size$time, tz="GMT"), z=Size$freq.dist,  xlab="Size distribution", ylab="Frequency", zlab="Time", color=jet.colors(100)[percentile], pch=16, type='h', scale.y=2, grid=F, angle=80)

par.set <-list(axis.line = list(col = "transparent"),clip = list(panel = "off"))
pro <- cloud(Size$freq.dist  ~ Size$stages * as.numeric(Size$time),xlab=list(label="Size distribution",cex=size),ylab=list(label="Frequency",cex=size),zlab=list(label="Time",cex=size),pch=21,type='p',cex=1, lwd=3, scales=list(arrows=T,col='black') ,col='black',fill="skyblue3", screen=list( x=-60, y=-40,z=-20), main="Prochlorococcus", par.settings= par.set,key=list(text=list(title="A"),corner=c(0,1),cex=2))
print(pro, position=c(0, 0.5, 0.5, 1), more=TRUE)

