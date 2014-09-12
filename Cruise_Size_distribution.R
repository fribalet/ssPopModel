# [ribalet@bloom Cell_Division]
# for i in $(seq 1 1 6); do echo "Rscript Cruise_Size_distribution.R $i crypto Rhodomonas_Feb2014" | qsub -lwalltime=00:30:00,nodes=1:ppn=1 -N crypto_dist$i -d.; done
# for i in $(seq 1 1 11); do echo "Rscript Cruise_Size_distribution.R $i prochloro Med4_TimeCourse_July2012" | qsub -lwalltime=00:30:00,nodes=1:ppn=1 -N pro_dist$i -d.; done


args <- commandArgs(TRUE)
d <- as.numeric(args[1])
phyto <- as.character(args[2])
cruise <- as.character(args[3])

library(flowPhyto)
library(stats)


###################
### BATCH FILES ###
###################
home <- '~/Cell_Division/'
folder <- NULL
root <- "/misc/seaflow/"

# home <- "/Users/francois/Documents/DATA/SeaFlow/"
# folder <- "Cell_Division/"
# root <- "/Volumes/seaflow/"
# cruise <- "Med4_TimeCourse_July2012"
# phyto <- 'prochloro'
# d <- 1

para <- "fsc_small"
n.breaks <- 2^10 # 1024
concat <- 1 # 3-minute file
out <- NULL


print(paste(cruise," for day",d))
	df <- read.delim(paste(root,cruise,"/", "stats.tab",sep=""))
	df$time <- as.POSIXct(df[,"time"],tz="GMT")
	all.stat <- subset(df, flag == 0) 
	all.stat$time <- as.POSIXct(all.stat[,"time"],tz="GMT")
	all.pop <- subset(all.stat, pop == phyto)

	all.files <- paste(root,cruise,"/",all.pop[,"day"],"/",all.pop[,"file"],".evt.opp", sep="")
	julian.day <- unique(basename(dirname(all.files)))[d]
	filename <- all.files[which(basename(dirname(all.files)) == julian.day)]
	
	
	stat <- subset(all.pop, day == julian.day)


print(paste("Generating Kernel distribution of", phyto, "for cruise",cruise, "at day", julian.day))

		########################
		### INSTRUMENT DRIFT ###
		########################
		spar <- 0.45
		beads <- subset(all.stat, pop == "beads")

		#beads$fsc_small_median <- 10^((beads$fsc_small_mode/2^16)*3.5) # FOR OLD FILES ANALYZED WITHOUT LOG TRANSFORMATION


png(paste(home,folder,cruise,"/Beads_",julian.day,".png",sep=""), width=13.5, height=10, units='in',res=150)
	
		par(mfrow=c(2,1))
		plot(beads[,"time"], beads[,paste(para,"_median",sep="")],pch=1,col='grey', ylim=c(1,10^3.5), xlim=c(min(all.stat[,"time"],na.rm=T),max(all.stat[,"time"],na.rm=T)), xlab='time',ylab=paste("BEADS",para,"_median",sep=""), log='y')
			smooth <- smooth.spline(beads[,"time"], beads[,paste(para,"_median",sep="")], spar=spar)
			time.in.common <- subset(all.pop, time > min(beads[,'time'],na.rm=T) & time < max(beads[,'time'],na.rm=T))
			smooth.beads <- spline(as.POSIXct(smooth$x,origin="1970-01-01",tz="GMT"), smooth$y, xout=as.POSIXct(time.in.common[,"time"],tz="GMT"))
		lines(smooth.beads, col= 'red', lwd=3)	
			fsc_ref <- median(smooth.beads$y)
		abline(h=fsc_ref, lwd=3, lty=2)

		plot(beads[,"time"], beads[,paste(para,"_median",sep="")],pch=1,col='grey', ylim=c(1,10^3.5), xlim=c(min(stat[,"time"],na.rm=T),max(stat[,"time"],na.rm=T)), xlab='time',ylab=paste("BEADS",para,"_median",sep=""), main=paste(julian.day), log='y')
		lines(smooth.beads, col= 'red', lwd=3)	
		abline(h=fsc_ref, lwd=3, lty=2)

dev.off()



		#################################################
		## GENERATE SIZE STRUCTURE FOR EACH TIME STAMP ##
		#################################################

		size.class <- NULL

#f<- filename[1]

			for (f in filename) {
					print(paste("processing", phyto,f))
					opp <- readSeaflow(f, transform=T)
					vct <- try(read.delim(paste(f,".consensus.vct", sep="")),silent=FALSE)
					opp$pop <- vct[,"pop"]
					pop <- subset(opp, pop == phyto)
					time <- df[df[,'day'] == flowPhyto:::.getYearDay(f) & df[,'file'] == getFileNumber(f) & df[,'pop'] == phyto, 'time']
					id.match <- which(as.POSIXct(smooth.beads$x,origin="1970-01-01", tz="GMT") ==  time)
					fsc_beads <- smooth.beads$y[id.match]
						
						if(length(fsc_beads) == 0){
							fsc_beads <- fsc_ref
							print("no beads found")
							}
							
					id <- which(stat[,'time'] == time & stat[,"pop"] == phyto)
					tot <- stat[id, 'evt']
					o.p.p <- stat[id, 'opp']
					n <- stat[id, 'n']
					ntot <- n * (tot/o.p.p)
					time.class <- rep(time, n.breaks)
					
					dens <- density(log10(pop[,para]), n=n.breaks,from=0, to=3.5, kernel='gaussian',na.rm=T)
					freq.dist <- dens$y*diff(dens$x)[1]
					size.dist <- round(freq.dist * ntot)
					stages <- round(10^dens$x,3)
						
					class <- data.frame(cbind(stages=as.numeric(stages), freq.dist=as.numeric(freq.dist), size.dist=as.numeric(size.dist), time=as.character(as.POSIXct(time.class, format="%Y-%m-%d %H:%M:%S",tz="GMT")), fsc_beads=as.numeric(round(fsc_beads))), stringsAsFactors = FALSE)
					#class[,c(1:3,5,6)] <- apply(class[,c(1:3,5,6)],2, as.numeric)
					size.class <- rbind(size.class, class)
							}
					

write.csv(size.class, file=paste(home,folder,cruise,"/HD.size.class_",cruise,"_",phyto,"_",julian.day,".csv", sep=""), row.names=FALSE, quote=FALSE)
print("DONE")	

