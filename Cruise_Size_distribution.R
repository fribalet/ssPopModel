# [ribalet@bloom Cell_Division]
# echo "Rscript Cruise_Size_distribution.R $i synecho MBARI_2" | qsub -lwalltime=24:00:00,nodes=1:ppn=1 -N syn_size_dist$i -d.
# echo "Rscript Cruise_Size_distribution.R $i pico MBARI_2" | qsub -lwalltime=24:00:00,nodes=1:ppn=1 -N pico_size_dist$i -d.
# echo "Rscript Cruise_Size_distribution.R $i ultra MBARI_2" | qsub -lwalltime=24:00:00,nodes=1:ppn=1 -N ultra_size_dist$i -d.
# echo "Rscript Cruise_Size_distribution.R $i nano MBARI_2" | qsub -lwalltime=24:00:00,nodes=1:ppn=1 -N nano_size_dist$i -d.



args <- commandArgs(TRUE)
phyto <- as.character(args[1])
cruise <- as.character(args[2])

library(flowPhyto)


jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",	"#FF7F00", "red", "#7F0000"))

###################
### BATCH FILES ###
###################
#home <- "/Users/francois/Documents/DATA/SeaFlow/"
home <- '~/Cell_Division/'

#folder <- "Cell_Division/"
folder <- NULL

#root <- "/Volumes/seaflow/"
root <- "/misc/seaflow/"

para <- "fsc_small"
n.breaks <- 57
concat <- 1 # 3-minute file
min.phyto <- 30 # minimum number of individuals
spar <- 0.75
out <- NULL


	print(cruise)
	df <- read.delim(paste(root,cruise,"/", "stats.tab",sep=""))
	df$time <- as.POSIXct(df[,"time"],tz="GMT")
	stat <- subset(df, flag == 0) 
	stat$time <- as.POSIXct(stat[,"time"],tz="GMT")
	
	
		print(paste("Generating Kernel distribution of", phyto, "for cruise",cruise))
		ref <- subset(stat, pop == phyto)
		filename <- paste(root,cruise,"/",ref[,"day"],"/",ref[,"file"],".evt.opp", sep="")

		size.min <- min(smooth.spline(ref[,paste(para,"_mode",sep="")] - ref[,paste(para,"_width",sep="")])$y , na.rm=T); if(size.min < 0) size.min <- 0
		size.max <- max(smooth.spline(ref[,paste(para,"_mode",sep="")] + ref[,paste(para,"_width",sep="")])$y , na.rm=T);	if(size.max > 2^16) size.max <- 2^16
		
			plot(1:nrow(ref), 1:nrow(ref), ylim=c(0,2^16), pch=NA)
				points(ref[,paste(para,"_mode",sep="")] + ref[,paste(para,"_width",sep="")], col='grey')
				points(ref[,paste(para,"_mode",sep="")], col='red')
				points(ref[,paste(para,"_mode",sep="")] - ref[,paste(para,"_width",sep="")], col='grey')
				lines(smooth.spline(ref[,paste(para,"_mode",sep="")] + ref[,paste(para,"_width",sep="")])$y)
				lines(smooth.spline(ref[,paste(para,"_mode",sep="")])$y)
				lines(smooth.spline(ref[,paste(para,"_mode",sep="")] - ref[,paste(para,"_width",sep="")])$y)

		########################
		### INSTRUMENT DRIFT ###
		########################
		beads <- subset(stat, pop == "beads" & fsc_small_mode < 60000)
		plot(beads[,"time"], beads[,paste(para,"_mode",sep="")],pch=1,col='grey', ylim=c(0,2^16), xlim=c(min(stat[,"time"],na.rm=T),max(stat[,"time"],na.rm=T)), xlab='time',ylab=paste("BEADS",para,"_mode",sep=""))
		smooth <- smooth.spline(beads[,"time"], beads[,paste(para,"_mode",sep="")], spar=spar)
		time.in.common <- subset(stat, time > min(beads[,'time'],na.rm=T) & time < max(beads[,'time'],na.rm=T))
		
		smooth.beads <- spline(as.POSIXct(smooth$x,origin="1970-01-01",tz="GMT"), smooth$y, xout=as.POSIXct(time.in.common[,"time"],tz="GMT"))
				lines(smooth.beads, col= 'red', lwd=3)	

		fsc_ref <- median(smooth.beads$y)
		abline(h=fsc_ref, lwd=3, lty=2)



		#################################################
		## GENERATE SIZE STRUCTURE FOR EACH TIME STAMP ##
		#################################################

		size.class <- NULL
		norm.size.class <- NULL

#f<- filename[134]

			for (f in filename) {
					print(paste("processing", phyto,f))
					opp <- readSeaflow(f)
					vct <- try(read.delim(paste(f,".consensus.vct", sep="")),silent=FALSE)
					opp$pop <- vct[,"pop"]
					pop <- subset(opp, pop == phyto)
					time <- df[df[,'day'] == flowPhyto:::.getYearDay(f) & df[,'file'] == getFileNumber(f) & df[,'pop'] == phyto, 'time']
					if(nrow(pop) > min.phyto){
						id.match <- which(as.POSIXct(smooth.beads$x,origin="1970-01-01", tz="GMT") ==  time)
						fsc_beads <- smooth.beads$y[id.match]
						
						if(length(fsc_beads) == 0){
							fsc_beads <- fsc_ref
							print("no beads")
						}
					
					drift.beads <- fsc_ref/fsc_beads
					
					id <- which(stat[,'time'] == time & stat[,"pop"] == phyto)
					tot <- stat[id, 'evt']
					o.p.p <- stat[id, 'opp']
					n <- stat[id, 'n']
					ntot <- n * (tot/o.p.p)
					time.class <- rep(time, n.breaks)
					
					
					i <- 1
					for(drift in c(1, drift.beads)){
						dens <- density(pop[,para]*drift, n=n.breaks,from=size.min, to=size.max, bw="SJ", kernel='gaussian',na.rm=T)
						freq.dist <- dens$y*diff(dens$x)[1]
						size.dist <- round(freq.dist * ntot)
						stages <- round(dens$x)
						
						if(i ==1){
							class <- data.frame(cbind(stages=as.numeric(stages), freq.dist=as.numeric(freq.dist), size.dist=as.numeric(size.dist), time=as.character(as.POSIXct(time.class, format="%Y-%m-%d %H:%M:%S",tz="GMT")), normalization=as.numeric(drift)),stringsAsFactors = FALSE)
							class[,c(1:3,5)] <- apply(class[,c(1:3,5)],2, as.numeric)
							size.class <- rbind(size.class, class)
								}
						
						if(i ==2){
							norm.class <- data.frame(cbind(stages=as.numeric(stages), freq.dist=as.numeric(freq.dist), size.dist=as.numeric(size.dist), time=as.character(as.POSIXct(time.class, format="%Y-%m-%d %H:%M:%S",tz="GMT")), normalization=as.numeric(drift)),stringsAsFactors = FALSE)
							norm.class[,c(1:3,5)] <- apply(norm.class[,c(1:3,5)],2, as.numeric)
							norm.size.class <- rbind(norm.size.class, norm.class)
								}
						
						i <- i+1									
						}
						
					}else print(paste("no enough", phyto, "in",f))
					
				}
			
					
					write.csv(size.class, file=paste(home,folder,cruise,"/size.class_",cruise,"_",phyto,".csv", sep=""), row.names=FALSE, quote=FALSE)
					write.csv(norm.size.class, file=paste(home,folder,cruise,"/norm.size.class_",cruise,"_",phyto,".csv", sep=""), row.names=FALSE, quote=FALSE)