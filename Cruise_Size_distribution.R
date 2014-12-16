#[gwennm@bloom DeepDOM/Cell_Division]
##Run this script using the commented out code just below (paste directly into the command line), should be in the directory listed above
# Rscript ~/DeepDOM/scripts/Cruise_Size_distribution_sqlite_GMH.R d1 d2 prochloro DeepDOM 

#note edit line above according to changes to script: make day input an integer number range (ie: 1,2 for 1:2).


#allowing arguments to be input from the command line input commented out above
args <- commandArgs(TRUE)
d1 <- as.numeric(args[1]) #first argument day number to start
d2 <- as.numeric(args[2]) #2nd argument, day number to end
phyto <- as.character(args[3]) # 2nd arg phyto to calc size dist
cruise <- as.character(args[4]) # 3rd arg cruise name


library(popcycle)
library(stats)

#########################
### BATCH FILE inputs ###
#########################
#globals necessary for running on bloom
# home <- '~/DeepDOM/Cell_Division/' #change to take from input, so not hardcoded in
# folder <- NULL
# root <- "/misc/seaflow/"
#db.location <- "~/popcycle" #change to make variable?

#globals necessary for running on local machine connected to bloom
home <- "/Users/gwen/Desktop/Cruises/DeepDOM_2013/seaflow/"
folder <- "Cell_Division/"
root <- "/Volumes/seaflow/"
db.location <- "/Volumes/gwennm/popcycle"
d1 <- 1 #start day
d2 <- 2 # end day
phyto <- 'prochloro'
cruise <- "DeepDOM"

#Globals necessary for popcycle commands: consider making an input variable to the script
set.evt.location(paste(root, cruise, sep=""))
set.project.location(db.location)
set.cruise.id("march2013")

#parameters for making smoothed distributions
para <- "fsc_small"
n.breaks <- 2^10 # 1024
#concat <- 4 # note this does not work here, make as an option later?
out <- NULL


# new function to query sqlite database by parameter and population and file.name
get.param.by.pop <- function(file.name, para, phyto, db= db.name){
	#this sqlite query will return a table with columns: file, param, pop (only phyto)
	sql <- paste0("SELECT opp.file, opp.", para, ", vct.pop
		FROM opp, vct 
		WHERE opp.cruise == vct.cruise 
		AND opp.file == vct.file
		AND opp.file == '", file.name, "'
		AND opp.particle == vct.particle
		AND vct.pop == '", phyto,"'",
		";")
	#to only subset good files, should add flag txt file into the db as a separate table
	con <- dbConnect(SQLite(), dbname = db)
	opp.slice <- dbGetQuery(con, sql)
	dbDisconnect(con)
	return(opp.slice)
}

#new function to retrieve opp.evt.table from sqlite database, suggest adding to popcycle
get.opp.evt.ratio.table <- function(db=db.name){
	sql <- paste0("SELECT * FROM ", opp.evt.ratio.table.name)
	con <- dbConnect(SQLite(), dbname= db)
	table <- dbGetQuery(con, sql)
	dbDisconnect(con)
	return (table)
}


###########################
### Loading stats table ###
### keep only good data ###
###########################
#note: new input allows stats table to be loaded once and all days run in loop in R
df <- get.stat.table()
df$time <- strptime(df$time, "%Y-%m-%dT%H:%M:%S", tz="GMT")

#load flag file and pare down stats to only good files
flag <- read.table(paste(project.location,"/sqlite/flag_file.txt", sep=""), header=T)
flag.good <- subset(flag, flag ==0)
all.stat <- subset(df, df$file %in% flag.good$file) 

#pare down stats again to only include population of interest
all.pop <- subset(all.stat, pop == phyto)
all.files <- all.pop$file
julian.day <- unique(basename(dirname(all.files)))

#load opp.evt.ratios to correct for true cell density, so far not necessary because I use abundance from stats table
#oer <- get.opp.evt.ratio.table()


		########################
		### INSTRUMENT DRIFT ###
		########################
		#need this code to normalize to beads fsc_small across the whole cruise
		
		spar <- 0.45
		beads <- subset(all.stat, pop == "beads" & fsc_small < 60000)

tiff(paste(home,folder,"/Beads_",cruise,".tiff",sep=""),compression="lzw", width=13.5, height=10, units='in',res=150)
	
		plot(beads[,"time"], beads[,para],pch=1,col='grey', ylim=c(1,10^3.5), xlab='time',ylab=paste("BEADS",para,"_median",sep=""), log='y')
			smooth <- smooth.spline(beads[,"time"], beads[,para], spar=spar)
			time.in.common <- subset(all.pop, time > min(beads[,'time'],na.rm=T) & time < max(beads[,'time'],na.rm=T))
			smooth.beads <- spline(as.POSIXct(smooth$x,origin="1970-01-01",tz="GMT"), smooth$y, xout=as.POSIXct(time.in.common[,"time"],tz="GMT"))
		lines(smooth.beads, col= 'red', lwd=3)	
			fsc_ref <- median(smooth.beads$y)
		abline(h=fsc_ref, lwd=3, lty=2)

dev.off()


		#################################################
		## GENERATE SIZE STRUCTURE FOR EACH TIME STAMP ##
		## new file for all the days in range (d1,d2)  ##
		#################################################

for(n in d1:d2){
	filenames <- all.files[which(basename(dirname(all.files)) == julian.day[n])]#list all filenames from the day 'n'
	stat <- subset(all.pop, file %in% filenames)# subset stats for only day 'n'
	print(paste0("Generating size distribution for day ", n,": ", julian.day[n]))
	
	size.class <- NULL

	for(file in filenames){
	#still need to figure out how to implement concat to put together 12 min chunks
	#do we need to concat at this step? or can it wait until the model?
		slice <- get.param.by.pop(file, para, phyto)
		
		#find smoothed beads from same time stamp to normalize param of interest
		time <- stat[which(stat$file== file),"time"]
		id.match <- which(as.POSIXct(smooth.beads$x,origin="1970-01-01", tz="GMT") ==  time)
		fsc_beads <- smooth.beads$y[id.match]
		if(length(fsc_beads) == 0){
			fsc_beads <- fsc_ref
			print("no beads found")
		}
		
		#insert total abundance of this pop from the correct file, ** is this correct??
		ntot <- stat[which(stat$file == file), "abundance"]
		time.class <- rep(time, n.breaks) #make vector of time to go with dist
		
		#Note changed the start range of density from 1 to 0 for smaller pro
		#shouldn't this range be dependent on the phyto we are trying to model? maybe this should be a range determined by the max of the slice
		#could add an if statement, if phyto = prochloro then this range for dens
		dens <- density(log10(slice[,para]), n=n.breaks,from=0, to=3.5, bw="SJ", kernel='gaussian',na.rm=T)
		freq.dist <- dens$y*diff(dens$x)[1]
		size.dist <- round(freq.dist * ntot)
		stages <- round(10^dens$x,3)
						
		class <- data.frame(cbind(stages=as.numeric(stages), freq.dist=as.numeric(freq.dist), size.dist=as.numeric(size.dist), time=as.character(as.POSIXct(time.class, format="%Y-%m-%d %H:%M:%S",tz="GMT")), fsc_beads=as.numeric(round(fsc_beads))), stringsAsFactors = FALSE)
		size.class <- rbind(size.class, class)
				
	}
	
	write.csv(size.class, file=paste(home,folder,"/HD.size.class_",cruise,"_",phyto,"_",julian.day[n],".csv", sep=""), row.names=FALSE, quote=FALSE)
	print("DONE")
}

	

	

