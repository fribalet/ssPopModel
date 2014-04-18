library(rgl)
library(zoo)
library(plotrix)
library(oce)
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",	"#FF7F00", "red", "#7F0000"))
jet.colors <- colorRampPalette(c("blue4","royalblue4","deepskyblue3", "seagreen3", "yellow", "orangered2","darkred")) # Banas


#### PARAMETERS
cat <- 2^6
phyto <- 'crypto'
cruise <- 'Rhodomonas_Feb2014'


## PAR
Par.path <- paste("/Volumes/ribalet/Cell_Division/",cruise,"/Par_",cruise,sep="")
Par <- read.csv(Par.path, sep=",")
Par$time <- as.POSIXct(Par$time, tz='GMT')
Par$UNIXtime <- as.numeric(Par$time)
night <- Par[Par[,"par"] < 2,] ## select Dusk and Dawn time for Day i




## SEAFLOW
all.df <- read.delim(paste("/Volumes/seaflow/",cruise,"/stats.tab",sep=""))
df <- subset(all.df, flag == 0)
df$time <- as.POSIXct(df$time, tz='GMT')
df$UNIXtime <-  as.numeric(df$time)
if(cruise == "Thompson_9"){
	df[which(df$salinity < 1),"salinity"] <- NA
}
pop <- subset(df, pop == phyto) 
if(cruise == 'MBARI_2' & phyto == 'pico') pop <- subset(pop, fsc_small_mode < 30000)






## SIZE DISTRIBUTION
#Size <- read.csv(paste("/Volumes/ribalet/Cell_Division/",cruise,"/size.class_",cruise,"_",phyto,".csv", sep=""))
load(paste("/Volumes/ribalet/Cell_Division/",cruise,"/",phyto,"_dist_Ncat",cat,"_",cruise, sep=""))
Vhists <- distribution[[1]]
N_dist <- distribution[[2]]
volbins <- as.numeric(row.names(Vhists))
	time.numc <- as.numeric(colnames(Vhists))	
	time <- as.POSIXct(time.numc, origin="1970-01-01" ,tz="GMT")	
	n.day <- round(diff(range(time))); print(paste("Number of days in the dataset:",n.day))

# para <- Vhists; percentile <- cut(unlist(para), 100); plot3d(log(rep(as.numeric(row.names(para)), dim(para)[2])), rep(as.numeric(colnames(para)), each=dim(para)[1]) , Vhists , col=jet.colors(100)[percentile], type='l', lwd=6, xlab="size class", ylab="time", zlab="Frequency")
	

## MODEL

all.filelist <- list.files(paste("/Volumes/ribalet/Cell_Division/",cruise,sep=""),pattern=paste(phyto,"_modelHD_growth_",cruise,"_Ncat",cat,sep=""))
filelist <- all.filelist[grep(pattern=paste(phyto), all.filelist)]

n <- 1
Conc.all <- N.proj.all <- V.hist.all <- div.rate <- para.all <-  NULL
for(file in filelist){
	#file <- filelist[3]
	load(paste("/Volumes/ribalet/Cell_Division/",cruise,"/",file, sep=""))
	print(file)
	print(n)
		dim <- conc.proj.all <- n.proj.all <- v.hist.all <- dr.all <- p.all <- NULL
			for(i in seq(2,as.numeric(n.day),by=1)){
				n.proj <- model[4,i][[1]]
				n.proj.all <- cbind(n.proj.all, n.proj)			
							
				conc.proj <- cbind(as.numeric(colnames(n.proj)), as.numeric(colSums(n.proj)))
				conc.proj.all <- rbind(conc.proj.all, conc.proj)
				
				dr <- model[2,i][[1]]
				h.dr <- cbind(as.numeric(colnames(dr)), as.numeric(dr))
				dr.all <- rbind(dr.all, h.dr)
				
				v.proj <- model[3,i][[1]]	
				v.hist.all <- cbind(v.hist.all, v.proj)			

				para <- model[1,i][[1]]
				param <- cbind(time=as.numeric(colnames(n.proj)), para)
				p.all <- rbind(p.all, param)
			}
		
	
		div.rate <- rbind(div.rate, dr.all)
		N.proj.all <- cbind(N.proj.all, n.proj.all)
		Conc.all <- rbind(Conc.all, conc.proj.all)
		V.hist.all <- cbind(V.hist.all, v.hist.all)
		para.all <- rbind(para.all, p.all)

		par(mfrow=c(4,2))
		plot(div.rate, ylab="Div Rate (h-1)", main=paste("file",file))
		plot(Conc.all[,1], Conc.all[,2], ylab="Projected concentration" )
		plot(para.all[,"time"], para.all[,"gmax"], ylab="gmax")
		plot(para.all[,"time"], para.all[,"dmax"],ylab="dmax")
		plot(para.all[,"time"], para.all[,"a"],ylab="a")
		plot(para.all[,"time"], para.all[,"b"],ylab="b")
		plot(para.all[,"time"], para.all[,"E_star"],ylab="E_star")
		plot(para.all[,"time"], para.all[,"resnorm"],ylab="resnorm")

		# 
		# names(para) <- c("gmax","a","b","E_star","dmax","resnorm")
		# par(mfrow=c(4,2))
		# barplot(d.GR, col='grey', main="GR")
		# for(i in 1:6) barplot(para[,i], main=colnames(para)[i])
n <- n + 1	
}

	Div.rate <- div.rate[order(div.rate[,1]),]
	Nproj <- N.proj.all[,order(as.numeric(colnames(N.proj.all)))]
	Conc.proj <- Conc.all[order(Conc.all[,1]),]
	Vproj <- V.hist.all[,order(as.numeric(colnames(V.hist.all)))]
	Para.all <- para.all[order(para.all[,"time"]),]
		
	para <- Vproj 
	 percentile <- cut(unlist(para), 100); plot3d(rep(1:dim(para)[1], dim(para)[2]), rep(1:dim(para)[2], each=dim(para)[1]), z=matrix(para), col=jet.colors(100)[percentile], type='l', lwd=3, xlab="size class", ylab="time", zlab="Frequency")
	

###############
### BINNING ###
###############

		breaks <- seq(min(pop[,"UNIXtime"]),max(pop[,"UNIXtime"]),by=60*60)

	
	# DIVISION 
		h <- cut(Div.rate[,1], breaks=breaks, labels=F)
		h.time.numc <- as.vector(tapply(Div.rate[,1], h, function(x) mean(x, na.rm=T)))#; h.time.numc <- na.approx(h.time.numc, na.rm=F)
		h.dr.mean <- as.vector(tapply(Div.rate[,2], h, function(x) mean(x, na.rm=T)))#; h.dr.mean <- na.approx(h.dr.mean, na.rm=F)
		h.dr.sd <- as.vector(tapply(Div.rate[,2], h, function(x) sd(x, na.rm=T)))#; h.dr.sd <- na.approx(h.dr.sd, na.rm=F)
		#h.time <- as.POSIXct(h.time.numc,origin="1970-01-01",tz='GMT')
		h.time <- as.POSIXct(breaks[findInterval(h.time.numc, breaks)],origin="1970-01-01",tz='GMT')

		id <- findInterval(h.time, pop$time, rightmost.closed=F)
		h.lat <- pop[id,"lat"]
		h.lon<- pop[id,"long"]
		
		D <- data.frame(cbind(h.time, h.lat, h.lon, h.dr.mean, h.dr.sd))
		   
		    plot(h.time, h.dr.mean, ylim=c(0,max(h.dr.mean, na.rm=T)*1.3))		
			abline(v=night$UNIXtime,col='lightgrey')
		    plotCI(h.time, h.dr.mean, h.dr.sd, add=T)		


					
	# SEAFLOW
		if(cruise == "Thompson_9"){pop <- subset(pop, temperature > 10)}

		h2 <- cut(pop[,"UNIXtime"], breaks=breaks, labels=F)
		h2.time.numc <- as.vector(tapply(pop[,"UNIXtime"], h2, function(x) mean(x, na.rm=T)))#; h2.time.numc <- na.approx(h2.time.numc, na.rm=F)
		#h2.time <- as.POSIXct(h2.time.numc,origin="1970-01-01",tz='GMT')
		h2.time <- as.POSIXct(breaks[findInterval(h2.time.numc, breaks)],origin="1970-01-01",tz='GMT')

		h2.conc.mean <- as.vector(tapply(pop[,"conc"], h2, function(x) mean(x, na.rm=T)))#; h2.conc.mean <- na.approx(h2.conc.mean, na.rm=F)
		h2.conc.sd <- as.vector(tapply(pop[,"conc"], h2, sd))#; h2.conc.sd <- na.approx(h2.conc.sd, na.rm=F)
		
		h2.fsc.mean <- as.vector(tapply(pop[,"fsc_small_mode"], h2, function(x) mean(x, na.rm=T)))#; h2.fsc.mean <- na.approx(h2.fsc.mean, na.rm=F)
		h2.fsc.sd <- as.vector(tapply(pop[,"fsc_small_mode"], h2, function(x) sd(x, na.rm=T)))#; h2.fsc.sd <- na.approx(h2.fsc.sd, na.rm=F)
	
		h2.chl.mean <- as.vector(tapply(pop[,"chl_small_mode"], h2, function(x) mean(x, na.rm=T)))#; h2.chl.mean <- na.approx(h2.chl.mean, na.rm=F)
		h2.chl.sd <- as.vector(tapply(pop[,"chl_small_mode"], h2, function(x) sd(x, na.rm=T)))#; h2.chl.sd <- na.approx(h2.chl.sd, na.rm=F)
		
		
			h2.temp.mean <- as.vector(tapply(pop[,"temperature"], h2, function(x) mean(x, na.rm=T)))#; h2.temp.mean <- na.approx(h2.temp.mean, na.rm=F)
			h2.temp.sd <- as.vector(tapply(pop[,"temperature"], h2, function(x) sd(x, na.rm=T)))#; h2.temp.sd <- na.approx(h2.temp.sd, na.rm=F)
			h2.sal.mean <- as.vector(tapply(pop[,"salinity"], h2, function(x) mean(x, na.rm=T)))#; h2.sal.mean <- na.approx(h2.sal.mean, na.rm=F)
			h2.sal.sd <- as.vector(tapply(pop[,"salinity"], h2, function(x) sd(x, na.rm=T)))#; h2.sal.sd <- na.approx(h2.sal.sd, na.rm=F)
			h2.fluo.mean <- as.vector(tapply(pop[,"fluorescence"], h2, function(x) mean(x, na.rm=T)))#; h2.temp.mean <- na.approx(h2.temp.mean, na.rm=F)
			h2.fluo.sd <- as.vector(tapply(pop[,"fluorescence"], h2, function(x) sd(x, na.rm=T)))#; h2.temp.sd <- na.approx(h2.temp.sd, na.rm=F)
		

		S <- data.frame(cbind(h2.time, h2.conc.mean, h2.conc.sd, h2.fsc.mean, h2.fsc.sd, h2.chl.mean, h2.chl.sd, h2.temp.mean, h2.temp.sd ,h2.sal.mean, h2.sal.sd, h2.fluo.mean, h2.fluo.sd))

	par(mfrow=c(4,2))
		    plotCI(h2.time, h2.conc.mean, h2.conc.sd); lines(h2.time, h2.conc.mean,lwd=2,col=2)
		    plotCI(h2.time, h2.fsc.mean, h2.fsc.sd)	
			plotCI(h2.time, h2.chl.mean,h2.chl.sd)
	  		plotCI(h2.time, h2.temp.mean, h2.temp.sd)		
		    plotCI(h2.time, h2.sal.mean, h2.sal.sd)		
			plotCI(h2.time, h2.fluo.mean,h2.fluo.sd)

	# LIGHT
		h4 <- cut(as.numeric(Par$time), breaks=breaks, labels=F)
		h4.time.numc <- as.vector(tapply(as.numeric(Par$time), h4, mean))#; h4.time.numc <- na.approx(h4.time.numc, na.rm=F)
		h4.par.mean <- as.vector(tapply(Par[,"par"], h4, mean))#; h4.par.mean <- na.approx(h4.par.mean, na.rm=F)
		h4.par.sd <- as.vector(tapply(Par[,"par"], h4, sd))#; h4.par.sd <- na.approx(h4.par.sd, na.rm=F)
		#h4.time <- as.POSIXct(h4.time.numc,origin="1970-01-01",tz='GMT')
		h4.time <- as.POSIXct(breaks[findInterval(h4.time.numc, breaks)],origin="1970-01-01",tz='GMT')
		

		L <- data.frame(cbind(h4.time, h4.par.mean, h4.par.sd))

		    plotCI(h4.time, h4.par.mean, h4.par.sd)		
	
	
	DS <- merge(D, S, by.x=c("h.time"), by.y= c("h2.time"),all=T)
	DSL <- merge(DS, L, by.x=c("h.time"), by.y= c("h4.time"),all=T)

	

write.csv(DSL, paste("/Volumes/ribalet/Cell_Division/",cruise,"/",phyto,"_HD_",cruise, ".binned.csv",sep=""),quote=F, row.names=F)






















########################
### MODEL PARAMETERS ###
########################

	par(mfrow=c(3,2))
	for(p in c("h.dr","h2.conc","h2.gr","h2.fsc","h3.lr","h2.chl")){
		plotCI(DSL$h.time, DSL[,paste(p,'.mean',sep="")], DSL[,paste(p,'.sd',sep="")], col=NA, ylab=NA, main=paste(p))
		abline(v=night$UNIXtime,col='lightgrey')
		plotCI(DSL$h.time, DSL[,paste(p,'.mean',sep="")], DSL[,paste(p,'.sd',sep="")],add=T)
	}
	
	
	
	
	
	
	 # PRODUCTION
		# if(phyto == "prochloro" | phyto == "pico") mass <- 53*10^-12 # mg C per cell
		# if(phyto == "synecho") mass <- 250*10^-12
		# prod.mean <- h2.conc.mean*10^6 * (exp(h.dr.mean)-1) * mass * 10^3  / 12 # mg C m-3 h-1
		# prod.sd <- prod.mean * sqrt(((h2.conc.sd/h2.conc.mean)^2 + exp(1)*h.dr.sd))




	# PARAMETERS	
		h5 <- cut(Para.all[,"time"], breaks=breaks)
		h5.time.numc <- as.vector(tapply(Para.all[,"time"], h5, mean))#; h5.time.numc <- na.approx(h5.time.numc, na.rm=F)
		h5.gmax.mean <- as.vector(tapply(Para.all[,"gmax"], h5, mean))#; h5.gmax.mean <- na.approx(h5.gmax.mean, na.rm=F)
		h5.dmax.mean <- as.vector(tapply(Para.all[,"dmax"], h5, mean))#; h5.dmax.mean <- na.approx(h5.dmax.mean, na.rm=F)
		h5.a.mean <- as.vector(tapply(Para.all[,"a"], h5, mean))#; h5.a.mean <- na.approx(h5.a.mean, na.rm=F)
		h5.b.mean <- as.vector(tapply(Para.all[,"b"], h5, mean))#; h5.b.mean <- na.approx(h5.b.mean, na.rm=F)
		h5.E_star.mean <- as.vector(tapply(Para.all[,"E_star"], h5, mean))#; h5.E_star.mean <- na.approx(h5.E_star.mean, na.rm=F)
		h5.resnorm.mean <- as.vector(tapply(Para.all[,"resnorm"], h5, mean))#; h5.resnorm.mean <- na.approx(h5.resnorm.mean, na.rm=F)
		h5.gmax.sd <- as.vector(tapply(Para.all[,"gmax"], h5,sd))#; h5.gmax.sd <- na.approx(h5.gmax.sd, na.rm=F)
		h5.dmax.sd <- as.vector(tapply(Para.all[,"dmax"], h5, sd))#; h5.dmax.sd<- na.approx(h5.dmax.sd, na.rm=F)
		h5.a.sd <- as.vector(tapply(Para.all[,"a"], h5, sd))#; h5.a.sd <- na.approx(h5.a.sd, na.rm=F)
		h5.b.sd <- as.vector(tapply(Para.all[,"b"], h5, sd))#; h5.b.sd <- na.approx(h5.b.sd, na.rm=F)
		h5.E_star.sd <- as.vector(tapply(Para.all[,"E_star"], h5, sd))#; h5.E_star.sd <- na.approx(h5.E_star.sd, na.rm=F)
		h5.resnorm.sd <- as.vector(tapply(Para.all[,"resnorm"], h5, sd))#; h5.resnorm.sd <- na.approx(h5.resnorm.sd, na.rm=F)
		h5.time <- as.POSIXct(h5.time.numc,origin="1970-01-01",tz='GMT')
	
		id <- findInterval(h5.time, pop$time, rightmost.closed=F)
		h5.lat <- pop[id,"lat"]
		h5.lon <- pop[id,"long"]
	
		h5.gamma.mean <- h5.gmax.mean*(1-exp(-h4.par.mean/h5.E_star.mean))
		h5.gamma.sd <- h5.gmax.sd*(1-exp(-h4.par.sd/h5.E_star.sd))
	
	
		del <- matrix(nrow=length(h5.time.numc), ncol=cat)
		for(i in 1:cat){
			del[,i] <- h5.dmax.mean * (h5.a.mean*volbins[i])^h5.b.mean / (1 + (h5.a.mean*volbins[i])^h5.b.mean)
			}
		
	
	


		par(mfrow=c(1,2),mar=c(4,4,4,4))
		plot(volbins, del[1,], ylim=c(0,1), type='l', col="#00007F", lwd=2, xlab="Cell volume", ylab=paste("Delta (per",10,"min)"))
				for(i in 2:nrow(del))	points(volbins, del[i,], type='l', col=jet.colors(nrow(del))[cut(h5.time.numc,nrow(del))][i], lwd=2)
			ylim <- par('usr')[c(3,4)]
   			xlim <- par('usr')[c(1,2)]
  			color.legend(xlim[2]- diff(xlim)/40 , ylim[1], xlim[2], ylim[2], legend=pretty(h5.time.numc), rect.col=jet.colors(100), gradient='y',align='rb')


		plot(seq(0,2000,by=10),h5.gmax.mean[1]*(1-exp(-seq(0,2000,by=10)/h5.E_star.mean[1])), ylim=c(0,1),type='l', col="#00007F", lwd=2, xlab="Light Intensity", ylab=paste("Gamma (per",10,"min)"))
				for(i in 1:length(h5.time)) points(seq(0,2000,by=10),h5.gmax.mean[i]*(1-exp(-seq(0,2000,by=10)/h5.E_star.mean[i])),type='l',col=jet.colors(nrow(del))[cut(h5.time.numc,length(h5.time))][i],lwd=2)
					ylim <- par('usr')[c(3,4)]
   					xlim <- par('usr')[c(1,2)]
  			color.legend(xlim[2]- diff(xlim)/40 , ylim[1], xlim[2], ylim[2], legend=pretty(h5.time.numc), rect.col=jet.colors(100), gradient='y',align='rb')

	
	




