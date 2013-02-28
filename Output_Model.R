library(rgl)
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",	"#FF7F00", "red", "#7F0000"))
library(zoo)
library(plotrix)
library(oce)

cruise <- 'MBARI_2'

#cruise <- "med4_TimeCourse_July2012"

## PAR
Par.path <- paste("/Volumes/ribalet/Cell_Division/",cruise,"/Par_",cruise,sep="")
Par <- read.csv(Par.path, sep=",")
#unix <- paste(Par$date, Par$time)
#Par$UNIXtime <- as.numeric(strptime(unix, "%d-%m-%Y %H:%M:%S", tz="GMT"))
Par$time <- as.POSIXct(Par$time, tz='GMT')
Par$UNIXtime <- as.numeric(Par$time)
night <- Par[Par[,"par"] < 2,] ## select Dusk and Dawn time for Day i



phyto <- 'synecho'




## SEAFLOW
# all.df <- read.delim(paste("/Volumes/seaflow/",cruise,"/stats.tab",sep="")); df <- subset(all.df, flag == 0)
all.df <- read.delim(paste("~/Desktop/",cruise,"/stats.tab",sep="")); df <- subset(all.df, flag == 0)

df$time <- as.POSIXct(df$time, tz='GMT')
df$UNIXtime <-  as.numeric(df$time)
if(cruise == "Thompson_9"){
	cal <- 9.676766
	df$conc <- df$conc * cal
	df[which(df$salinity < 1),"salinity"] <- NA
}

pop <- subset(df, pop == phyto) 
if(cruise == 'MBARI_2' & phyto == 'pico') pop <- subset(pop, fsc_small_mode < 30000)
beads <- subset(df, pop == 'beads')

## SIZE DISTRIBUTION
Size <- read.csv(paste("/Volumes/ribalet/Cell_Division/",cruise,"/size.class_",cruise,"_",phyto,".csv", sep=""))
Size[Size[,"size.dist"] == 0,"freq.dist"] <- 0
Size$time <- as.POSIXct(Size$time, tz="GMT")
Size$num.time <- as.numeric(Size$time)
volbins <-  10^((unique(Size$stages)/2^16)*3.5)/10 # 2.0 to fit with Sosik size distribution



## MODEL
all.filelist <- list.files(paste("/Volumes/ribalet/Cell_Division/",cruise,sep=""),pattern=paste("model_growth_",cruise,sep=""))
filelist <- all.filelist[grep(pattern=paste(phyto), all.filelist)]

merge <- N.proj.all <- V.hist.all <- div.rate <- para.all <-  NULL
for(file in filelist){
	#file <- filelist[3]
	load(paste("/Volumes/ribalet/Cell_Division/",cruise,"/",file, sep=""))
	n.day <- ncol(model)
	print(file)
	
	d <- 5
	
	if(phyto == "synecho" & cruise == "Thompson_9") if(as.numeric(unlist(strsplit(file, split="_"))[d]) == 2 |
								 as.numeric(unlist(strsplit(file, split="_"))[d]) == 3 |
								 as.numeric(unlist(strsplit(file, split="_"))[d]) == 4 | 
								 as.numeric(unlist(strsplit(file, split="_"))[d]) == 5 |
								 as.numeric(unlist(strsplit(file, split="_"))[d]) == 6 |
							     as.numeric(unlist(strsplit(file, split="_"))[d]) == 15) next

	if(phyto == "pico" & cruise == "MBARI_2") if(as.numeric(unlist(strsplit(file, split="_"))[d]) == 16 |
								as.numeric(unlist(strsplit(file, split="_"))[d]) == 8 | 
								as.numeric(unlist(strsplit(file, split="_"))[d]) == 7 | 
								as.numeric(unlist(strsplit(file, split="_"))[d]) == 9 | 
								as.numeric(unlist(strsplit(file, split="_"))[d]) == 18) next
	
	if(phyto == "synecho" & cruise == "MBARI_2") if(as.numeric(unlist(strsplit(file, split="_"))[d]) == 19 |
								as.numeric(unlist(strsplit(file, split="_"))[d]) == 15) next
	
	if(phyto == "nano" & cruise == "MBARI_2") if(as.numeric(unlist(strsplit(file, split="_"))[d]) == 19	) next

	

	#c(15,10,1,2,3,4,5,6,11,12)) next
		
		dim <- conc.proj.all <- v.hist.all <- dr.all <- p.all <- NULL
			for(i in 1:(n.day)){
				n.proj <- model[4,][[i]]
				conc.proj <- cbind(as.numeric(colnames(n.proj)), as.numeric(colSums(n.proj)))
				conc.proj.all <- rbind(conc.proj.all, conc.proj)
				
				dr <- model[2,][[i]]
				h.dr <- cbind(as.numeric(colnames(dr)[-c(1,2)]), as.numeric(dr[1,-c(1,2)]))
				dr.all <- rbind(dr.all, h.dr)
				
				v.proj <- t(model[3,][[i]])	
				v.hist.all <- rbind(v.hist.all, v.proj)			
			
				para <- model[1,][[i]]
				param <- cbind(time=as.numeric(colnames(n.proj)), para)
				p.all <- rbind(p.all, param)
			}
		
	
		div.rate <- rbind(div.rate, dr.all)
		N.proj.all <- rbind(N.proj.all, conc.proj.all)
		V.hist.all <- rbind(V.hist.all, v.hist.all)
		para.all <- rbind(para.all, p.all)

		par(mfrow=c(4,2))
		plot(N.proj.all[,1], N.proj.all[,2]/1000, ylab="Nproj", main=paste("file", file))
		plot(div.rate, ylab="Div Rate (h-1)", main=paste("file",file))
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
	
}

	Div.rate <- div.rate[order(div.rate[,1]),]
	Vproj.all <- V.hist.all[order(rownames(V.hist.all)),]
	Nproj.all <- N.proj.all[order(N.proj.all[,1]),]
	Para.all <- para.all[order(para.all[,"time"]),]
	
	# para <- Vproj.all
	# percentile <- cut(unlist(para), 100); plot3d(rep(1:nrow(para), 25), rep(1:ncol(para), each=nrow(para)), z=matrix(para), col=jet.colors(100)[percentile], type='l', lwd=3, xlab="size class", ylab="time", zlab="Frequency")




###############
### BINNING ###
###############

		breaks <- seq(min(pop[,"UNIXtime"]),max(pop[,"UNIXtime"]),by=60*60)

	
	# DIVISION 
		h <- cut(Div.rate[,1], breaks=breaks, labels=F)
		h.time.numc <- as.vector(tapply(Div.rate[,1], h, mean))#; h.time.numc <- na.approx(h.time.numc, na.rm=F)
		h.dr.mean <- as.vector(tapply(Div.rate[,2], h, mean))#; h.dr.mean <- na.approx(h.dr.mean, na.rm=F)
		h.dr.sd <- as.vector(tapply(Div.rate[,2], h, sd))#; h.dr.sd <- na.approx(h.dr.sd, na.rm=F)
		#h.time <- as.POSIXct(h.time.numc,origin="1970-01-01",tz='GMT')
		h.time <- as.POSIXct(breaks[findInterval(h.time.numc, breaks)],origin="1970-01-01",tz='GMT')

		id <- findInterval(h.time, pop$time, rightmost.closed=F)
		h.lat <- pop[id,"lat"]
		h.lon<- pop[id,"long"]
		
		D <- data.frame(cbind(h.time, h.lat, h.lon, h.dr.mean, h.dr.sd))
		   
		    plotCI(h.time, h.dr.mean, h.dr.sd)		


		
	# # PROJECTED HIST
		# h0 <- cut(as.numeric(rownames(Vproj.all)), breaks=breaks)
		# h0.time.numc <- as.vector(tapply(as.numeric(rownames(Vproj.all)), h0, mean))#; h0.time.numc <- na.approx(h0.time.numc, na.rm=F)
		# h0.time <- as.POSIXct(h0.time.numc,origin="1970-01-01",tz='GMT')

			# h0.vhist.mean <- NULL
				# for(i in 1:ncol(Vproj.all)){
					# h0.vhist <- t(tapply(Vproj.all[,i], h0, mean))
					# h0.vhist.mean <- rbind(h0.vhist.mean, h0.vhist)
					# }
		# id <- findInterval(h0.time, pop$time, rightmost.closed=F)
		# h0.lat <- pop[id,"lat"]
		# h0.lon <- pop[id,"long"]
		
		

	# PROJECTED CONC
		h1 <- cut(Nproj.all[,1], breaks=breaks, labels=F)
		h1.time.numc <- as.vector(tapply(Nproj.all[,1], h1, mean))#; h1.time.numc <- na.approx(h1.time.numc, na.rm=F)
		h1.nproj.mean <- as.vector(tapply(Nproj.all[,2]/1000, h1, mean))#; h1.nproj.mean <- na.approx(h1.nproj.mean, na.rm=F)
		h1.nproj.sd <- as.vector(tapply(Nproj.all[,2]/1000, h1, sd))#; h1.nproj.sd <- na.approx(h1.nproj.sd, na.rm=F)
		#h1.time <- as.POSIXct(h1.time.numc,origin="1970-01-01",tz='GMT')
		h1.time <- as.POSIXct(breaks[findInterval(h1.time.numc, breaks)],origin="1970-01-01",tz='GMT')
		
		id <- findInterval(h1.time, pop$time, rightmost.closed=F)
		h1.lat <- pop[id,"lat"]
		h1.lon <- pop[id,"long"]

		P <- data.frame(cbind(h1.time, h1.lat, h1.lon, h1.nproj.mean, h1.nproj.sd))

		    plotCI(h1.time, h1.nproj.mean, h1.nproj.sd)		


		
	# SEAFLOW
	
		h2 <- cut(pop[,"UNIXtime"], breaks=breaks, labels=F)
		h2.time.numc <- as.vector(tapply(pop[,"UNIXtime"], h2, mean))#; h2.time.numc <- na.approx(h2.time.numc, na.rm=F)
		#h2.time <- as.POSIXct(h2.time.numc,origin="1970-01-01",tz='GMT')
		h2.time <- as.POSIXct(breaks[findInterval(h2.time.numc, breaks)],origin="1970-01-01",tz='GMT')

		h2.conc.mean <- as.vector(tapply(pop[,"conc"], h2, mean))#; h2.conc.mean <- na.approx(h2.conc.mean, na.rm=F)
		h2.conc.sd <- as.vector(tapply(pop[,"conc"], h2, sd))#; h2.conc.sd <- na.approx(h2.conc.sd, na.rm=F)
		
		h2.fsc.mean <- as.vector(tapply(pop[,"fsc_small_mode"], h2, mean))#; h2.fsc.mean <- na.approx(h2.fsc.mean, na.rm=F)
		h2.fsc.sd <- as.vector(tapply(pop[,"fsc_small_mode"], h2, sd))#; h2.fsc.sd <- na.approx(h2.fsc.sd, na.rm=F)
	
		h2.chl.mean <- as.vector(tapply(pop[,"chl_small_mode"], h2, mean))#; h2.chl.mean <- na.approx(h2.chl.mean, na.rm=F)
		h2.chl.sd <- as.vector(tapply(pop[,"chl_small_mode"], h2, sd))#; h2.chl.sd <- na.approx(h2.chl.sd, na.rm=F)
		
		
			h2.temp.mean <- as.vector(tapply(pop[,"temperature"], h2, mean))#; h2.temp.mean <- na.approx(h2.temp.mean, na.rm=F)
			h2.temp.sd <- as.vector(tapply(pop[,"temperature"], h2, sd))#; h2.temp.sd <- na.approx(h2.temp.sd, na.rm=F)
			h2.sal.mean <- as.vector(tapply(pop[,"salinity"], h2, mean))#; h2.sal.mean <- na.approx(h2.sal.mean, na.rm=F)
			h2.sal.sd <- as.vector(tapply(pop[,"salinity"], h2, sd))#; h2.sal.sd <- na.approx(h2.sal.sd, na.rm=F)
		

	
		# NET GROWTH
		degf <- 24
			s.conc <- smooth.spline(pop[,"time"], pop[,"conc"],df=degf, all.knots=T); plot(pop[,"time"], pop[,"conc"]); lines(s.conc,col='red',lwd=2)
			h2.gr.mean <- as.vector(tapply(20*180*diff(log(s.conc$y))/(as.numeric(diff(s.conc$x))),h2[-1],mean))#; h2.gr.mean <- na.approx(h2.gr.mean, na.rm=F)
			h2.gr.sd <- as.vector(tapply(20*180*diff(log(s.conc$y))/(as.numeric(diff(s.conc$x))),h2[-1],sd))#; h2.gr.sd <- na.approx(h2.gr.sd, na.rm=F) 
			#error.smooth <- abs(residuals(s.conc))/pop[,"conc"]
			#h2.gr.sd <- abs(h2.gr.mean) * as.vector(tapply(error.smooth,h2,mean)); h2.gr.sd <- na.approx(h2.gr.sd, na.rm=F) 

		id <- findInterval(h2.time, pop$time, rightmost.closed=F)
		h2.lat <- pop[id,"lat"]
		h2.lon <- pop[id,"long"]
	
		S <- data.frame(cbind(h2.time, h2.lat, h2.lon, h2.conc.mean, h2.conc.sd, h2.fsc.mean, h2.fsc.sd, h2.chl.mean, h2.chl.sd, h2.temp.mean, h2.temp.sd ,h2.sal.mean, h2.sal.sd, h2.gr.mean, h2.gr.sd))

	
		    plotCI(h2.time, h2.conc.mean, h2.conc.sd)		
		    plotCI(h2.time, h2.fsc.mean, h2.fsc.sd)		
			plotCI(h2.time, h2.gr.mean,h2.gr.sd)

	# LIGHT
		h4 <- cut(Par[,"UNIXtime"], breaks=breaks, labels=F)
		h4.time.numc <- as.vector(tapply(Par[,"UNIXtime"], h4, mean))#; h4.time.numc <- na.approx(h4.time.numc, na.rm=F)
		h4.par.mean <- as.vector(tapply(Par[,"par"], h4, mean))#; h4.par.mean <- na.approx(h4.par.mean, na.rm=F)
		h4.par.sd <- as.vector(tapply(Par[,"par"], h4, sd))#; h4.par.sd <- na.approx(h4.par.sd, na.rm=F)
		#h4.time <- as.POSIXct(h4.time.numc,origin="1970-01-01",tz='GMT')
		h4.time <- as.POSIXct(breaks[findInterval(h4.time.numc, breaks)],origin="1970-01-01",tz='GMT')
		
		id <- findInterval(h4.time, pop$time, rightmost.closed=F)
		h4.lat <- pop[id,"lat"]
		h4.lon <- pop[id,"long"]


		L <- data.frame(cbind(h4.time, h4.lat, h4.lon, h4.par.mean, h4.par.sd))

		    plotCI(h4.time, h4.par.mean, h4.par.sd)		
	
	
	DP <- merge(D, P, by.x=c("h.time","h.lat",'h.lon'), by.y= c("h1.time","h1.lat",'h1.lon'),all=T)
	DPS <- merge(DP, S, by.x=c("h.time","h.lat",'h.lon'), by.y= c("h2.time","h2.lat",'h2.lon'),all=T)
	DPSL <- merge(DPS, L, by.x=c("h.time","h.lat",'h.lon'), by.y= c("h4.time","h4.lat",'h4.lon'),all=T)

	
	
	
	
	# LOSE
		DPSL$h3.lr.mean <- DPSL$h.dr.mean - DPSL$h2.gr.mean
		DPSL$h3.lr.sd <- sqrt((DPSL$h2.gr.sd)^2 + (DPSL$h.dr.sd)^2) # propagation error
		

write.csv(DPSL, paste(phyto,"_",cruise, ".binned.csv",sep=""),quote=F, row.names=F)



################
### PLOTTING ###
################

	par(mfrow=c(3,2))
	for(p in c("h.dr","h2.conc","h2.gr","h2.fsc","h3.lr","h2.chl")){
		plotCI(DPSL$h.time, DPSL[,paste(p,'.mean',sep="")], DPSL[,paste(p,'.sd',sep="")], col=NA, ylab=NA, main=paste(p))
		abline(v=night$UNIXtime,col='lightgrey')
		plotCI(DPSL$h.time, DPSL[,paste(p,'.mean',sep="")], DPSL[,paste(p,'.sd',sep="")],add=T)
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	# # PRODUCTION
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
	
	
		del <- matrix(nrow=length(h5.time.numc), ncol=57)
		for(i in 1:57){
			del[,i] <- h5.dmax.mean * (h5.a.mean*volbins[i])^h5.b.mean / (1 + (h5.a.mean*volbins[i])^h5.b.mean)
			}
		
	
	


		par(mfrow=c(1,2),mar=c(4,4,4,4))
		plot(volbins, del[1,], ylim=c(0,0.1), type='l', col="#00007F", lwd=2, xlab="Cell volume", ylab=paste("Delta (per",10,"min)"))
				for(i in 2:nrow(del))	points(volbins, del[i,], type='l', col=jet.colors(nrow(del))[cut(h5.time.numc,nrow(del))][i], lwd=2)
			ylim <- par('usr')[c(3,4)]
   			xlim <- par('usr')[c(1,2)]
  			color.legend(xlim[2]- diff(xlim)/40 , ylim[1], xlim[2], ylim[2], legend=pretty(h5.time.numc), rect.col=jet.colors(100), gradient='y',align='rb',cex=size)


		plot(seq(0,500,by=10),h5.gmax.mean[1]*(1-exp(-seq(0,500,by=10)/h5.E_star.mean[1])), ylim=c(0,1),type='l', col="#00007F", lwd=2, xlab="Light Intensity", ylab=paste("Gamma (per",10,"min)"))
				for(i in 1:length(h5.time)) points(seq(0,500,by=10),h5.gmax.mean[i]*(1-exp(-seq(0,500,by=10)/h5.E_star.mean[i])),type='l',col=jet.colors(nrow(del))[cut(h5.time.numc,length(h5.time))][i],lwd=2)
					ylim <- par('usr')[c(3,4)]
   					xlim <- par('usr')[c(1,2)]
  			color.legend(xlim[2]- diff(xlim)/40 , ylim[1], xlim[2], ylim[2], legend=pretty(h5.time.numc), rect.col=jet.colors(100), gradient='y',align='rb',cex=size)

	
	



###### PARA

par(mfrow=c(6,1), mar=c(2,4,3,2))

	id <- which(is.na(h5.gmax.mean+h5.gmax.sd))	
	plot(h5.time, h5.gmax.mean,type='n',xlab=NA, ylab=NA,xlim=c(min(pop[,"time"]),max(pop[,"time"])))
	abline(v=night$UNIXtime,col='lightgrey')
	polygon(x=c(h5.time[-id],rev(h5.time[-id])),y=c(h5.gmax.mean[-id]+h5.gmax.sd[-id],rev(h5.gmax.mean[-id]-h5.gmax.sd[-id])), col='darkgrey',border=NA)
	lines(h5.time, h5.gmax.mean,col='black',lwd=2, type='l')
	box(col='black')
	mtext(substitute(paste("gmax")), side=3, line=1, cex=size, las=0)

	id <- which(is.na(h5.dmax.mean+h5.dmax.sd))	
	plot(h5.time, h5.dmax.mean,type='n',xlab=NA, ylab=NA,xlim=c(min(pop[,"time"]),max(pop[,"time"])))
	abline(v=night$UNIXtime,col='lightgrey')
	polygon(x=c(h5.time[-id],rev(h5.time[-id])),y=c(h5.dmax.mean[-id]+h5.dmax.sd[-id],rev(h5.dmax.mean[-id]-h5.dmax.sd[-id])), col='darkgrey',border=NA)
	lines(h5.time, h5.dmax.mean,col='black',lwd=2, type='l')
	box(col='black')
	mtext(substitute(paste("dmax")), side=3, line=1, cex=size, las=0)
	
	id <- which(is.na(h5.a.mean+h5.a.sd))	
	plot(h5.time, h5.a.mean,type='n',xlab=NA, ylab=NA,xlim=c(min(pop[,"time"]),max(pop[,"time"])))
	abline(v=night$UNIXtime,col='lightgrey')
	polygon(x=c(h5.time[-id],rev(h5.time[-id])),y=c(h5.a.mean[-id]+h5.a.sd[-id],rev(h5.a.mean[-id]-h5.a.sd[-id])), col='darkgrey',border=NA)
	lines(h5.time, h5.a.mean,col='black',lwd=2, type='l')
	box(col='black')
	mtext(substitute(paste("a")), side=3, line=1, cex=size, las=0)
	
	id <- which(is.na(h5.b.mean+h5.b.sd))	
	plot(h5.time, h5.b.mean,type='n',xlab=NA, ylab=NA,xlim=c(min(pop[,"time"]),max(pop[,"time"])))
	abline(v=night$UNIXtime,col='lightgrey')
	polygon(x=c(h5.time[-id],rev(h5.time[-id])),y=c(h5.b.mean[-id]+h5.b.sd[-id],rev(h5.b.mean[-id]-h5.b.sd[-id])), col='darkgrey',border=NA)
	lines(h5.time, h5.b.mean,col='black',lwd=2, type='l')
	box(col='black')
	mtext(substitute(paste("b")), side=3, line=1, cex=size, las=0)
	
	id <- which(is.na(h5.E_star.mean+h5.E_star.sd))	
	plot(h5.time, h5.E_star.mean,type='n',xlab=NA, ylab=NA,xlim=c(min(pop[,"time"]),max(pop[,"time"])))
	abline(v=night$UNIXtime,col='lightgrey')
	polygon(x=c(h5.time[-id],rev(h5.time[-id])),y=c(h5.E_star.mean[-id]+h5.E_star.sd[-id],rev(h5.E_star.mean[-id]-h5.E_star.sd[-id])), col='darkgrey',border=NA)
	lines(h5.time, h5.E_star.mean,col='black',lwd=2, type='l')
	box(col='black')
	mtext(substitute(paste("E_star")), side=3, line=1, cex=size, las=0)

	id <- which(is.na(h5.resnorm.mean+h5.resnorm.sd))	
	plot(h5.time, h5.resnorm.mean,type='n',xlab=NA, ylab=NA,xlim=c(min(pop[,"time"]),max(pop[,"time"])))
	abline(v=night$UNIXtime,col='lightgrey')
	polygon(x=c(h5.time[-id],rev(h5.time[-id])),y=c(h5.resnorm.mean[-id]+h5.resnorm.sd[-id],rev(h5.resnorm.mean[-id]-h5.resnorm.sd[-id])), col='darkgrey',border=NA)
	lines(h5.time, h5.resnorm.mean,col='black',lwd=2, type='l')
	box(col='black')
	mtext(substitute(paste("resnorm")), side=3, line=1, cex=size, las=0)