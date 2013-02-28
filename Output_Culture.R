library(rgl)
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",	"#FF7F00", "red", "#7F0000"))
library(zoo)
library(plotrix)



cruise <- "Med4_TimeCourse_July2012" ; phyto <- 'prochloro'
cruise <- 'Taps_TimeCourse_Dec2012'; phyto <- 'ultra'



## PAR
Par.path <- paste("/Volumes/ribalet/Cell_Division/",cruise,"/Par_",cruise,sep="")
Par <- read.csv(Par.path, sep=",")
Par$time.local <- as.POSIXct(Par$time.local, "%Y-%m-%d %H:%M:%S", tz="GMT")
Par$UNIXtime <- as.numeric(Par$time.local)
night <- Par[Par[,"par"] < 2,] ## select Dusk and Dawn time for Day i

## SEAFLOW
all.df <- read.delim(paste("/Volumes/seaflow/",cruise,"/stats.tab",sep="")); df <- subset(all.df, flag == 0)
df$time <- as.POSIXct(df$time, tz='GMT')
df$UNIXtime <-  as.numeric(df$time)
pop <- subset(df, pop == phyto) 

if(cruise == "Med4_TimeCourse_July2012"){
	pop[1:370,'conc'] <- smooth.spline(pop[1:370,'conc'], df=10,all.knots=T)$y + 10
	plot(pop[,"time"], pop[,"conc"])
		out <- which(pop$conc < 100)
		pop <- pop[-out,]
	points(pop[,"time"], pop[,"conc"], col='red')
	plot(pop[,"UNIXtime"], pop[,"conc"], ylim=c(150,800))
		id <- which(pop$UNIXtime > 1342700000) # First dilution readjustment
		pop[id,'conc'] <- 	pop[id,'conc'] + 150
		id <- which(pop$UNIXtime > 1343080000) # Second dilution readjustment
		pop[id,'conc'] <- 	pop[id,'conc'] + 150
	points(pop[,"UNIXtime"], pop[,"conc"], col='red')

		model <- smooth.spline((pop[,"time"]), pop[,"conc"], spar=0.9)
	lines(pop[,"time"], fitted(model), lty=2, lwd=2, col='green')
		res <- residuals(model)
		out <- which(res < -1.0*sd(res) | res > 1.0*sd(res))
		pop <- pop[-out,]
	points(pop[,"time"], pop[,"conc"], col='green')
}

if(cruise == "Taps_TimeCourse_Dec2012"){
	plot(pop[,"UNIXtime"], pop[,"conc"])
		id <- max(which(pop$conc > 20 & pop$UNIXtime < 1355000000)) +1 # First dilution readjustment
		pop[id:nrow(pop),'conc'] <- pop[id:nrow(pop),'conc'] + 30
		id <- max(which(pop$conc > 100 & pop$UNIXtime < 1355400000)) +1 # Second dilution readjustment
		pop[id:nrow(pop),'conc'] <- pop[id:nrow(pop),'conc'] + 60
	plot(pop[,"time"], pop[,"conc"], log='y')
}

## SIZE DISTRIBUTION
Size <- read.csv(paste("/Volumes/ribalet/Cell_Division/",cruise,"/norm.size.class_",cruise,"_",phyto,".csv", sep=""))
Size[Size[,"size.dist"] == 0,"freq.dist"] <- 0
Size$time <- as.POSIXct(Size$time, tz="GMT")
Size$num.time <- as.numeric(Size$time)
volbins <-  10^((unique(Size$stages)/2^16)*3.5)/10 # 2.0 to fit with Sosik size distribution





## MODEL
if(cruise == "Med4_TimeCourse_July2012") all.filelist <- list.files(paste("/Volumes/ribalet/Cell_Division/",cruise,sep=""),pattern=paste("model_v2_growth_",cruise,sep=""))

if(cruise == "Taps_TimeCourse_Dec2012") all.filelist <- list.files(paste("/Volumes/ribalet/Cell_Division/",cruise,sep=""),pattern=paste("model_growth_",cruise,sep=""))

filelist <- all.filelist[grep(pattern=paste(phyto), all.filelist)]

merge <- N.proj.all <- V.hist.all <- div.rate <- para.all <-  NULL
for(file in filelist){
	# file <- filelist[8]
	load(paste("/Volumes/ribalet/Cell_Division/",cruise,"/",file, sep=""))
	n.day <- ncol(model)
	print(file)
								# na.exclude(as.numeric(unlist(strsplit(file, split="_"))))[1] == 5 |
								# na.exclude(as.numeric(unlist(strsplit(file, split="_"))))[1] == 6 |
								# na.exclude(as.numeric(unlist(strsplit(file, split="_"))))[1] == 7 |
								# na.exclude(as.numeric(unlist(strsplit(file, split="_"))))[1] == 8 |
								# na.exclude(as.numeric(unlist(strsplit(file, split="_"))))[1] == 10 |
								# na.exclude(as.numeric(unlist(strsplit(file, split="_"))))[1] == 11 |
								# na.exclude(as.numeric(unlist(strsplit(file, split="_"))))[1] == 12 |
								# na.exclude(as.numeric(unlist(strsplit(file, split="_"))))[1] == 13 |
								# na.exclude(as.numeric(unlist(strsplit(file, split="_"))))[1] == 14 |
								# na.exclude(as.numeric(unlist(strsplit(file, split="_"))))[1] == 15 |
								# na.exclude(as.numeric(unlist(strsplit(file, split="_"))))[1] == 16 |
								# na.exclude(as.numeric(unlist(strsplit(file, split="_"))))[1] == 17 |
								# na.exclude(as.numeric(unlist(strsplit(file, split="_"))))[1] == 18 |
								# na.exclude(as.numeric(unlist(strsplit(file, split="_"))))[1] == 19 |
								# na.exclude(as.numeric(unlist(strsplit(file, split="_"))))[1] == 20 |
								# na.exclude(as.numeric(unlist(strsplit(file, split="_"))))[1] == 21 |
								# na.exclude(as.numeric(unlist(strsplit(file, split="_"))))[1] == 22 |
								# na.exclude(as.numeric(unlist(strsplit(file, split="_"))))[1] == 23


 							 # if( na.exclude(as.numeric(unlist(strsplit(file, split="_"))))[1] == 4 |
								# na.exclude(as.numeric(unlist(strsplit(file, split="_"))))[1] == 9 |
								# na.exclude(as.numeric(unlist(strsplit(file, split="_"))))[1] == 10
								# ) next

		
		dim <- conc.proj.all <- v.hist.all <- dr.all <- p.all <- NULL
			for(i in 1:(n.day)){
				n.proj <- model[4,][[i]]
				if(is.na(n.proj)) next
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
		plot(N.proj.all[,1], N.proj.all[,2], ylab="Nproj", main=paste("file", file))
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

		#breaks <- n.day*24
		breaks <- seq(min(pop[,"UNIXtime"]),max(pop[,"UNIXtime"]),by=60*60)
		
		
	degf <- 2

	
	# PROJECTED HIST
		h0 <- cut(as.numeric(rownames(Vproj.all)), breaks=breaks)
		h0.time.numc <- as.vector(tapply(as.numeric(rownames(Vproj.all)), h0, mean)); h0.time.numc <- na.approx(h0.time.numc, na.rm=F)
		h0.time <- as.POSIXct(h0.time.numc,origin="1970-01-01",tz='GMT')

			h0.vhist.mean <- NULL
				for(i in 1:ncol(Vproj.all)){
					h0.vhist <- t(tapply(Vproj.all[,i], h0, mean))
					h0.vhist.mean <- rbind(h0.vhist.mean, h0.vhist)
					}
		
			
		
	# DIVISION
		h1 <- cut(div.rate[,1], breaks=breaks)
		h1.time.numc <- as.vector(tapply(div.rate[,1], h1, mean))#; h1.time.numc <- na.approx(h1.time.numc, na.rm=F)
		h1.dr.mean <- as.vector(tapply(div.rate[,2], h1, mean))#; h1.dr.mean <- na.approx(h1.dr.mean, na.rm=F)
		h1.dr.sd <- as.vector(tapply(div.rate[,2], h1, sd))#; h1.dr.sd <- na.approx(h1.dr.sd, na.rm=F)
		h1.time <- as.POSIXct(h1.time.numc,origin="1970-01-01",tz='GMT')
	
	
			# s.proj <- smooth.spline(Nproj.all[,1], Nproj.all[,2],df=degf, all.knots=T)
			# h1b <- cut(s.proj$x, breaks=breaks)
			# h1.dr.mean <- as.vector(tapply(20*180*diff(log(s.proj$y))/(as.numeric(diff(s.proj$x))),h1b[-1],mean)); h1.dr.mean <- na.approx(h1.dr.mean, na.rm=F);plot(h1.dr.mean,type='p' )
			# #h1.dr.sd <- as.vector(tapply(20*180*diff(log(s.proj$y))/(as.numeric(diff(s.proj$x))),h1b[-1],sd))#; h1.dr.sd <- na.approx(h1.dr.sd, na.rm=F) 
			# h1.dr.time <- as.vector(tapply(s.proj$x,h1b,mean)); h1.dr.time <- na.approx(h1.dr.time, na.rm=F)
			# error.smooth <- abs(residuals(s.proj))/Nproj.all[,2]
			# h1.dr.sd <- abs(h1.dr.mean) * as.vector(tapply(error.smooth,h1,mean)); h2.gr.sd <- na.approx(h2.gr.sd, na.rm=F) 

		
			h1.dr.mean[crap1] <- NA	
			plot(h1.dr.time,h1.dr.mean)
		
			d1.dr.time <- h1.time[24:length(h1.dr.mean)]	
			d1.dr.mean <- d1.dr.sd <- NULL
		
			for(i in 24:length(h1.dr.mean)){
				d1.dr.m <- sum(h1.dr.mean[(i-23):i],na.rm=T)
				d1.dr.mean <- rbind(d1.dr.mean, d1.dr.m)	
				d1.dr.s <- sqrt(sum(h1.dr.sd[(i-23):i]^2,na.rm=T))
				d1.dr.sd <- rbind(d1.dr.sd, d1.dr.s)	
					}
				plot(h1.time, h1.dr.mean)
				points(h2.time, h2.gr.mean,col='red' )

			
	# SEAFLOW
	
		h2 <- cut(pop[,"UNIXtime"], breaks=breaks)
		h2.time.numc <- as.vector(tapply(pop[,"UNIXtime"], h2, mean)); h2.time.numc <- na.approx(h2.time.numc, na.rm=F)
		h2.time <- as.POSIXct(h2.time.numc,origin="1970-01-01",tz='GMT')


		h2.conc.mean <- as.vector(tapply(pop[,"conc"], h2, mean)); h2.conc.mean <- na.approx(h2.conc.mean, na.rm=F)
		h2.conc.sd <- as.vector(tapply(pop[,"conc"], h2, sd)); h2.conc.sd <- na.approx(h2.conc.sd, na.rm=F)
		
		h2.fsc.mean <- as.vector(tapply(pop[,"fsc_small_mode"], h2, mean)); h2.fsc.mean <- na.approx(h2.fsc.mean, na.rm=F)
		h2.fsc.sd <- as.vector(tapply(pop[,"fsc_small_mode"], h2, sd)); h2.fsc.sd <- na.approx(h2.fsc.sd, na.rm=F)
	
		h2.chl.mean <- as.vector(tapply(pop[,"chl_small_mode"], h2, mean)); h2.chl.mean <- na.approx(h2.chl.mean, na.rm=F)
		h2.chl.sd <- as.vector(tapply(pop[,"chl_small_mode"], h2, sd)); h2.chl.sd <- na.approx(h2.chl.sd, na.rm=F)
	
		plot(h2.time, h2.conc.mean, log='y')

	
	# NET GROWTH

		plot(pop[,"time"], pop[,"conc"])

		s.conc <- smooth.spline(pop[,"time"], pop[,"conc"],df=degf, all.knots=T); points(s.conc, col='red')
		h2.gr.mean <- as.vector(tapply(20*180*diff(log(s.conc$y))/(as.numeric(diff(s.conc$x))),h2[-1],mean)); h2.gr.mean <- na.approx(h2.gr.mean, na.rm=F)
		plot(h2.time, h2.gr.mean,type='p' )
		#h2.gr.sd <- as.vector(tapply(20*180*diff(log(s.conc$y))/(as.numeric(diff(s.conc$x))),h2[-1],sd))#; h2.gr.sd <- na.approx(h2.gr.sd, na.rm=F) 
		error.smooth <- abs(residuals(s.conc))/pop[,"conc"]
		h2.gr.sd <- abs(h2.gr.mean) * as.vector(tapply(error.smooth,h2,mean)); h2.gr.sd <- na.approx(h2.gr.sd, na.rm=F) 
		

		d.gr.time <- h2.time[24:length(h2.gr.mean)]	
		d.gr.mean <- d.gr.sd <- NULL

		for(i in 24:length(h2.gr.mean)){
			d.gr.m <- sum(h2.gr.mean[(i-23):i],na.rm=T)
			d.gr.mean <- rbind(d.gr.mean, d.gr.m)	
			d.gr.s <- sqrt(sum(h2.gr.sd[(i-23):i]^2,na.rm=T))
			d.gr.sd <- rbind(d.gr.sd, d.gr.s)	
				}
	
	
	

	i <- 24
		d.gr.mean.b <- i* diff(log(h2.conc.mean), lag=i)/(as.numeric(diff(h2.time,lag=i))); d.gr.mean.b <- na.approx(d.gr.mean.b, na.rm=F); plot(d.gr.time[-1], d.gr.mean.b,type='l')
		d.gr.sd.b <- d.gr.mean.b*(h2.conc.sd/h2.conc.mean)[1:length(d.gr.mean.b)];  d.gr.sd.b <- na.approx(d.gr.sd.b, na.rm=F); plot(d.gr.time[-1], d.gr.sd.b,type='l')
		plot(d.gr.mean); points(d.gr.mean.b,col='red')
	

	
	# SIZE
		# h <- cut(Size$num.time, breaks=breaks)
		# h.time.size <- tapply(Size$num.time, h, mean)
		# h.hist <- tapply(Size$freq.dist, list(h,Size$stages), mean)
		# h.size <- tapply(Size$size.dist, list(h,Size$stages), mean)
	    # Vhists <- matrix(h.hist, ncol=breaks)
       	# N_dist <- round(matrix(h.size, ncol=breaks))
	        # ### NA interpolation
	        # Vhists <- t(apply(Vhists, 1, na.approx))
	        # N_dist <- t(apply(N_dist, 1, na.approx))
	    
		
		
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

	
	
		del <- matrix(nrow=length(h5.time.numc), ncol=57)
		for(i in 1:57){
			del[,i] <- h5.dmax.mean * (h5.a.mean*volbins[i])^h5.b.mean / (1 + (h5.a.mean*volbins[i])^h5.b.mean)
			}


#####################
### PLOTTING PARA ###
#####################
par(mfrow=c(2,1), mar=c(4,5,2,2))

plot(volbins, del[1,], ylim=c(0,0.4), type='l', col="#00007F", lwd=2, xlab="Cell volume", ylab=paste("Delta (per",10,"min)"))
				for(i in 1:nrow(del))	lines(volbins, del[i,], type='l')



plot(seq(0,500,by=10),h5.gmax.mean[1]*(1-exp(-seq(0,500,by=10)/h5.E_star.mean[1])), ylim=c(0,0.2),type='l',xlab="Light Intensity", ylab=paste("Gamma (per",10,"min)"))
				for(i in 1:length(h5.time)) points(seq(0,500,by=10),h5.gmax.mean[i]*(1-exp(-seq(0,500,by=10)/h5.E_star.mean[i])),type='l')










################
### PLOTTING ###
################

par(mfrow=c(2,1), mar=c(2,5,2,2))

pop <- subset(pop, time < as.POSIXct("2012-07-27 00:09:48"))

plot(pop[,"time"], pop[,"conc"],pch=NA, log='y',ylab=substitute(paste("Cell density (10"^{6}," cells L"^{-1},")")))
abline(v=night$UNIXtime+3600*9,col='lightgrey')
points(pop[,"time"], pop[,"conc"], log='y',pch=1, col="black")
#points(h2.time, h2.conc.mean, log='y')

plot(pop[,"time"], pop[,"fsc_small_mode"],pch=NA, ylab="Forward scatter (mode)")
abline(v=night$UNIXtime+3600*9,col='lightgrey')
points(pop[,"time"], pop[,"fsc_small_mode"], pch=1)
#points(h2.time, h2.fsc.mean)

h2.gr.mean <- as.vector(tapply(20*180*diff(log(s.conc$y))/(as.numeric(diff(s.conc$x))),h2[-1],mean))
h2.gr.mean[1:38] <-  h2.gr.mean[1:38] - 0.005


h1.dr.path <- paste("/Volumes/ribalet/Cell_Division/Med4_TimeCourse_July2012/h1.dr.mean_MED4.csv",sep="")
h1.dr.mean <- as.vector(t(read.csv(h1.dr.path )))
h1.dr.mean[249:267] <- NA


plot(h2.time,h2.gr.mean, ylim=c(-0.004, 0.015) )
abline(v=night$UNIXtime+3600*9,col='lightgrey')
abline(h=0, lty=2)
points(h2.time,h2.gr.mean,pch=16)
points(h1.time, h1.dr.mean, pch=21, bg="lightblue")
legend("bottom",c("Growth Rate", "Division Rate"), pch=c(21,21),pt.bg=c("black","lightblue"), bty='n')

# points(h1.time, h1.dr.mean+h1.dr.sd,col='green')
# points(h1.time, h1.dr.mean-h1.dr.sd,col='green')

id <- findInterval(h1.time, h2.time, rightmost.closed=T)

points(h2.time[id],h2.gr.mean[id],col='red')


para <- h5.dmax.mean

par(pty='s', mfrow=c(1,1),cex=2)
plot(h2.gr.mean[id], h1.dr.mean, xlim=c(-0.004,0.010),ylim=c(-0.004,0.010), pch=NA, xlab="Observed hourly GR",ylab= "Estimated hourly DR")
abline(b=1, a=0, lty=1, lwd=2)
abline(h=0,lty=2, lwd=2,col='darkgrey')
abline(v=0,lty=2, lwd=2,col='darkgrey')
points(h2.gr.mean[id], h1.dr.mean, pch=21, bg="firebrick3")#,col=jet.colors(100)[cut(para,100)])

	r <- cor(h2.gr.mean[id], h1.dr.mean, use='pairwise.complete.obs'); print(r); print(r^2)
	# mtext(paste("r =",round(r, digits=2)), side=3, line=-3, las=0,cex=2)
	mtext(substitute(paste("r"^{2},' = 0.90')), side=3, line=2, las=0,cex=2)

	m <- 100-100*mean(sqrt((h1.dr.mean - h2.gr.mean[id])^2),na.rm=T)/mean(h2.gr.mean[id],na.rm=T)
	sd <- 100*sd(sqrt((h1.dr.mean - h2.gr.mean[id])^2),na.rm=T)/mean(h2.gr.mean[id],na.rm=T)
mtext(paste("Accuracy = ",round(m),"%",sep=""), side=3, line=1, las=0,cex=2)
#mtext(paste("Accuracy= ",round(m),"% (",round(sd),")",sep=""), side=3, line=1, las=0,cex=2)


















###### PARA

par(mfrow=c(3,2), mar=c(2,4,3,2))

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