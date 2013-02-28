library(rgl)
library(R.matlab)
library(plotrix)

source('/Volumes/ribalet/Cell_Division/size.class.model_functions.R', chdir = TRUE);mat <- readMat("~/Documents/DATA/SeaFlow/Cell_Division/Matlab/day733320data.mat")

volbins <- mat$volbins
Edata <- mat$Edata
Vhists <- mat$Vhists
N_dist <- mat$N.dist

	resol <- 60 # number of minutes per 1h-interval
	breaks <- 24*60/resol+1
		t <- 0
		dt <- 1/(resol/10)	
		ti <- seq(min(Edata[,1],na.rm=T),max(Edata[,1],na.rm=T), length.out=breaks/dt)
		ep <- data.frame(spline(Edata[,1], Edata[,2], xout=ti)) #interpolate E data according to dt resolution
		Einterp <- ep$y
		Einterp[Einterp < 0] <- 0


####MODEL PARAMETERS

filelist <- list.files("/Volumes/ribalet/Cell_Division/Simulations/",pattern="simulation_");print(paste(length(filelist), "files"))

i <- 1
all.output <- NULL
for(file in filelist){
	Sys.sleep(.1) 
	cat(paste(round(100-100*(length(filelist) - i)/(length(filelist)),2),"% completed\r"))
	i <- i + 1
			
	t <- try(load(paste("/Volumes/ribalet/Cell_Division/Simulations/",file,sep="")))

	if(is.null(dim(proj[[1]]))){ 
			next
		}
		
	re.num <- unlist(strsplit(file, '_'))
	para <- substring(re.num, regexpr('[0-9]',re.num)[1:6])[-1]

	a_ref <- as.numeric(para[1])
	b_ref <- as.numeric(para[2])
	gmax_ref <- as.numeric(para[3])
	E_star_ref <- as.numeric(para[4])
	dmax_ref <- as.numeric(para[5])
	para_i <- data.frame(cbind(gmax_ref,a_ref,b_ref,E_star_ref,dmax_ref))
		
	GR<- sum(proj$mu_N)	
	
		Vproj <- Vhists
		Nproj <- N_dist
		timepoint <- which(apply(N_dist,2,function(x)all(!is.na(x)))) #determine the number of time points with data
		pre.hr <- 0
			for(hr in timepoint[-length(timepoint)]){
					B <- matrix.conct.fast(hr=hr-1, Einterp=Einterp, volbins=volbins, gmax=gmax_ref, a=a_ref, b=b_ref, E_star=E_star_ref,dmax=dmax_ref)						
					Nproj[,hr+1] <- as.vector(round(B %*% Nproj[,hr])) # calculate numbers of individuals
					Vproj[,hr+1] <- as.vector(B %*% Vproj[,hr]) # calculate the projected size-frequency distribution
					Vproj[,hr+1] <- Vproj[,hr+1]/sum(Vproj[,hr+1]) # normalize distribution
					pre.hr <- hr
				}
			mu_N <- diff(colSums(Nproj))/(colSums(Nproj))[-ncol(Nproj)]
			GR_ref <- sum(mu_N, na.rm=T)*(24/length(mu_N))
	
	out <- cbind(para_i, proj$modelresults, GR ,GR_ref)
	all.output <- rbind(all.output, out)
}
	all.output$error <- sqrt((abs(all.output[,1]-all.output[,1+5])/all.output[,1])^2 + (abs(all.output[,2]-all.output[,2+5])/all.output[,2])^2 + (abs(all.output[,3]-all.output[,3+5])/all.output[,3])^2 + (abs(all.output[,4]-all.output[,4+5])/all.output[,4])^2 + (abs(all.output[,5]-all.output[,5+5])/all.output[,5])^2)
	d <- which(all.output$dmax_ref == 0 | all.output$gmax_ref == 0 | all.output$a_ref == 0 |  all.output$b_ref == 0 | all.output$E_star_ref == 0); all.output[d,]
	all.output.1 <- all.output[-d,]
	#nrow(all.output.1)

print(paste(nrow(all.output.1), "valid simulations"))



output <- all.output.1
output <- subset(all.output.1, GR_ref < 2 & GR_ref > 0 & a < 15 & b < 15); print(paste(nrow(output), "realistic simulations"))


d <- which(output$error > 1 | output$resnorm > 2500); output[d,] ### 
output <- output[-d,]


write.csv(output, "Output_simulations1.csv", quote=F, row.names=F)


output <- read.csv("/Volumes/ribalet/Cell_Division/Output_simulationsNoiseF0.csv")
output <- output[-8,]

#jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",	"#FF7F00", "red", "#7F0000"))
jet.colors <- colorRampPalette(c("lightgrey", "pink", "red", "#7F0000","purple"))



#z <- output$resnorm; zlab <- "Difference between Observed and Predicted size distribution"
z <- abs((output$GR_ref-output$GR)/output$GR_ref); zlab <- "GR estimation error"
output <- output[order(z),]


# z <- output$error; zlab <- "Parameter estimation error"
# z <- output$GR_ref; zlab <- "Growth Rate"

par(mfrow=c(3,2), pty='s')
for(p in c("GR","gmax","a","b","E_star","dmax")){
	plot(output[,paste(p,"_ref",sep="")], output[,p], main=paste(p), xlab="true value", ylab="estimated value", pch=16, cex=3,col=jet.colors(100)[cut(z, 100)])#, log='xy')
	points(output[,paste(p,"_ref",sep="")], output[,p], cex=0.5)
	abline(b=1,a=0)
	ylim <- par('usr')[c(3,4)]
    xlim <- par('usr')[c(1,2)]
    
	color.legend(xlim[2]- diff(xlim)/40 , ylim[1], xlim[2], ylim[2], legend=pretty(z), rect.col=jet.colors(100), gradient='y',align='rb',cex=0.75)
	mtext(zlab, side=4, line=3,cex=0.75)  
		}
	


par(mfrow=c(2,1), pty='s')

### DELTA
del <- matrix(nrow=nrow(output), ncol=57)
		for(i in 1:57){
			del[,i] <- output[,"dmax"] * (output[,"a"]*volbins[i])^ output[,"b"] / (1 + ( output[,"a"]*volbins[i])^ output[,"b"])
			}

plot(volbins, del[1,], ylim=c(0,0.02), type='l', col="#00007F", lwd=2, xlab="Cell volume", ylab=paste("Delta (per",10,"min)"))
				for(i in 1:nrow(del))	lines(volbins, del[i,], type='l', col=jet.colors(100)[cut(z, 100)][i], lwd=2)
										ylim <- par('usr')[c(3,4)]
   										xlim <- par('usr')[c(1,2)]
  										color.legend(xlim[2]- diff(xlim)/40 , ylim[1], xlim[2], ylim[2], legend=pretty(z), rect.col=jet.colors(100), gradient='y',align='rb')
												mtext(zlab, side=4, line=3)  

### GAMMA
gamma <- NULL
		for(i in seq(0,1000,10)){
			gam <- output[,"gmax"]*(1-exp(-i/output[,"E_star"]))
			gamma <- cbind(gamma,gam)
			}


		plot(seq(0,1000,by=10),gamma[1,], ylim=c(0,1),type='l', col="#00007F", lwd=2, xlab="Light Intensity", ylab=paste("Gamma (per",10,"min)"))
				for(i in 1:nrow(gamma)) points(seq(0,1000,by=10),gamma[i,],type='l',col=jet.colors(100)[cut(z, 100)][i],lwd=2)
										ylim <- par('usr')[c(3,4)]
   										xlim <- par('usr')[c(1,2)]
  										color.legend(xlim[2]- diff(xlim)/40 , ylim[1], xlim[2], ylim[2], legend=pretty(z), rect.col=jet.colors(100), gradient='y',align='rb')
										mtext(zlab, side=4, line=3)  
		




