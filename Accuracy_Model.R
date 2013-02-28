# for i in $(seq 0.001 0.01 0.1); do echo "Rscript Accuracy_Model.R $i 1.23 3.77 124 0.03" | qsub -lwalltime=02:00:00,nodes=1:ppn=1 -N Gmax$i -d.; done
# for i in $(seq 0.01 0.1 1); do echo "Rscript Accuracy_Model.R 0.14 $i 3.77 124 0.03" | qsub -lwalltime=02:00:00,nodes=1:ppn=1 -N a$i -d.; done
# for i in $(seq 0.01 0.1 1); do echo "Rscript Accuracy_Model.R 0.14 1.23 $i 124 0.03" | qsub -lwalltime=02:00:00,nodes=1:ppn=1 -N b$i -d.; done
# for i in $(seq 0.01 1 10); do echo "Rscript Accuracy_Model.R 0.14 1.23 3.77 $i 0.03" | qsub -lwalltime=02:00:00,nodes=1:ppn=1 -N Estar$i -d.; done
# for i in $(seq 0.001 0.01 0.1); do echo "Rscript Accuracy_Model.R 0.14 1.23 3.77 124 $i" | qsub -lwalltime=02:00:00,nodes=1:ppn=1 -N Dmax$i -d.; done

# for i in $(seq 500) ;do echo "Rscript Accuracy_Model.R" | qsub -lwalltime=02:00:00,nodes=1:ppn=1 -N modelPara$i -d.; done



library(R.matlab)
library(DEoptim)
library(Matrix)
source('~/Cell_Division/size.class.model_functions.R', chdir = TRUE); mat <- readMat("~/Cell_Division/day733320data.mat")

#source('/Volumes/ribalet/Cell_Division/size.class.model_functions.R', chdir = TRUE);mat <- readMat("~/Documents/DATA/SeaFlow/Cell_Division/Matlab/day733320data.mat")


# args <- commandArgs(TRUE)
# gmax <- as.numeric(args[1])
# a <- as.numeric(args[2])
# b <- as.numeric(args[3])
# E_star <- as.numeric(args[4])
# dmax <- as.numeric(args[5])


	### Matlab
		# gmax <- 0.14
		# a <-1.23
		# b <- 3.77
		# E_star <- 124
		# dmax <- 0.03


volbins <- mat$volbins
Edata <- mat$Edata
Vhists <- mat$Vhists
N_dist <- mat$N.dist
		
		timepoint <- which(apply(N_dist,2,function(x)all(!is.na(x)))) #determine the number of time points with data
	
	resol <- 60 # number of minutes per 1h-interval
	breaks <- 24*60/resol+1
		t <- 0
		dt <- 1/(resol/10)	
		ti <- seq(min(Edata[,1],na.rm=T),max(Edata[,1],na.rm=T), length.out=breaks/dt)
		ep <- data.frame(spline(Edata[,1], Edata[,2], xout=ti)) #interpolate E data according to dt resolution
		Einterp <- ep$y
		Einterp[Einterp < 0] <- 0


d.mu_N <- 3
while(d.mu_N > 2){

    gmax <- runif(1, 10^-6, 1)
    a <- runif(1, 10^-6, 15)
    b <- runif(1, 10^-6, 15)
    E_star <- runif(1, 10^-6, 500)
    dmax <- runif(1, 10^-6, 1)


	y <- (1-exp(-(Einterp/E_star)))*gmax #;y[which(Einterp > (E_star * log(2)))] <- gmax/2
	del <- dmax * (a*volbins)^b / (1 + (a*volbins)^b)
        
        Vproj <- Vhists
		Nproj <- N_dist

		pre.hr <- 0
			for(hr in timepoint[-length(timepoint)]){
					B <- matrix.conct.fast(hr=hr-1, Einterp=Einterp, volbins=volbins, gmax=gmax, a=a, b=b, E_star=E_star,dmax=dmax)						
					Nproj[,hr+1] <- as.vector(round(B %*% Nproj[,hr])) # calculate numbers of individuals
					Vproj[,hr+1] <- as.vector(B %*% Vproj[,hr]) # calculate the projected size-frequency distribution
					Vproj[,hr+1] <- Vproj[,hr+1]/sum(Vproj[,hr+1]) # normalize distribution
					pre.hr <- hr
				}
				
						
			mu_N <- diff(colSums(Nproj))/(colSums(Nproj))[-ncol(Nproj)]
			d.mu_N <- sum(mu_N, na.rm=T)*(24/length(mu_N))
			print(paste("daily growth rate=",round(d.mu_N,2)))
			
			}
			
			Vproj.noise <- jitter(Vproj, factor=0.5, amount=0);	para <- Vproj.noise
			Nproj.noise <- jitter(Nproj,factor=0.5, amount=0);	para <- Nproj.noise

	    	proj <- try(determine.opt.para(Vhists=Vproj,N_dist=Nproj,Edata=Edata,volbins=volbins)); proj$GR <- d.mu_N
			save(proj, file=paste("Simulations/simulation_a",round(a,2),"_b",round(b,2),"_gmax",round(gmax,2),"_Estar",round(E_star,1),"_dmax",round(dmax,2), sep=""))
	
	
	# library(rgl);jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",	"#FF7F00", "red", "#7F0000"))

	# percentile <- cut(unlist(para), 100); plot3d(rep(1:nrow(para), 24), rep(1:ncol(para), each=nrow(para)), z=matrix(para), col=jet.colors(100)[percentile], type='l', lwd=6, xlab="size class", ylab="time", zlab="Frequency")
