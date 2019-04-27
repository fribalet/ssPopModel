
library(R.matlab)
df <- readMat("/Users/francois/Documents/DATA/SeaFlow/Cell_Division/Matlab_V2003/day733320data.mat")
plot(df$volbins, df$Vhists[,1], log='x')


#######################
## matrix.conct.fast ##
#######################
#Construct matrix A(t) for each time step within an hour based on delta and gamma at each 10 minute time intervals then construct B(t) which is A(t)'s multiplied for the given hour
#multiply B(t)*w(t) to get the projection to next hour dist if desired



matrix.conct.fast <- function(hr, Einterp, volbins, gmax, a, b, E_star, dmax){

		########################
		## INITIAL PARAMETERS ##
		########################
		nodiv <- 0 # # no division during the first X hours after dawn
		t.nodiv <- nodiv * (ncol(N_dist)-1)/24
		dt <- 1/(resol/10)

		# v.min <- rowSums(N_dist[,which(apply(N_dist,2,function(x)all(!is.na(x))))]) # find the size range of the population of interest
		# v.min.size <- volbins[min(which(v.min > 0))]
		v.min.size <- volbins[1]
		j <- findInterval(2 * v.min.size, volbins)
		m <- length(volbins) ## dimensions of the squared matrix

		####################
		## GAMMA FUNCTION ## fraction of cells that grow into next size class between t and t + dt
		####################
		y <- (1-exp(-(Einterp/E_star)))*gmax
		#y[which(Einterp > (E_star * log(2)))] <- gmax/2

		####################
		## DELTA FUNCTION ## fraction of cells that divide between t and t + dt
		####################
		del <- dmax * (a*volbins)^b / (1 + (a*volbins)^b)

				if(hr <= t.nodiv){delta <- matrix(data=0, 1, m)}else{delta <- matrix(del, 1,m)}


		### PLOT GAMMA AND DELTA
		# par(mfrow=c(2,1))
		# plot(Einterp, y, type='l', col='red', lwd=4, xlab="Radiations", ylab=paste("Gamma (per",60*dt,"min)"))
		# plot(volbins, del, type='l', col='red', lwd=4, xlab="Cell volume", ylab=paste("Delta (per",60*dt,"min)"))

		################################
		## CONSTRUCTION SPARSE MATRIX ##
		################################

		A <- matrix(data=0,nrow=m, ncol=m)

		stasis_ind <- seq(1,m^2,by=m+1) # indexing for matrix -goes down a column and then right to next column
		growth_ind <- seq(2,m^2,by=m+1) # corresponds to growth assignments
		div_ind <- seq((((j-1)*m)+1), m^2, by=m+1) #division assignments

		for(t in 1:(1/dt)){

		# Stasis (main diagonal)
		A[stasis_ind] <- (1-delta)*(1-y[t+hr/dt])	# the hr/dt part in the indexing is because each hour is broken up into dt segments for the irradiance spline
		A[m,m] <- 1-delta[,m]

		# Cell growth (subdiagonal region of the matrix)
		A[growth_ind] <- y[t+hr/dt]*(1-delta[1,1:m-1])

		# Division (first row and superdiagonal j-1)
		A[1,1:j-1] <- A[1,1:j-1] + 2 * delta[1:j-1]

		A[div_ind] <- 2 * delta[j:m]

		if(t == 1){B <- A}else{B <- A %*% B}

			}


		# for(t in 1:(1/dt)){

		# # Stasis (main diagonal)
		# i.stasis_ind <- 1:m; j.stasis_ind <- 1:m # indexing for matrix -goes down a column and then right to next column
		# stasis.diag.1 <- (1-y[t+hr/dt])*(1-delta[1])+2*delta[1]
		# stasis.diag.2.m1 <- (1-y[t+hr/dt])*(1-delta[2:(m-1)])
		# stasis.diag.m <- 1 - delta[m]
		# fraction.stasis <- c(stasis.diag.1, stasis.diag.2.m1, stasis.diag.m)

		# # Cell growth (subdiagonal region of the matrix)
		# i.growth_ind <- 2:m; j.growth_ind <- 1:(m-1) # corresponds to growth assignments

		# fraction.growth <- y[t+hr/dt]*(1-delta[1:(m-1)])

		# # Division (first row and superdiagonal j-1)
		# i.div_ind <- c(rep(1,j-1),2:(length(j:m))); j.div_ind <- 2:m # division assignments
		# growth.top.row <- 2 * delta[2:(j-1)]
		# growth.diag <- 2 * delta[j:m]
		# fraction.division <- c(growth.top.row, growth.diag)

		# A <- sparseMatrix(i=c(i.stasis_ind ,i.growth_ind ,i.div_ind),
							# j=c(j.stasis_ind ,j.growth_ind ,j.div_ind),
							# x=c(fraction.stasis,fraction.growth,fraction.division))
					# # the hr/dt part in the indexing is because each hour is broken up into dt segments for the irradiance spline

		# if(t == 1){B <- A}
			# else{B <- A %*% B}

		# }



		return(B)
}


###############
## sigma.lsq ##
###############
# This function calculates the sum of squares of the of the differences between the hourly observations and the model given the specified parameters
# This function returns a column vector - called by "determine.opt.para" for the optimization.
		# gmax <- 0.14
		# a <-1.23
		# b <- 3.77
		# E_star <- 124
		# dmax <- 0.03
		# params <- data.frame(cbind(gmax, a, b, E_star, dmax))
# sigma.lsq(params, Einterp, N_dist, Vhists, TotN, volbins)



	sigma.lsq <- function(params, Einterp, N_dist, Vhists, TotN, volbins){

				dim <- dim(N_dist)
				sigma <- matrix(NA, dim[1], dim[2]-1) # preallocate sigma

				timepoint <- which(apply(N_dist,2,function(x)all(!is.na(x)))) #determine the number of time points with data

			for(hr in timepoint[-length(timepoint)]){
					B <- matrix.conct.fast(hr=hr-1, Einterp=Einterp, volbins=volbins, gmax=as.numeric(params[1]), a=as.numeric(params[2]), b=as.numeric(params[3]), E_star=as.numeric(params[4]),dmax=as.numeric(params[5]))
					wt1 <- B %*% Vhists[,hr] # calculate the projected size-frequency distribution
					wt1.norm <- wt1/sum(wt1) # normalize distribution
					sigma[,hr] <- (N_dist[, hr+1] - TotN[hr+1]*wt1.norm) #observed value - fitted value
					}
			sigma <- sigma[,which(apply(sigma,2,function(x)all(!is.na(x))))]
			return(sum(sigma^2))
	}




########################
## determine.opt.para ##
########################





determine.opt.para <- function(Vhists,N_dist,Edata,volbins){

		dt <- 1/(resol/10)

		# dt <- 1/6; breaks <- 25 ## MATLAB
		TotN <- matrix(colSums(N_dist), ncol=breaks)

		ti <- seq(min(Edata[,1],na.rm=T),max(Edata[,1],na.rm=T), length.out=breaks/dt)
		ep <- data.frame(spline(Edata[,1], Edata[,2], xout=ti)) #interpolate E data according to dt resolution
		Einterp <- ep$y
		Einterp[Einterp < 0] <- 0

		##################
		## OPTIMIZATION ##
		##################


		print("Optimizing model parameters")

		f <- function(params) sigma.lsq(params, Einterp, N_dist, Vhists, TotN, volbins)

		opt <- DEoptim(f, lower=c(1e-6,1e-6,1e-6,0,1e-6), upper=c(1,15,15,500,1), control=DEoptim.control(itermax=1000, reltol=1e-6, steptol=100, trace=10))

		params <- opt$optim$bestmem
		gmax <- params[1]
		a <- params[2]
		b <- params[3]
		E_star <- params[4]
		dmax <- params[5]
		resnorm <- opt$optim$bestval

		####################################################
		## Calculate projections from best fit parameters ##
		####################################################
		print(params)

		Vproj <- Vhists
		Nproj <- N_dist

		timepoint <- which(apply(N_dist,2,function(x)all(!is.na(x)))) #determine the number of time points with data

			for(hr in timepoint[-length(timepoint)]){
					B <- matrix.conct.fast(hr=hr-1, Einterp=Einterp, volbins=volbins, gmax=gmax, a=a, b=b, E_star=E_star,dmax=dmax)
					Nproj[,hr+1] <- round(B %*% Nproj[,hr]) # calculate numbers of individuals
					Vproj[,hr+1] <- B %*% Vproj[,hr] # calculate the projected size-frequency distribution
					Vproj[,hr+1] <- Vproj[,hr+1]/sum(Vproj[,hr+1]) # normalize distribution
				}


		##############################
		## Growth rate calculations ##
		##############################

		mu_N <- diff(colSums(Nproj))/(colSums(Nproj))[-ncol(Nproj)]
		d.mu_N <- sum(mu_N, na.rm=T)*(24/length(mu_N))

		modelresults <- data.frame(cbind(gmax,a,b,E_star,dmax,resnorm), row.names=NULL)

		print(paste("daily growth rate=",round(d.mu_N,2)))
		modelproj <- list(modelresults, mu_N, Vproj, Nproj)
		names(modelproj) <- c("modelresults", "mu_N","Vproj","Nproj")
		return(modelproj)
}
