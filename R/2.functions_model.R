#######################
## matrix.conct.fast ##
#######################
#Construct matrix A(t) for each time step within an hour based on delta and gamma at each 10 minute time intervals then construct B(t) which is A(t)'s multiplied for the given hour
#multiply B(t)*w(t) to get the projection to next hour dist if desired

# Einterp <- sin(seq(3,9, length.out=100));Einterp[Einterp < 0] <- 0
# volbins <- 2^seq(log2(15), log2(100), by=0.07)

.matrix.conct.fast <- function(hr, Einterp, volbins, gmax, dmax, b, E_star, resol){

		########################
		## INITIAL PARAMETERS ##
		########################
		t.nodiv <- 0 # # no division during the first X hours after dawn
		dt <- resol/60
		j <- findInterval(2 * volbins[1], volbins)
		m <- length(volbins) ## dimensions of the squared matrix

		####################
		## GAMMA FUNCTION ## fraction of cells that grow into next size class between t and t + dt
		####################
		y <- (gmax/E_star) * Einterp # NEW VERSION
		y[which(Einterp >= E_star)] <- gmax


		##########################
		## Respiration FUNCTION ## fraction of cells that shrink between t and t + dt
		##########################
				# Assumptions:
				# 1) rate of respiration is a function of growth rate
				# 2) rate of respiration constant (h-1) constant over the nighttime period (polysaccharide is drawn down linearly over the nighttime period)
		resp <- 1.35294 * y + 0.05018 # from Zavrel et al. 2019 eLife. Simplificatiojn: reg <- lm(seq(0.084, 0.199, length.out=10) ~ seq(0.025, 0.11, length.out=10))

		####################
		## DELTA FUNCTION ## fraction of cells that divide between t and t + dt
		####################
		a <- 1 # in Hunter-Cevera et al. 2014
		del <- dmax * (a*volbins)^b / (1 + (a*volbins)^b) # SOSIK et al. 2003

		# NOTE: most values of volbinske need tp be < 1 for compatibility issue with the Delta function # based on HYNES et al. 2015


		del[1:(j-1)] <- 0
				if(hr <= t.nodiv){delta <- matrix(data=0, 1, m)
					}else{delta <- matrix(del, 1, m)}

		# ### PLOT GAMMA AND DELTA
		# par(mfrow=c(2,1))
		# plot(Einterp, y, type='p', col='red', lwd=4, xlab="Radiations", ylab=paste("Gamma (per",60*dt,"min)"))
		# plot(volbins, del, type='p', col='red', lwd=4, xlab="Cell volume", ylab=paste("Delta (per",60*dt,"min)"))



		################################
		## CONSTRUCTION SPARSE MATRIX ##
		################################
		stasis_ind <- seq(1,m^2,by=m+1) # Diagonal stasis (0)
		growth_ind <- seq(2,m^2,by=m+1) # Subdiagonal growth (-1)
		div_ind <- seq((((j-1)*m)+1), m^2, by=m+1) # Superdiagonal division (j-1)
		resp_ind <- seq(((m)+1), m^2, by=m+1) # superdiagonal (+1)

		for(t in 1:(1/dt)){
			A <- matrix(data=0,nrow=m, ncol=m)

			# Cell growth (subdiagonal)
			A[growth_ind] <- y[t+hr/dt]*(1-delta[1:(m-1)])*(1-resp[t+hr/dt])

			# Division (first row and superdiagonal j-1)
			A[1,1:(j-1)] <- A[1,1:(j-1)]  + 2 * delta[1:(j-1)] # Top row; Small phytoplanktoin (i=1,..., j-1) are less than twice as big as the smallest size class, and so newly divided are put in the smallest size class.
			A[div_ind] <- 2 * delta[j:m] # The cell division terms for large (i > = j) phytoplankton

			# Stasis (main diagonal)
			A[stasis_ind] <- (1-delta)*(1-y[t+hr/dt])*(1-resp[t+hr/dt]) # the hr/dt part in the indexing is because each hour is broken up into dt segments for the irradiance spline
			A[1,1] <- (1-delta[1])*(1-y[t+hr/dt])*(1-resp[t+hr/dt]) + 2 * delta[1]
			A[m,m] <- (1-delta[m])*(1-resp[t+hr/dt])

			# Respiration (superdiagonal)
			A[1,2] <- A[1,2] + resp[t+hr/dt]
			A[resp_ind] <- resp[t+hr/dt]



					if(t == 1){B <- A}else{B <- A %*% B}
			}

		return(B)
}


###############
## sigma.lsq ##
###############
# This function calculates the sum of squares of the of the differences between the hourly observations and the model given the specified parameters
# This function returns a column vector - called by "determine.opt.para" for the optimization.
		gmax <- 0.14
		# a <-1.23
		b <- 3.77
		E_star <- 124
		dmax <- 1
		params <- data.frame(cbind(gmax, dmax, b, E_star))



	.sigma.lsq <- function(params, Einterp, N.dist, V.hists, resol){

				time.interval <- median(diff(as.numeric(colnames(V.hists))))
				res <- which(diff(as.numeric(colnames(V.hists))) == time.interval)# select time that have at least 2 consecutive time points, required for comparing the projection to the next time point
				dim <- dim(N.dist)
				sigma <- matrix(NA, dim[1], dim[2]-1) # preallocate sigma
				TotN <- as.matrix(colSums(N.dist))
				volbins <- as.numeric(row.names(V.hists))
				#n <- length(volbins)

			for(hr in res){
					B <- .matrix.conct.fast(hr=hr-1, Einterp=Einterp, volbins=volbins, gmax=as.numeric(params[1]), dmax=as.numeric(params[2]), b=as.numeric(params[3]), E_star=as.numeric(params[4]), resol=resol)
					wt <- B %*% V.hists[,hr] # calculate the projected size-frequency distribution
					wt.norm <- wt/sum(wt, na.rm=T) # normalize distribution
					sigma[,hr] <- abs(N.dist[, hr+1] - round(TotN[hr+1]*wt.norm))^2 # observed value - fitted value
					# mean.size.pred <- sum(volbins*wt.norm)
					# mean.size.obs <-  sum(volbins*V.hists[, hr+1])
					# sd.pred <- sqrt(sum((volbins - mean.size.pred)^2)/n)
					# sd.obs <- sqrt(sum((volbins - mean.size.obs)^2)/n)
					# sigma[,hr] <- abs(mean.size.obs - mean.size.pred)/sqrt(sd.obs/TotN[hr+1] + sd.pred/TotN[hr+1])
					}
			sigma <- colSums(sigma)/sum(N.dist) # sum of least squared deviations
			sigma <- sum(sigma, na.rm=T)
			return(sigma)

}



########################
## determine.opt.para ##
########################





.determine.opt.para <- function(V.hists,N.dist,Edata,resol){

		require(DEoptim)

		dt <- resol/60
		time.interval <- median(diff(as.numeric(colnames(V.hists))))
		# dt <- 1/6; breaks <- 25 ## MATLAB
		ti <- as.numeric(colnames(V.hists))

		# create Light data with 'dt' time interval.
			seq <- NULL
			for(i in 1:(length(ti)-1)){
				s <- seq(ti[i], ti[i+1], length.out=1/dt)
				seq <- c(seq, s)
			}

		ep <- data.frame(spline(Edata[,1], Edata[,2], xout=seq)) #interpolate E data according to dt resolution
		Einterp <- ep$y
		Einterp[Einterp < 0] <- 0


		##################
		## OPTIMIZATION ##
		##################
		print("Optimizing model parameters")

		f <- function(params) .sigma.lsq(params=params, Einterp=Einterp, N.dist=N.dist, V.hists=V.hists, resol=resol)

		opt <- DEoptim(f, lower=c(1e-6,1e-6,1e-6,1), upper=c(1,1,15,max(Einterp)), control=DEoptim.control(itermax=1000, reltol=1e-6, trace=10, steptol=100, strategy=2, parallelType=1))

		params <- opt$optim$bestmem
		gmax <- params[1]
		dmax <- params[2]
		b <- params[3]
		E_star <- params[4]
		resnorm <- opt$optim$bestval

		####################################################
		## Calculate projections from best fit parameters ##
		####################################################
		print(params)

		res <- which(diff(as.numeric(colnames(V.hists))) == time.interval) # select time that have at least 2 consecutive time points, required for comparing the projection to the next time point
		Vproj <- V.hists
		Nproj <- N.dist
		mu_N <- matrix(nrow=1,ncol=dim(V.hists)[2])
		volbins <- as.numeric(row.names(V.hists))

		for(hr in res){
					B <- .matrix.conct.fast(hr=hr-1, Einterp=Einterp, volbins=volbins, gmax=gmax, b=b, E_star=E_star,dmax=dmax, resol=resol)
					Nproj[,hr+1] <- round(B %*% Nproj[,hr]) # calculate numbers of individuals
					Vproj[,hr+1] <- B %*% Vproj[,hr] # calculate the projected size-frequency distribution
					Vproj[,hr+1] <- Vproj[,hr+1]/sum(Vproj[,hr+1]) # normalize distribution
					mu_N[,hr+1] <- log(sum(Nproj[,hr+1])/sum(Nproj[,hr]))/
								((as.numeric(colnames(Nproj)[hr+1])-as.numeric(colnames(Nproj)[hr]))/(time.interval))
						}
		Nproj <- colSums(Nproj)
		colnames(mu_N) <- colnames(Vproj)

		##############################
		## Growth rate calculation ##
		##############################
		d.mu_N <- 24*mean(mu_N, na.rm=T)
		print(paste("daily growth rate=",round(d.mu_N,2)))

		modelresults <- data.frame(cbind(gmax,dmax,b,E_star,resnorm), row.names=NULL)

		modelproj <- list(modelresults, mu_N, Vproj, Nproj)
		names(modelproj) <- c("modelresults", "mu_N","Vproj","Nproj")
		return(modelproj)
}
