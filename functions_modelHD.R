#######################
## matrix.conct.fast ##
#######################
#Construct matrix A(t) for each time step within an hour based on delta and gamma at each 10 minute time intervals then construct B(t) which is A(t)'s multiplied for the given hour
#multiply B(t)*w(t) to get the projection to next hour dist if desired



matrix.conct.fast <- function(hr, Einterp, volbins, gmax, dmax, a, b, E_star){
	
		########################
		## INITIAL PARAMETERS ##
		########################
		nodiv <- 0 # # no division during the first X hours after dawn 
		t.nodiv <- nodiv * (ncol(N_dist)-1)/24 
		dt <- (resol/10)^-1 
		j <- findInterval(2 * volbins[1], volbins)
		m <- length(volbins) ## dimensions of the squared matrix

		####################		
		## GAMMA FUNCTION ## fraction of cells that grow into next size class between t and t + dt
		#################### 
		# y <- gmax*(1-exp(-Einterp)/E_star)
		y <- gmax*(1-exp(-Einterp/(E_star*gmax)))

		# y[which(Einterp > (E_star * log(2)))] <- gmax/2
		
		####################		
		## DELTA FUNCTION ## fraction of cells that divide between t and t + dt
		#################### 
		del <- dmax * a*(volbins)^b / (1 + a*(volbins)^b)
		del[1:(j-1)] <- 0		
				# if(hr <= t.nodiv){delta <- matrix(data=0, 1, m)
					# }else{delta <- matrix(del, 1, m)}
		delta <- matrix(del, 1, m)
		
		# ### PLOT GAMMA AND DELTA
		# par(mfrow=c(2,1))
		# plot(Einterp, y, type='p', col='red', lwd=4, xlab="Radiations", ylab=paste("Gamma (per",60*dt,"min)"))
		# plot(volbins, del, type='p', col='red', lwd=4, xlab="Cell volume", ylab=paste("Delta (per",60*dt,"min)"))

		################################
		## CONSTRUCTION SPARSE MATRIX ##
		################################
		stasis_ind <- seq(m+2,m^2-1,by=m+1) # Diagonal stasis (0)
		growth_ind <- seq(2,(m-1)^2,by=m+1) # Subdiagonal growth (-1)
		div_ind <- seq((((j-1)*m)+1), m^2, by=m+1) # Superdiagonal division (j-1)
		
		for(t in 1:(1/dt)){
			A <- matrix(data=0,nrow=m, ncol=m)
			
			# Cell growth (subdiagonal region of the matrix)
			A[growth_ind] <- y[t+hr/dt]*(1-delta[1:(m-2)])	

			# Division (first row and superdiagonal j-1)
			A[1,2:(j-1)] <- 2 * delta[2:(j-1)] # Top row; Small phytoplanktoin (i=1,..., j-1) are less than twice as big as the smallest size class, and so newly divided are put in the smallest size class.
			A[div_ind] <- 2 * delta[j:m] # The cell division terms for large (i > = j) phytoplankton
		
			# Stasis (main diagonal)
			A[stasis_ind] <- (1-delta[2:(m-1)])*(1-y[t+hr/dt])	# the hr/dt part in the indexing is because each hour is broken up into dt segments for the irradiance spline
			A[1,1] <- (1-delta[1])*(1-y[t+hr/dt]) + 2 * delta[1]
			A[m,m] <- 1-delta[m]


			# A[m,m] <- 1-colSums(A)[m]
			# A[1,1] <- 1-colSums(A)[1]
			# A[stasis_ind] <- 1-colSums(A)[-c(1,m)]	
		
					if(t == 1){B <- A}else{B <- A %*% B}
			}

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
		# params <- data.frame(cbind(gmax, dmax, a, b, E_star))
		# params <- as.numeric(proj$modelresults)

	

	sigma.lsq <- function(params, Einterp, N.dist, V.hists, TotN, volbins){
				dim <- dim(N.dist)
				sigma <- matrix(NA, dim[1], dim[2]-1) # preallocate sigma
			
					
			for(hr in 1:24){
					B <- matrix.conct.fast(hr=hr-1, Einterp=Einterp, volbins=volbins, gmax=as.numeric(params[1]), dmax=as.numeric(params[2]), a=as.numeric(params[3]),b=as.numeric(params[4]), E_star=as.numeric(params[5]))
						
					wt <- B %*% V.hists[,hr] # calculate the projected size-frequency distribution 
					wt.norm <- wt/sum(wt, na.rm=T) # normalize distribution
					sigma[,hr] <- (round(N.dist[, hr+1] - TotN[hr+1]*wt.norm)^2) #observed value - fitted value
					#sigma[,hr] <- abs(V.hists[, hr+1] - wt.norm) #observed value - fitted value
					}
			sigma <- colSums(sigma)/colSums(na.omit(N.dist[,-1]))
			sigma <- sum(sigma, na.rm=T)
			return(sigma)

}



########################
## determine.opt.para ##
########################




	
determine.opt.para <- function(V.hists,N.dist,Edata,volbins){
		
		dt <- 1/(resol/10)	
			
		# dt <- 1/6; breaks <- 25 ## MATLAB
		
		TotN <- matrix(colSums(na.omit(N.dist)), ncol=breaks)
		ti <- seq(min(Edata[,1],na.rm=T),max(Edata[,1],na.rm=T), length.out=breaks/dt)
		ep <- data.frame(spline(Edata[,1], Edata[,2], xout=ti)) #interpolate E data according to dt resolution
		Einterp <- ep$y
		Einterp[Einterp < 0] <- 0
		

		##################
		## OPTIMIZATION ##
		##################
		print("Optimizing model parameters")
		
		f <- function(params) sigma.lsq(params=params, Einterp=Einterp, N.dist=N.dist, V.hists=V.hists, TotN=TotN, volbins=volbins)
			
		opt <- DEoptim(f, lower=c(1e-6,1e-6,1e-6,1e-6,1), upper=c(1,1,15,15,2000), control=DEoptim.control(itermax=1000, reltol=1e-6, trace=10, steptol=100))
		
		params <- opt$optim$bestmem
		gmax <- params[1]
		dmax <- params[2]
		a <- params[3]
		b <- params[4]
		E_star <- params[5]
		resnorm <- opt$optim$bestval
										
		####################################################
		## Calculate projections from best fit parameters ##	
		####################################################
		print(params)

		Vproj <- V.hists
		Nproj <- N.dist
		mu_N <- matrix(nrow=1,ncol=dim(V.hists)[2])

			for(hr in 1:24){
					B <- matrix.conct.fast(hr=hr-1, Einterp=Einterp, volbins=volbins, gmax=gmax, a=a, b=b, E_star=E_star,dmax=dmax)
					Nproj[,hr+1] <- round(B %*% Nproj[,hr]) # calculate numbers of individuals
					Vproj[,hr+1] <- B %*% Vproj[,hr] # calculate the projected size-frequency distribution
					Vproj[,hr+1] <- Vproj[,hr+1]/sum(Vproj[,hr+1]) # normalize distribution
					mu_N[,hr+1] <- log(sum(Nproj[,hr+1])/sum(Nproj[,hr]))
				
				}

		colnames(mu_N) <- colnames(Nproj)
		
		##############################
		## Growth rate calculation ##
		##############################
		d.mu_N <- sum(mu_N, na.rm=T)
		print(paste("daily growth rate=",round(d.mu_N,2)))
			
		modelresults <- data.frame(cbind(gmax,dmax,a,b,E_star,resnorm), row.names=NULL)
		
		modelproj <- list(modelresults, mu_N, Vproj, Nproj)
		names(modelproj) <- c("modelresults", "mu_N","Vproj","Nproj")
		return(modelproj)
}

