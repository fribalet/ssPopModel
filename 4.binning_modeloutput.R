## MODEL
cruise <- "DeepDOM"
model.output <- "/Volumes/ribalet/Cell_Division/"
phyto <- 'prochloro'
cat <- 57# number of size bin

## PAR
Par.path <- paste(model.output ,"/",cruise,"/Par_",cruise,sep="")
Par <- read.csv(Par.path, sep=",")
Par$time <- as.POSIXct(Par$time, tz='GMT')
Par$UNIXtime <- as.numeric(Par$time)
night <- Par[Par[,"par"] < 2,] ## select Dusk and Dawn time for Day i




all.filelist <- list.files(paste(model.output,cruise,sep="/"),pattern=paste(phyto,"_modelHD_growth_",cruise,"_Ncat",cat,sep=""))
filelist <- all.filelist[grep(pattern=paste(phyto), all.filelist)]

n <- c <- 1
Conc.all <-  div.rate <- para.all <- Col <- NULL
N.proj.all <- V.hist.all  <- matrix(nrow=cat)

for(file in filelist){
    #file <- filelist[11]
    load(paste(model.output,cruise,file, sep="/"))
    print(file)
    print(n)
        n.proj.all <- v.hist.all  <- matrix(nrow=cat)
        conc.proj.all <- dr.all <- p.all <- NULL

            for(i in seq(2,dim(model)[2],by=1)){
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

        col <- rep(c, nrow(dr.all))
        Col <- c(Col,col)

        leg <- unlist(list(strsplit(filelist,"_t")))[seq(2,length(filelist[1:n])*2,2)]

        layout(matrix(c(1,1,2:7),4,2, byrow=T)) 
        par(pty='m')    
        plot(div.rate, ylab="Div Rate", xlab="time",col=Col)
            #abline(v=night$UNIXtime,col='lightgrey');points(div.rate,col=Col)
            legend("topleft",legend=leg, col=1:c, ncol=length(leg), pch=1)
        plot(para.all[,"time"], para.all[,"gmax"], ylab="gmax", xlab="time",col = Col)
        plot(para.all[,"time"], para.all[,"dmax"],ylab="dmax", xlab="time",col = Col)
        plot(para.all[,"time"], para.all[,"b"],ylab="b", xlab="time",col = Col)
        plot(para.all[,"time"], para.all[,"E_star"],ylab="E_star", xlab="time",col = Col)
        plot(para.all[,"time"], para.all[,"resnorm"],ylab="resnorm", xlab="time",col = Col)

        # # 
        # names(para) <- c("gmax","a","b","E_star","dmax","resnorm")
        # par(mfrow=c(4,2))
        # barplot(d.GR, col='grey', main="GR")
        # for(i in 1:6) barplot(para[,i], main=colnames(para)[i])
n <- n + 1  
c <- c + 1
}

### ORDERING DATA by TIME

    Div.rate <- div.rate[order(div.rate[,1]),]
    Nproj <- N.proj.all[,order(as.numeric(colnames(N.proj.all)))]
    Conc.proj <- Conc.all[order(Conc.all[,1]),]
    Vproj <- V.hist.all[,order(as.numeric(colnames(V.hist.all)))]
    Para.all <- para.all[order(para.all[,"time"]),]
        
  ### VISUALIZATION      
    para <- Vproj
    percentile <- cut(unlist(para), 100); plot3d(rep(1:dim(para)[1], dim(para)[2]), rep(1:dim(para)[2], each=dim(para)[1]), z=as.matrix(para), col=jet.colors(100)[percentile], type='l', lwd=3, xlab="size class", ylab="time", zlab="Frequency")




    

###############
### BINNING ###
###############

        breaks <- seq(as.numeric(colnames(Vproj)[1]),as.numeric(colnames(Vproj)[dim(Vproj)[2]]),by=60*60)

    
    # DIVISION 
        h <- cut(Div.rate[,1], breaks=breaks, labels=F)
        h.time.numc <- as.vector(tapply(Div.rate[,1], h, function(x) mean(x, na.rm=T)))#; h.time.numc <- na.approx(h.time.numc, na.rm=F)
        h.dr.mean <- as.vector(tapply(Div.rate[,2], h, function(x) mean(x, na.rm=T)))#; h.dr.mean <- na.approx(h.dr.mean, na.rm=F)
        h.dr.sd <- as.vector(tapply(Div.rate[,2], h, function(x) sd(x, na.rm=T)))#; h.dr.sd <- na.approx(h.dr.sd, na.rm=F)
        #h.time <- as.POSIXct(h.time.numc,origin="1970-01-01",tz='GMT')
        h.time <- as.POSIXct(breaks[findInterval(h.time.numc, breaks)],origin="1970-01-01",tz='GMT')

        D <- data.frame(cbind(h.time, h.dr.mean, h.dr.sd))
           
            par(mfrow=c(3,1))
            plot(h.time, h.dr.mean, ylim=c(0,max(h.dr.mean, na.rm=T)*1.3))      
            abline(v=night$UNIXtime,col='lightgrey')
            plotCI(h.time, h.dr.mean, h.dr.sd, add=T)       


   # PAR
        h4 <- cut(as.numeric(Par$time), breaks=breaks, labels=F)
        h4.time.numc <- as.vector(tapply(as.numeric(Par$time), h4, mean))#; h4.time.numc <- na.approx(h4.time.numc, na.rm=F)
        h4.par.mean <- as.vector(tapply(Par[,"par"], h4, mean))#; h4.par.mean <- na.approx(h4.par.mean, na.rm=F)
        h4.par.sd <- as.vector(tapply(Par[,"par"], h4, sd))#; h4.par.sd <- na.approx(h4.par.sd, na.rm=F)
        #h4.time <- as.POSIXct(h4.time.numc,origin="1970-01-01",tz='GMT')
        h4.time <- as.POSIXct(breaks[findInterval(h4.time.numc, breaks)],origin="1970-01-01",tz='GMT')
 
         L <- data.frame(cbind(h4.time, h4.par.mean, h4.par.sd))

            plotCI(h4.time, h4.par.mean, h4.par.sd)     
    
            DL <- merge(D, L, by.x=c("h.time"), by.y= c("h4.time"),all=T)
        
        write.csv(D, paste(model.output,"/",cruise,"/",phyto,"_HD_",cruise, ".division_rate.csv",sep=""),quote=F, row.names=F)





    # PARAMETERS    
        h5 <- cut(Para.all[,"time"], breaks=breaks)
        h5.time <- as.vector(tapply(Para.all[,"time"], h5, mean))#; h5.time.numc <- na.approx(h5.time.numc, na.rm=F)
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
        h5.time <- as.POSIXct(h5.time,origin="1970-01-01",tz='GMT')
    
        id <- findInterval(h5.time, pop$time, rightmost.closed=F)
        h5.lat <- pop[id,"lat"]
        h5.lon <- pop[id,"long"]
    
        h5.gamma.mean <- h5.gmax.mean*(1-exp(-h4.par.mean/h5.E_star.mean))
        h5.gamma.sd <- h5.gmax.sd*(1-exp(-h4.par.sd/h5.E_star.sd))
    
        P <- data.frame(h5.time,h5.resnorm.mean,h5.resnorm.sd,h5.E_star.mean,h5.E_star.sd, h5.gmax.mean,h5.gmax.sd, h5.dmax.mean,h5.dmax.sd, h5.a.mean,h5.a.sd, h5.b.mean,h5.b.sd )


        volbins <- unique(as.numeric(row.names(Vproj)))

        par(mfrow=c(3,2))
            for(p in c("h5.resnorm","h5.E_star","h5.gmax","h5.dmax","h5.a","h5.b")){
                plotCI(P$h5.time, P[,paste(p,'.mean',sep="")], P[,paste(p,'.sd',sep="")], col=NA, ylab=NA, main=paste(p))
                abline(v=night$UNIXtime,col='lightgrey')
                plotCI(P$h5.time, P[,paste(p,'.mean',sep="")], P[,paste(p,'.sd',sep="")],add=T)
            }
    

        del <- matrix(nrow=length(h5.time), ncol=cat)
        for(i in 1:cat){
            del[,i] <- h5.dmax.mean * (h5.a.mean*volbins[i])^h5.b.mean/ (1 + (h5.a.mean*volbins[i])^h5.b.mean)

            }
        
      par(mfrow=c(2,1),mar=c(4,4,4,4), las=1)
        plot(volbins, del[1,], ylim=c(0,0.6), type='l', col="#00007F", lwd=2, xlab="Cell volume", ylab=paste("Delta (per",10,"min)"))
                for(i in 2:nrow(del))   points(volbins, del[i,], type='l', col=jet.colors(nrow(del))[cut(as.numeric(h5.time),nrow(del))][i], lwd=2)
            ylim <- par('usr')[c(3,4)]
            xlim <- par('usr')[c(1,2)]
            color.legend(xlim[2]- diff(xlim)/40 , ylim[1], xlim[2], ylim[2], legend=format(as.POSIXct(pretty(h5.time),origin="1970-01-01"),"%d %b"), rect.col=jet.colors(100), gradient='y',align='rb')


        plot(seq(0,1000,by=10),h5.gmax.mean[1]*(1-exp(-seq(0,1000,by=10)/h5.E_star.mean[1])), ylim=c(0,0.3),type='l', col="#00007F", lwd=2, xlab="Light Intensity", ylab=paste("Gamma (per",10,"min)"))
                for(i in 1:length(h5.time)) points(seq(0,1000,by=10),h5.gmax.mean[i]*(1-exp(-seq(0,1000,by=10)/h5.E_star.mean[i])),type='l',col=jet.colors(nrow(del))[cut(as.numeric(h5.time),length(h5.time))][i],lwd=2)
                    ylim <- par('usr')[c(3,4)]
                    xlim <- par('usr')[c(1,2)]
            color.legend(xlim[2]- diff(xlim)/40 , ylim[1], xlim[2], ylim[2], legend=format(as.POSIXct(pretty(h5.time),origin="1970-01-01"),"%d %b"), rect.col=jet.colors(100), gradient='y',align='rb')

    
    


        write.csv(P, paste(model.output,"/",cruise,"/",phyto,"_HD_",cruise, ".parameters.csv",sep=""),quote=F, row.names=F)










