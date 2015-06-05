# ## MODEL
# cruise <- "DeepDOM"
# model.output <- "/Volumes/ribalet/Cell_Division/"
# phyto <- 'prochloro'
# cat <- 57# number of size bin

# all.filelist <- list.files(paste(model.output,cruise,sep="/"),pattern=paste(phyto,"_modelHD_growth_",cruise,"_Ncat",cat,sep=""))
# filelist <- all.filelist[grep(pattern=paste(phyto), all.filelist)]



binning.model.output <- function(model.output, plot.raw=TRUE){

    require(plotrix)


    n <- c <- 1
    Conc.all <-  div.rate <- para.all <- Col <- NULL
    N.proj.all <- V.hist.all  <- NULL

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
            V.hist.all <- cbind(V.hist.all, v.hist.all)
            para.all <- rbind(para.all, p.all)

            col <- rep(c, nrow(dr.all))
            Col <- c(Col,col)

            leg <- unlist(list(strsplit(filelist,"_t")))[seq(2,length(filelist[1:n])*2,2)]

            if(plot.raw){
                    layout(matrix(c(1,1,2:7),4,2, byrow=T)) 
                    par(pty='m')    
                    plot(div.rate, ylab="Div Rate", xlab="time",col=Col)
                        legend("topleft",legend=leg, col=1:c, ncol=length(leg), pch=1)
                    plot(para.all[,"time"], para.all[,"gmax"], ylab="gmax", xlab="time",col = Col)
                    plot(para.all[,"time"], para.all[,"dmax"],ylab="dmax", xlab="time",col = Col)
                    plot(para.all[,"time"], para.all[,"b"],ylab="b", xlab="time",col = Col)
                    plot(para.all[,"time"], para.all[,"E_star"],ylab="E_star", xlab="time",col = Col)
                    plot(para.all[,"time"], para.all[,"resnorm"],ylab="resnorm", xlab="time",col = Col)
                    }
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
        Vproj <- V.hist.all[,order(as.numeric(colnames(V.hist.all)))]
        Para.all <- para.all[order(para.all[,"time"]),]
        
     

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
           
   
    # PARAMETERS    
        h2 <- cut(Para.all[,"time"], breaks=breaks)
        h2.time <- as.vector(tapply(Para.all[,"time"], h2, mean))#; h2.time.numc <- na.approx(h2.time.numc, na.rm=F)
        h2.gmax.mean <- as.vector(tapply(Para.all[,"gmax"], h2, mean))#; h2.gmax.mean <- na.approx(h2.gmax.mean, na.rm=F)
        h2.dmax.mean <- as.vector(tapply(Para.all[,"dmax"], h2, mean))#; h2.dmax.mean <- na.approx(h2.dmax.mean, na.rm=F)
        h2.a.mean <- as.vector(tapply(Para.all[,"a"], h2, mean))#; h2.a.mean <- na.approx(h2.a.mean, na.rm=F)
        h2.b.mean <- as.vector(tapply(Para.all[,"b"], h2, mean))#; h2.b.mean <- na.approx(h2.b.mean, na.rm=F)
        h2.E_star.mean <- as.vector(tapply(Para.all[,"E_star"], h2, mean))#; h2.E_star.mean <- na.approx(h2.E_star.mean, na.rm=F)
        h2.resnorm.mean <- as.vector(tapply(Para.all[,"resnorm"], h2, mean))#; h2.resnorm.mean <- na.approx(h2.resnorm.mean, na.rm=F)
        h2.gmax.sd <- as.vector(tapply(Para.all[,"gmax"], h2,sd))#; h2.gmax.sd <- na.approx(h2.gmax.sd, na.rm=F)
        h2.dmax.sd <- as.vector(tapply(Para.all[,"dmax"], h2, sd))#; h2.dmax.sd<- na.approx(h2.dmax.sd, na.rm=F)
        h2.a.sd <- as.vector(tapply(Para.all[,"a"], h2, sd))#; h2.a.sd <- na.approx(h2.a.sd, na.rm=F)
        h2.b.sd <- as.vector(tapply(Para.all[,"b"], h2, sd))#; h2.b.sd <- na.approx(h2.b.sd, na.rm=F)
        h2.E_star.sd <- as.vector(tapply(Para.all[,"E_star"], h2, sd))#; h2.E_star.sd <- na.approx(h2.E_star.sd, na.rm=F)
        h2.resnorm.sd <- as.vector(tapply(Para.all[,"resnorm"], h2, sd))#; h2.resnorm.sd <- na.approx(h2.resnorm.sd, na.rm=F)
        h2.time <- as.POSIXct(h2.time,origin="1970-01-01",tz='GMT')
     
         P <- data.frame(h2.time,h2.resnorm.mean,h2.resnorm.sd,h2.E_star.mean,h2.E_star.sd, h2.gmax.mean,h2.gmax.sd, h2.dmax.mean,h2.dmax.sd, h2.a.mean,h2.a.sd, h2.b.mean,h2.b.sd )

        DP <- merge(D, P, by.x=c("h.time"), by.y= c("h2.time"),all=T)

        binned.model.output <- matrix(cbind(as.array(Nproj, Vproj), DP,), nrow=3,ncol=1)

        return(binned.model.output )

        if(!(plot.raw)){

            plot.






        }








        return(DP)

}














plot.model.parameters <- function(binned.model.output){

        volbins <- unique(as.numeric(row.names(Vproj)))

        par(mfrow=c(3,2))
            for(p in c("h2.resnorm","h2.E_star","h2.gmax","h2.dmax","h2.a","h2.b")){
                plotCI(P$h2.time, P[,paste(p,'.mean',sep="")], P[,paste(p,'.sd',sep="")], col=NA, ylab=NA, main=paste(p))
                abline(v=night$UNIXtime,col='lightgrey')
                plotCI(P$h2.time, P[,paste(p,'.mean',sep="")], P[,paste(p,'.sd',sep="")],add=T)
            }
    

        del <- matrix(nrow=length(h2.time), ncol=cat)
        for(i in 1:cat){
            del[,i] <- h2.dmax.mean * (h2.a.mean*volbins[i])^h2.b.mean/ (1 + (h2.a.mean*volbins[i])^h2.b.mean)

            }
        
      par(mfrow=c(2,1),mar=c(4,4,4,4), las=1)
        plot(volbins, del[1,], ylim=c(0,0.6), type='l', col="#00007F", lwd=2, xlab="Cell volume", ylab=paste("Delta (per",10,"min)"))
                for(i in 2:nrow(del))   points(volbins, del[i,], type='l', col=jet.colors(nrow(del))[cut(as.numeric(h2.time),nrow(del))][i], lwd=2)
            ylim <- par('usr')[c(3,4)]
            xlim <- par('usr')[c(1,2)]
            color.legend(xlim[2]- diff(xlim)/40 , ylim[1], xlim[2], ylim[2], legend=format(as.POSIXct(pretty(h2.time),origin="1970-01-01"),"%d %b"), rect.col=jet.colors(100), gradient='y',align='rb')


        plot(seq(0,1000,by=10),h2.gmax.mean[1]*(1-exp(-seq(0,1000,by=10)/h2.E_star.mean[1])), ylim=c(0,0.3),type='l', col="#00007F", lwd=2, xlab="Light Intensity", ylab=paste("Gamma (per",10,"min)"))
                for(i in 1:length(h2.time)) points(seq(0,1000,by=10),h2.gmax.mean[i]*(1-exp(-seq(0,1000,by=10)/h2.E_star.mean[i])),type='l',col=jet.colors(nrow(del))[cut(as.numeric(h2.time),length(h2.time))][i],lwd=2)
                    ylim <- par('usr')[c(3,4)]
                    xlim <- par('usr')[c(1,2)]
            color.legend(xlim[2]- diff(xlim)/40 , ylim[1], xlim[2], ylim[2], legend=format(as.POSIXct(pretty(h2.time),origin="1970-01-01"),"%d %b"), rect.col=jet.colors(100), gradient='y',align='rb')

    }
    head
    










