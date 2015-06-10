# cruise <- "DeepDOM"
# model.output <- "/Volumes/ribalet/Cell_Division/"
# phyto <- 'prochloro'
# cat <- 57# number of size bin
# filelist <- list.files(paste(model.output,cruise,sep="/"),pattern=paste(phyto,"_modelHD_growth_",cruise,"_Ncat",cat,sep=""), full.names=T)

merge.model.output <- function(output.files, plot.raw=TRUE){

    require(plotrix)

    Col <- NULL
    c <- 1

    for(file in filelist){
        #file <- filelist[1]
        load(file)
        print(file)

            n.proj.all <- v.hist.all  <- dr.all <- p.all <- NULL

             for(i in seq(2,dim(model)[2],by=1)){
                    if(i == 2){
                        v.hist.all <- model[3,i][[1]]   
                        n.proj.all <- model[4,i][[1]]
                          para <- model[1,i][[1]]
                        p.all <- cbind(time=as.numeric(colnames(n.proj.all)), para)
                            dr <- model[2,i][[1]]
                        dr.all <- cbind(as.numeric(colnames(dr)), as.numeric(dr))


                    }else{
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
             }   
            
            if(c ==1){
                div.rate <- dr.all
                N.proj.all <- n.proj.all
                V.hist.all <- v.hist.all
                para.all <- p.all
                  
                }else{
                div.rate <- rbind(div.rate, dr.all)
                N.proj.all <- cbind.data.frame(N.proj.all, n.proj.all)
                V.hist.all <- cbind(V.hist.all, v.hist.all)
                para.all <- rbind(para.all, p.all)
                }

            col <- rep(c, nrow(dr.all))
            Col <- c(Col,col)

            if(plot.raw){
                    par(pty='m', mfrow=c(4,2))    
                    plot(as.POSIXct(div.rate[,1], origin='1970-01-01', tz='GMT'),div.rate[,2], ylab="Div Rate", xlab=NA,col=Col)
                    plot(1,1, pch=NA, ylab=NA, xlab=NA, bty='n', xaxt='n', yaxt='n')
                    legend("topleft",legend=filelist, col=1:c, ncol=2, pch=1, bty='n')
                    plot(as.POSIXct(para.all$time, origin='1970-01-01', tz='GMT'), para.all[,"gmax"], ylab="gmax", xlab=NA,col = Col)
                    plot(as.POSIXct(para.all$time, origin='1970-01-01', tz='GMT'), para.all[,"dmax"],ylab="dmax", xlab=NA,col = Col)
                    plot(as.POSIXct(para.all$time, origin='1970-01-01', tz='GMT'), para.all[,"b"],ylab="b", xlab=NA,col = Col)
                    plot(as.POSIXct(para.all$time, origin='1970-01-01', tz='GMT'), para.all[,"E_star"],ylab="E_star", xlab=NA,col = Col)
                    plot(as.POSIXct(para.all$time, origin='1970-01-01', tz='GMT'), para.all[,"resnorm"],ylab="resnorm", xlab=NA,col = Col)
                    }
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
        print('merging model outputs...')

        time.interval <- median(diff(as.numeric(colnames(v.proj))))
        time.range <- range(as.numeric(colnames(Vproj)))
        breaks <- seq(time.range[1],time.range[2], by= time.interval)

    
    # DIVISION 
        h <- cut(Div.rate[,1], breaks=breaks, labels=F)
        h.time <- round(as.vector(tapply(Div.rate[,1], h, function(x) mean(x, na.rm=T))))
        h.dr.mean <- as.vector(tapply(Div.rate[,2], h, function(x) mean(x, na.rm=T)))
        h.dr.sd <- as.vector(tapply(Div.rate[,2], h, function(x) sd(x, na.rm=T)))

        D <- data.frame(cbind(h.time, h.dr.mean, h.dr.sd))
           
    # SIZE DISTRIBUTION
         Vproj.mean <- Vproj.sd <- NULL
                for(i in 1:nrow(Vproj)){
                    vhist <- cbind(as.numeric(colnames(Vproj[i,])), unlist(Vproj[i,]))
                    h1 <- cut(vhist[,1],breaks=breaks, labels=F)
                    h1.time <- round(as.vector(tapply(vhist [,1], h1, function(x) mean(x, na.rm=T))))
                    h1.vhist.mean <- as.vector(tapply(vhist[,2], h1, function(x) mean(x, na.rm=T)))
                    h1.vhist.sd <- as.vector(tapply(vhist[,2], h1, function(x) sd(x, na.rm=T)))
                    Vproj.mean <- rbind(Vproj.mean, h1.vhist.mean)
                    Vproj.sd <- rbind(Vproj.sd, h1.vhist.sd)
                    }
            row.names(Vproj.mean) <-  row.names(Vproj.sd) <- row.names(Vproj)
            colnames(Vproj.mean) <-  colnames(Vproj.sd) <- h1.time
 
         Nproj.mean <- Nproj.sd <- NULL
                 for(i in 1:nrow(Nproj)){
                    nhist <- cbind(as.numeric(colnames(Nproj[i,])), unlist(Nproj[i,]))
                    h1 <- cut(nhist[,1],breaks=breaks, labels=F)
                    h1.time <- round(as.vector(tapply(nhist [,1], h1, function(x) mean(x, na.rm=T))))
                    h1.nhist.mean <- as.vector(tapply(nhist[,2], h1, function(x) mean(x, na.rm=T)))
                    h1.nhist.sd <- as.vector(tapply(nhist[,2], h1, function(x) sd(x, na.rm=T)))
                    Nproj.mean <- rbind(Nproj.mean, h1.nhist.mean)
                    Nproj.sd <- rbind(Nproj.sd, h1.nhist.sd)
                    }
            row.names(Nproj.mean) <-  row.names(Nproj.sd) <- row.names(Nproj)
            colnames(Nproj.mean) <-  colnames(Nproj.sd) <- h1.time





    # PARAMETERS    
        h2 <- cut(Para.all[,"time"], breaks=breaks)
        h2.time <- round(as.vector(tapply(Para.all[,"time"], h2, function(x) mean(x, na.rm=T))))
        h2.gmax.mean <- as.vector(tapply(Para.all[,"gmax"], h2, function(x) mean(x, na.rm=T)))
        h2.dmax.mean <- as.vector(tapply(Para.all[,"dmax"], h2, function(x) mean(x, na.rm=T)))
        h2.b.mean <- as.vector(tapply(Para.all[,"b"], h2, function(x) mean(x, na.rm=T)))
        h2.E_star.mean <- as.vector(tapply(Para.all[,"E_star"], h2, function(x) mean(x, na.rm=T)))
        h2.resnorm.mean <- as.vector(tapply(Para.all[,"resnorm"], h2, function(x) mean(x, na.rm=T)))
        h2.gmax.sd <- as.vector(tapply(Para.all[,"gmax"], h2,function(x) sd(x, na.rm=T)))
        h2.dmax.sd <- as.vector(tapply(Para.all[,"dmax"], h2, function(x) sd(x, na.rm=T)))
        h2.b.sd <- as.vector(tapply(Para.all[,"b"], h2, function(x) sd(x, na.rm=T)))
        h2.E_star.sd <- as.vector(tapply(Para.all[,"E_star"], h2, function(x) sd(x, na.rm=T)))
        h2.resnorm.sd <- as.vector(tapply(Para.all[,"resnorm"], h2, function(x) sd(x, na.rm=T)))
        h2.time <- as.POSIXct(h2.time,origin="1970-01-01",tz='GMT')
     
         P <- data.frame(h2.time,h2.resnorm.mean,h2.resnorm.sd,h2.E_star.mean,h2.E_star.sd, 
                            h2.gmax.mean,h2.gmax.sd, h2.dmax.mean,h2.dmax.sd, h2.b.mean,h2.b.sd)

        DP <- merge(D, P, by.x=c("h.time"), by.y= c("h2.time"),all=T)



    

        if(!(plot.raw)){
            par(mfrow=c(3,2))
            DP$h.time <- as.POSIXct(DP$h.time, origin="1970-01-01", tz='GMT')
                 for(i in seq(2,12,by=2))    plotCI(DP$h.time, DP[,i], uiw=DP[,i+1], sfrac=0, xlab=NA, ylab=paste(colnames(DP)[i]))
                   }




    merged.model.output <- list(DP, Vproj.mean, Nproj.mean, Vproj.sd, Nproj.sd)
    names(merged.model.output) <- c("estimates","Vproj","Nproj","Vproj.sd", "Nproj.sd")
    
    return(merged.model.output)


}





plot.parameters <- function(merged.model.output){

    jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

        volbins <- unique(as.numeric(row.names(merged.model.output$Vproj)))
        para <- binned.model.output$estimates
        h2.time <- para$h.time
       del <- matrix(nrow=length(h2.time), ncol=cat)
        
        for(i in 1:cat)  del[,i] <- para$h2.dmax.mean * (volbins[i]/max(volbins))^para$h2.b.mean/ (1 + (volbins[i]/max(volbins))^para$h2.b.mean)

            
        
      par(mfrow=c(2,1),mar=c(4,4,4,4), las=1)
        plot(volbins, del[1,], ylim=c(0,max(del, na.rm=T)), type='l', col="#00007F", lwd=2, xlab="Cell volume", ylab=paste("Delta (per",10,"min)"))
                for(i in 2:nrow(del))   points(volbins, del[i,], type='l', col=jet.colors(nrow(del))[cut(as.numeric(h2.time),nrow(del))][i], lwd=2)
            ylim <- par('usr')[c(3,4)]
            xlim <- par('usr')[c(1,2)]
            color.legend(xlim[2]- diff(xlim)/40 , ylim[1], xlim[2], ylim[2], legend=format(as.POSIXct(pretty(h2.time),origin="1970-01-01"),"%d %b"), rect.col=jet.colors(100), gradient='y',align='rb')

         max <- max(para$h2.gmax.mean*(1-exp(-1000)/para$h2.E_star.mean), na.rm=T)
        plot(seq(0,1000,by=10),para$h2.gmax.mean[1]*(1-exp(-seq(0,1000,by=10)/para$h2.E_star.mean[1])), ylim=c(0,max),type='l', col="#00007F", lwd=2, xlab="Light Intensity", ylab=paste("Gamma (per",10,"min)"))
                for(i in 1:length(h2.time)) points(seq(0,1000,by=10),para$h2.gmax.mean[i]*(1-exp(-seq(0,1000,by=10)/para$h2.E_star.mean[i])),type='l',col=jet.colors(nrow(del))[cut(as.numeric(h2.time),length(h2.time))][i],lwd=2)
                    ylim <- par('usr')[c(3,4)]
                    xlim <- par('usr')[c(1,2)]
            color.legend(xlim[2]- diff(xlim)/40 , ylim[1], xlim[2], ylim[2], legend=format(as.POSIXct(pretty(h2.time),origin="1970-01-01"),"%d %b"), rect.col=jet.colors(100), gradient='y',align='rb')

    }











