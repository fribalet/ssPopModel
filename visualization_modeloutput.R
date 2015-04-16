## MODEL
cruise <- "Thompson_9"
location.model <- "/Volumes/ribalet/Cell_Division/"
phyto <- 'prochloro'
cat <- 64 # number of size bin


all.filelist <- list.files(paste(location.model,cruise,sep="/"),pattern=paste(phyto,"_modelHD_growth_",cruise,"_Ncat",cat,sep=""))
filelist <- all.filelist[grep(pattern=paste(phyto), all.filelist)]

n <- c <- 1
Conc.all <- N.proj.all <- V.hist.all <- div.rate <- para.all <- Col <- NULL

for(file in filelist){
    #file <- filelist[11]
    load(paste(location.model ,cruise,file, sep="/"))
    print(file)
    print(n)
        dim <- conc.proj.all <- n.proj.all <- v.hist.all <- dr.all <- p.all <- NULL
            for(i in seq(2,as.numeric(n.day),by=1)){
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
        plot(para.all[,"time"], para.all[,"a"],ylab="a", xlab="time",col = Col)
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
    percentile <- cut(unlist(para), 100); plot3d(rep(1:dim(para)[1], dim(para)[2]), rep(1:dim(para)[2], each=dim(para)[1]), z=matrix(para), col=jet.colors(100)[percentile], type='l', lwd=3, xlab="size class", ylab="time", zlab="Frequency")

