##############################
### PLOT size distribution ###
##############################
plot.size.distribution <- function(freq.distribution, mode = c('log', 'lin'), ...){

    require(rgl)
    jet.colors <- viridis::viridis
    mode <- as.character(mode[1])

    param <- freq.distribution
    percentile <- cut(unlist(param), 100)

    # in linear scale
    if(mode =='lin'){
        plot3d(rep(as.numeric(row.names(param)), dim(param)[2]),
                rep(as.numeric(colnames(param)), each=dim(param)[1]) ,
                unlist(param),
                col=jet.colors(100)[percentile], xlab="size class", ylab="time", zlab="Frequency",...)
     }


    # in log scale
    if(mode =='log'){
        plot3d(log2(rep(as.numeric(row.names(param)), dim(param)[2])),
                rep(as.numeric(colnames(param)), each=dim(param)[1]) ,
                unlist(param),
                col=jet.colors(100)[percentile],  xlab="size class", ylab="time", zlab="Frequency",...)
    }
}
