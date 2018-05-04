MeanHeatMap<-function(bdhmm_fitted, bdhmm_yeast){
########################
## Heatmap for the means
########################
means_fitted = do.call("rbind",bdhmm_fitted$hmm@emission@parameters$mean, )
require(gplots)
heat = c("dark green", "dark green", "dark green", "gold", "orange", "red")
colfct = colorRampPalette(heat)
colpal = colfct(200)
if (class(bdhmm_fitted$hmm) =='bdHMM'){
colnames(means_fitted) = names(bdhmm_yeast@emission@parameters$mean[[1]])
rownames(means_fitted) = bdhmm_fitted$hmm@stateLabel
}
return(heatmap.2(means_fitted, col=colpal, trace="none", cexCol=0.9, cexRow=0.9,
cellnote=round(means_fitted,1), notecol="black", dendrogram="row",
Rowv=TRUE, Colv=FALSE, notecex=0.9,ylab="States", xlab="Tracks", density.info="none"))
}
