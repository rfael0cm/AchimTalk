## requires: gplots, MASS, Gviz
## work with vignette es
### Course in bdHMM
## 1. load a bdHMM object and see its characteristics
## 2. Generate 10^4 time points from the model
## 3. Fit a standar HMM model to the generated Data
## 4. See some characteristics
## 5. Fit a bdHMM
## 6. Compare it with the HMM and bdHMM

# 1. LOAD A bdHMM OBJECT AND SEE ITS CHARACTERISTICS

library(STAN)
library(gplots)

load('/home/campos/Desktop/Achim_talk/bdhmm_fitted_class.RData')

bdhmm_fitted # This bdHMM model has been fitted to a set of 10 different tracks, 2 of them directional.

bdhmm_fitted$loglik # We can see the last 10 iteration during Baum-Welch Method

bdhmm_fitted$stateLabel # Name of the states, this model was fitted for 12 states 5 directional (F/R 1:5) and 2 undirectional (U1, U2)

bdhmm_fitted$hmm@directedObs # Directionality of the Data we have 8 tracks undirectional (TF, histone modification,...) coded by 0 and 2 directional (strand specific RNA-seq data) coded by 1

bdhmm_fitted$hmm@initProb # Initial Probabilities 

bdhmm_fitted$hmm@transMat # transition Matrix between states


source('/home/campos/Desktop/Achim_talk/MeanHeatMap_function.R')


MeanHeatMap(bdhmm_fitted, bdhmm_yeast)


heat=c('dark grey','steelblue', 'yellow','gold', 'orange', 'dark orange',rep('red',10),rep('white',200))
colfct = colorRampPalette(heat)
colpal = colfct(200)
heatmap.2(bdhmm_fitted$hmm@transMat, col=colpal,trace="none", cexCol=0.9, cexRow=0.9, notecol="black", dendrogram="none",
Rowv=F, Colv=FALSE, notecex=0.9)


## GENERATE 10^3 OBSERVATIONS 

source('/home/campos/Desktop/CreateData/Function/GenerateData_function.R')

Data_G<-GenerateData(1299, bdhmm_fitted)
colnames(Data_G$observation)<-names(bdhmm_yeast@emission@parameters$mean[[1]])
table(Data_G$viterbi)

ChrMatrix<-list()
ChrMatrix$chr<-Data_G$observation


## GENERATE HMM MODEL

nStates=7

myMat = ChrMatrix$chr[apply(ChrMatrix$chr, 1, function(x) all(! is.na(x))),]
myMat[, c("YPDexprW", "YPDexprC")] = t(apply(myMat[, c("YPDexprW", "YPDexprC")], 1,
                                             sort, decreasing=TRUE))
km = kmeans(myMat, centers=7, iter.max=1000, nstart=100)$centers

Mean_hmm <- lapply(1:nrow(km), function(x)km[x,])


Covs_hmm<-cov(myMat[complete.cases(myMat),])

Covs_hmm<-lapply(1:nStates, function(x)Covs_hmm)


gaussEmission = HMMEmission(type="Gaussian",
        parameters=list(mean=Mean_hmm, cov=new_Colapsed$hmm@emission@parameters$cov),
        nStates=nStates)

transMat = matrix(1/nStates, nrow=nStates, ncol=nStates)
initProb = rep(1/nStates, nStates)
hmm = HMM(initProb=initProb, transMat=transMat,
        emission=gaussEmission, nStates=nStates,
        status="initial")

hmm_fitted<-fitHMM(ChrMatrix, hmm)

hmm_viterbi<-getViterbi(hmm_fitted$hmm,ChrMatrix )

table(hmm_viterbi$chr, Data_G$viterbi)
names(hmm@emission@parameters$mean[[1]])<-colnames(Data_G$observation)
MeanHeatMap(hmm_fitted, hmm)

heat=c('dark grey','steelblue', 'yellow','gold', 'orange', 'dark orange',rep('red',10),rep('white',200))
colfct = colorRampPalette(heat)
colpal = colfct(200)
heatmap.2(hmm_fitted$hmm@transMat, col=colpal,trace="none", cexCol=0.9, cexRow=0.9, notecol="black", dendrogram="none",
Rowv=F, Colv=FALSE, notecex=0.9)

## GENERATE bdHMM MODEL

nStates<- 11
stateLabel<-bdhmm_fitted$hmm@stateLabel

myMat = ChrMatrix$chr[apply(ChrMatrix$chr, 1, function(x) all(! is.na(x))),]
myMat[, c("YPDexprW", "YPDexprC")] = t(apply(myMat[, c("YPDexprW", "YPDexprC")], 1,
                                             sort, decreasing=TRUE))
km = kmeans(myMat, centers=7, iter.max=1000, nstart=100)$centers
km<-rbind(km, km[c(1:4),])
km<-km[c(1:4,8:11,5,6,7),]

Mean_bdhmm <- lapply(1:nrow(km), function(x)km[x,])

for (i in 5:8){
Mean_bdhmm[[i]]<-Mean_bdhmm[[i]][c(1:8,10,9)]
names(Mean_bdhmm[[i]])<-names(Mean_bdhmm[[1]])
}



Covs_hmm<-cov(myMat[complete.cases(myMat),])

Covs_hmm<-lapply(1:nStates, function(x)Covs_hmm)



gaussEmission <- HMMEmission(type="Gaussian",
parameters=list(mean=bdhmm_fitted$hmm@emission@parameters$mean, cov=bdhmm_fitted$hmm@emission@parameters$cov),
nStates=nStates)
dirobs = as.integer(c(rep(0,8), 1, 1))

transMat <- bdhmm_fitted$hmm@transMat
initProb <-bdhmm_fitted$hmm@initProb
bdhmm_GD = bdHMM(initProb=initProb, transMat=transMat,
        emission=gaussEmission, nStates=nStates,
        status="initial", stateLabel=stateLabel,
        transitionsOptim="analytical", directedObs=dirobs)
bdhmm_fittedGD = fitHMM(ChrMatrix, bdhmm_yeast, maxIters=100, verbose=FALSE)
viterbiGD = getViterbi(bdhmm_fittedGD$hmm, ChrMatrix)


table(viterbiGD$chr,Data_G$viterbi )

MeanHeatMap(bdhmm_fittedGD, bdhmm_GD)

## PLOT WITH TRACKS

library(Gviz)
ucscChromosomeNames=FALSE
gtrack<-GenomeAxisTrack()

names(faccols) = colnames(Data_G$observation)
chr = "chrIV"
gen = "sacCer3"
gtrack <- GenomeAxisTrack()

faccols = hcl(h = seq(15, 375 - 360/dim(ChrMatrix$chr)[2],
length = dim(ChrMatrix$chr)[2])%%360, c = 100, l = 65)
names(faccols) = colnames(ChrMatrix$chr)
dlist=list()
for(n in colnames(ChrMatrix$chr)) {
dlist[[n]] = DataTrack(data = ChrMatrix$chr[,n],
start = yeastTF_probeAnno_ex$chr04,
end = yeastTF_probeAnno_ex$chr04+8,
chromosome = "chrIV", genome=gen,
name = n, type="h", col=faccols[n])
}



library(GenomicRanges)
library(IRanges)
myViterbiDirs = list(F=c("F1", "F2", "F3", "F4"), U=c("U1", "U2","U3"),
R=c("R1", "R2", "R3", "R4"))
myViterbiPanels = list()
cols = rainbow(7)
cols = cols[c(1:5,1:5,6:7)]
names(cols) = stateLabel

myHiddenStates = list()

for(dir in c("F", "U", "R")) {
myPos = yeastTF_probeAnno_ex$chr04 >= 1217060 & yeastTF_probeAnno_ex$chr04 <= 1225000
myRle = Rle(viterbiGD$chr[myPos])
currItems = which(myRle@values %in% myViterbiDirs[[dir]])

start = yeastTF_probeAnno_ex$chr04[myPos][start(myRle)][currItems]
width = myRle@lengths[currItems]
ids = as.character(myRle@values[currItems])
values = as.character(myRle@values[currItems])
myViterbiPanels[[dir]] = AnnotationTrack(range=GRanges(seqnames=rep("chrIV",
length(currItems)), ranges=IRanges(start=start, width=width*8, names=values)),
genome=gen, chromosome=chr, name=paste("Viterbi\n", "(", dir, ")", sep=""),
id=ids[order(start)], shape="box",fill=cols[values[order(start)]], col="black",
stacking="dense")
}

for(dir in c("F", "U", "R")) {
myPos = yeastTF_probeAnno_ex$chr04 >= 1217060 & yeastTF_probeAnno_ex$chr04 <= 1225000
myRle = Rle(Data_G$viterbi[myPos])
currItems = which(myRle@values %in% myViterbiDirs[[dir]])

start = yeastTF_probeAnno_ex$chr04[myPos][start(myRle)][currItems]
width = myRle@lengths[currItems]
ids = as.character(myRle@values[currItems])
values = as.character(myRle@values[currItems])
myHiddenStates[[dir]] = AnnotationTrack(range=GRanges(seqnames=rep("chrIV",
length(currItems)), ranges=IRanges(start=start, width=width*8, names=values)),
genome=gen, chromosome=chr, name=paste("Hidden States\n", "(", dir, ")", sep=""),
id=ids[order(start)], shape="box",fill=cols[values[order(start)]], col="black",
stacking="dense")
}



sizes = rep(1,16)

sizes[12:16] = 0.7

plotTracks(c( dlist, myViterbiPanels, myHiddenStates),
from=1217060, to=1225000, sizes=sizes, showFeatureId=TRUE, featureAnnotation="id",
fontcolor.feature="black", cex.feature=0.7, background.title="darkgrey", showId=TRUE)

