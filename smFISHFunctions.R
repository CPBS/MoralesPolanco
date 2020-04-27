#Rules - This detailed dataframe should be the output from FISHquant that has been scraped with the python script.
#Rules - Currently can compare up to 3 different channels

require(gridExtra)
require(grid)
require(ggplot2)
require(plyr)
require(dplyr)
require(stringr)
require(plotly)
require(ggbeeswarm)
require(colorspace)
require(RColorBrewer)

mround <- function(x,base){
  base*round(x/base)
}

plotFISHfilt <- function(smFISH, type="plotly"){     #####Can plot the data as histograms based on several aspects that can be used to infer spot quality - use filtFISH to filter based on these
  if(type=="plotly"){
    p1 = plotly_build(ggplot(smFISH, aes(BGD)) + geom_density())
    p2 = plotly_build(ggplot(smFISH, aes(RES)) + geom_density())
    p3 = plotly_build(ggplot(smFISH, aes(SC_det_norm)) + geom_density())
    p4 = plotly_build(ggplot(smFISH, aes(INT_filt)) + geom_density())
    return(subplot(list(p1,p2,p3,p4), shareX = FALSE, shareY = FALSE, nrows = 2))
  }
  else if(type=="ggplot"){
    p1 = ggplot(smFISH, aes(BGD)) + geom_density() + theme_Publication_GridlinesTesting()+theme(plot.margin = unit(c(0.2, 0.3, 0.2, 0.3), "cm"))
    p2 = ggplot(smFISH, aes(RES)) + geom_density() + theme_Publication_GridlinesTesting()+theme(plot.margin = unit(c(0.2, 0.3, 0.2, 0.3), "cm"))
    p3 = ggplot(smFISH, aes(SC_det_norm)) + geom_density() + theme_Publication_GridlinesTesting()+theme(plot.margin = unit(c(0.2, 0.2, 0.3, 0.3), "cm"))
    p4 = ggplot(smFISH, aes(INT_filt)) + geom_density() + theme_Publication_GridlinesTesting()+theme(plot.margin = unit(c(0.2, 0.3, 0.2, 0.3), "cm"))
    g1 = ggplotGrob(p1)
    g2 = ggplotGrob(p2)
    g3 = ggplotGrob(p3)
    g4 = ggplotGrob(p4)
    g = rbind(g1, g2,g3,g4, size = "first")
    g$widths <- unit.pmax(g1$widths, g2$widths,g3$widths, g4$widths)
    grid.newpage()
    return(grid.arrange(g,top=sprintf("%s",deparse(substitute(smFISH)))))
  }
}
########change to plotly for better plot interaction!

filtFISH <- function(smFISH, intensity=0, score=0, background=Inf, residual=Inf){ ####Filter the data to show only 'real' spots
  return(smFISH[smFISH$INT_filt > intensity  & smFISH$SC_det_norm >score & smFISH$BGD<background &smFISH$RES<residual, ])
}
previewfiltFISH <- function(smFISH, intensity=0, score=0, background=Inf, residual=Inf){  ####Preview the effects of filtering on the data.
  smFISHfilt = smFISH[smFISH$INT_filt > intensity  & smFISH$SC_det_norm >score & smFISH$BGD<background &smFISH$RES<residual, ]
  plotFISHfilt(smFISHfilt)
}
'%!in%' <- function(x,y)!('%in%'(x,y)) ####Function for using inverse of in.

getSpotsXYZ <- function(ch01, ch02, ch03, distDim) {  ####Simplify the data so its just a list of X,Y,Z co-ordinates for the distance calculation - distDim can be 2 or 3 so
  if(distDim == 3) {
    if(!is.null(ch01) && length(ch01[,1])>0){
      ch01Pos = data.frame(X = ch01$Pos_X, Y = ch01$Pos_Y, Z = ch01$Pos_Z)
    }
    else{
      ch01Pos = data.frame(X=NA, Y=NA, Z=NA)
    }
    if(!is.null(ch02) && length(ch02[,1])>0){
      ch02Pos = data.frame(X = ch02$Pos_X, Y = ch02$Pos_Y, Z = ch02$Pos_Z)
    }
    else{
      ch02Pos = data.frame(X=NA, Y=NA, Z=NA)
    }
    if(!is.null(ch03) && length(ch01[,1])>0){
      ch03Pos = data.frame(X = ch03$Pos_X, Y = ch03$Pos_Y, Z = ch03$Pos_Z)
    }
    else{
      ch03Pos = data.frame(X=NA, Y=NA, Z=NA)
    }
  }
  if(distDim == 2) {
    if(!is.null(ch01) && length(ch01[,1])>0){
        ch01Pos = data.frame(X = ch01$Pos_X, Y = ch01$Pos_Y)
    }
    else{
      ch01Pos = data.frame(X=NA, Y=NA)
    }
    if(!is.null(ch02) && length(ch02[,1])>0){
      ch02Pos = data.frame(X = ch02$Pos_X, Y = ch02$Pos_Y)
    }
    else{
      ch02Pos = data.frame(X=NA, Y=NA)
    }
    if(!is.null(ch03) && length(ch03[,1])>0){
      ch03Pos = data.frame(X = ch03$Pos_X, Y = ch03$Pos_Y)
    }
    else{
      ch03Pos = data.frame(X=NA, Y=NA)
    }
  }
  rownames(ch02Pos) = paste(rownames(ch02),"ch02",sep = "")
  rownames(ch03Pos) = paste(rownames(ch03),"ch03",sep = "")
  rownames(ch01Pos) = paste(rownames(ch01),"ch01",sep = "")
  spotPosAll = rbind(ch01Pos, ch02Pos, ch03Pos)
  return(spotPosAll)
}

#getSpotXYZ <- function(spots, distDim, title=""){
#	if(distDim == 3) {
#		Pos = data.frame(X = spots$Pos_X, Y = spots$Pos_Y, Z = spots$Pos_Z)
#	}#
#	if(distDim == 2) {
#		Pos = data.frame(X = spots$Pos_X, Y = spots$Pos_Y)
#	}
#	rownames(ch02Pos) = paste(rownames(Pos),title,sep = "")
#	return(Pos)
#}

closestSpotGrabber <- function(dm = dm, i, numAgainst = numAgainst, numAll = numAll,spotsCompAgainst = spotsCompAgainst, spotsComp=spotsComp){
  #distanceMatrix is the big matrix of every spot vs every spot
  #spotPosAll is the dataframe with the X, Y, Z pos of both of the sets of spots to be compared - e.g. spotPos1v2
  #spotPosAgainst is the dataframe with the X,Y,Z pos of only the set of spots that act as the reference e.g if doing 1v2, this will be spotPosCh01 we will get this in the function
  #comp will be the comparison e.g. "1v2" or "2v3"
  closestSpotsDF = data.frame()
  if(length(dm[(numAgainst+1):numAll, i])>1){
  closestSpot = dm[(numAgainst+1):numAll, i][which(dm[(numAgainst+1):numAll, i]==min(dm[(numAgainst+1):numAll, i]))]
  }
  if(length(dm[(numAgainst+1):numAll, i])==1){
    closestSpot = as.data.frame(dm[(numAgainst+1):numAll, i][which(dm[(numAgainst+1):numAll, i]==min(dm[(numAgainst+1):numAll, i]))])
    colnames(closestSpot) = colnames(dm)[length(colnames(dm))]
  }
  print(closestSpot)
  spotNumber = as.numeric(substr(names(closestSpot),-4,nchar(names(closestSpot))-4))
  print(spotNumber)
  df = data.frame(spotsCompAgainst=str_remove(rownames(dm), "ch0[1-3]")[1], spotsComp=spotNumber, Distance = closestSpot[[1]])
  print(df)
  #Change the column names of the df to the variable e.g. before the df would have a column called spotsComp, but what
  #we actually want is the column name to be the string stored in the variable spotsComp (e.g. Ch02)
  names(df)[names(df)=="spotsCompAgainst"]  <- spotsCompAgainst
  names(df)[names(df)=="spotsComp"]  <- spotsComp
  closestSpotsDF = rbind(closestSpotsDF, df)
  return(closestSpotsDF)
}



getClosestSpots <- function(dm, comp){       ############Returns a df containing each spot paired with its closest spot in the opposite channel - crude output e.g. ch01 spot1, closest spot in ch02 spot45. Use getDetails to get detailed spot info for paired spots
  spotsCompAgainst = paste("ch0", substr(comp, 1, 1), sep = "")
  spotsComp=paste("ch0", substr(comp, 3, 3), sep = "")
  numAgainst=length(subset(dm, subset = str_sub(rownames(dm), -4) %in% spotsCompAgainst)[,1])
  numAll=length(dm[,1])
  return(do.call(rbind.data.frame,lapply(1:numAgainst,closestSpotGrabber, dm =dm, numAgainst = numAgainst, numAll = numAll, spotsCompAgainst = spotsCompAgainst, spotsComp=spotsComp)))
}

getDetails <- function(pairedSpots, comp,distDim) {   ###############Gets all spot info for paired spots e.g. intensity, x,y,z, size, etc.
  #1st, make 2 variables based on the comp, that can be used to get the appropriate dataframes without having to supply a huge list of vairables each time the function is called
  #This does require thought that the name of the detailed dataframe is simply ch01 or ch02 etc.
  #Paired spots is the output from function getClosestSpots (it is just a list of closest spots and distance between them)
  #i just needs to be the entire length of the pairedSpots
  i = 1:length(pairedSpots[,1])
  spotsCompAgainst = paste("ch0", substr(comp, 1, 1), sep = "")
  spotsComp = paste("ch0", substr(comp, 3, 3), sep = "")
  #Now get these dfs and assign them to the appropriate variable
  detailedAgainst = get(spotsCompAgainst)
  detailedComp = get(spotsComp)
  #Take each spot number - this is just listed in each row of the pairedSpots df
  #and then find the matching spot in the original detailed dataframes and get ALL of the data for that spot
  compAgainst = detailedAgainst[pairedSpots[i,1], ]
  comp = detailedComp[pairedSpots[i,2], ]
  if(distDim ==3){
    compAgainst$size = (compAgainst$SigmaX+compAgainst$SigmaY+compAgainst$SigmaZ)/3
    comp$size = (comp$SigmaX+comp$SigmaY+comp$SigmaZ)/3
  }
  if(distDim ==2){
    compAgainst$size = (compAgainst$SigmaX+compAgainst$SigmaY)/2
    comp$size = (comp$SigmaX+comp$SigmaY)/2
  }
  #Change the column names of these so we can easily figure out which spot attributes are associated with which channel
  colnames(compAgainst) = paste(colnames(compAgainst), spotsCompAgainst, sep = "")
  colnames(comp) = paste(colnames(comp), spotsComp, sep = "")
  return(cbind(compAgainst, comp, pairedSpots))
}

colocCalc<- function(pairedDetailed, step=5,maxColoc=NULL){ 
  #####Calculate co-localisation over a sliding window. 
  #####Function will return proportion of co-localising spots based on a range of distances between 0nm and maxColoc nm. 
  #####If left blank, maxColoc will default to the max distance between spots (giving 100% co-localisation)
  ##### step is the step size, e.g. increment co-localisation distance by 5nm, 10nm etc. Defaults to 5nm. - if set to "max/100" then step will be max distance/100.
  if(step=="max/100"){
    step=(max(pairedDetailed$Distance)/100)
  }
  if(is.null(maxColoc)){
    coloc = lapply(seq(from=0, to = mround(max(pairedDetailed$Distance),base=step), by = step), function(i, distData=pairedDetailed$Distance) length(distData[distData<i])/length(distData))
    dfColoc = (as.data.frame(do.call(rbind, coloc)))
    return(data.frame(colocalisation = dfcoloc$V1, distance = seq(from=0, to = mround(max(pairedDetailed$Distance),base=step), by = step)))
  }
  else{
    coloc = lapply(seq(from=0, to = mround(maxColoc,base=step), by = step), function(i, distData=pairedDetailed$Distance) length(distData[distData<i])/length(distData))
    dfColoc = (as.data.frame(do.call(rbind, coloc)))
    return(data.frame(colocalisation = dfcoloc$V1, distance = seq(from=0, to = mround(maxColoc,base=step),by=step)))
  }
}

compSpots <- function(ch01=NULL, ch02=NULL, ch03=NULL, pixelsize = 160, distDim = 2, step=5){ ######Compares spots from channels in pairwise comparisons e.g. 1v2, 1v3, 2v3 etc. V slow unless data has been filtered to remove 'low quality' spots.
  #ch01 etc is the output from FISHQuant that has been scraped ending in scraped detailed - it can be filtered beforehand using the filtFISH functions
  #pixelsize is only used if the correct pixelsize was not set in matlab
  #distDim is the number of dimensions used to calculate distance, 2 = xy, 3 =xyz
  #step is the size of the step between distances used to calculate colocalisation, default = 5, may be more appropriate to use "max/100".
  ch01<<-ch01
  ch02<<-ch02
  ch03<<-ch03
  if(!is.null(ch01)){
    ch01$Pos_X = ch01$Pos_X/(160/pixelsize)
    ch01$Pos_Y = ch01$Pos_Y/(160/pixelsize)
    ch01$Pos_Z = ch01$Pos_Z/(160/pixelsize)
  }
  if(!is.null(ch02)){
    ch02$Pos_X = ch02$Pos_X/(160/pixelsize)
    ch02$Pos_Y = ch02$Pos_Y/(160/pixelsize)
    ch02$Pos_Z = ch02$Pos_Z/(160/pixelsize)
  }
  if(!is.null(ch03)){
    ch03$Pos_X = ch03$Pos_X/(160/pixelsize)
    ch03$Pos_Y = ch03$Pos_Y/(160/pixelsize)
    ch03$Pos_Z = ch03$Pos_Z/(160/pixelsize)
  }
  spotPosAll = getSpotsXYZ(ch01, ch02, ch03, distDim)
  if(is.null(ch01)){
    spotPosAll = na.omit(spotPosAll)
    dist2v3 = as.matrix(dist(spotPosAll))
    paired2v3 = getClosestSpots(dm = dist2v3, comp = "2v3")
    pairedDetailed2v3 = getDetails(pairedSpots = paired2v3, comp = "2v3",distDim=distDim)
    coloc2v3 = colocCalcStep(pairedDetailed = pairedDetailed2v3, step=step)
    combColocStep = coloc2v3
    #coloc2v3 = data.frame(colocalisation = length(which(pairedDetailed2v3$Distance<((pairedDetailed2v3$sizech02+pairedDetailed2v3$sizech03)/2))), comp="2v3")
    coloc2v3spot = data.frame(colocalisation = length(which(pairedDetailed2v3$Distance<((pairedDetailed2v3$sizech02+pairedDetailed2v3$sizech03)/2)))/length(pairedDetailed2v3[,1]), comp="2v3")
    combColoc = coloc2v3spot
    return(list(combColoc, combColocStep,pairedDetailed2v3))
  }
  if(is.null(ch02)){
    spotPosAll = na.omit(spotPosAll)
    dist1v3 = as.matrix(dist(spotPosAll))
    paired1v3 = getClosestSpots(dm = dist1v3, comp = "1v3")
    pairedDetailed1v3 = getDetails(pairedSpots = paired1v3, comp = "1v3",distDim=distDim)
    coloc1v3 = colocCalcStep(pairedDetailed = pairedDetailed1v3, step=step)
    combColocStep = coloc1v3
    #coloc1v3 = data.frame(colocalisation = length(which(pairedDetailed1v3$Distance<((pairedDetailed1v3$sizech01+pairedDetailed1v3$sizech03)/2))), comp="1v3")
    coloc1v3spot = data.frame(colocalisation = length(which(pairedDetailed1v3$Distance<((pairedDetailed1v3$sizech01+pairedDetailed1v3$sizech03)/2)))/length(pairedDetailed1v3[,1]), comp="1v3")
    combColoc = coloc1v3spot
    return(list(combColoc, combColocStep,pairedDetailed1v3))
  }
  if(is.null(ch03)){
    spotPosAll = na.omit(spotPosAll)
    dist1v2 = as.matrix(dist(spotPosAll))
    paired1v2 = getClosestSpots(dm = dist1v2, comp = "1v2")
    pairedDetailed1v2 = getDetails(pairedSpots = paired1v2, comp = "1v2",distDim=distDim)
    coloc1v2 = colocCalcStep(pairedDetailed = pairedDetailed1v2, step=step)
    combColocStep = coloc1v2
    #coloc1v2 = data.frame(colocalisation = length(which(pairedDetailed1v2$Distance<((pairedDetailed1v2$sizech01+pairedDetailed1v2$sizech02)/2))), comp="1v2")
    coloc1v2spot = data.frame(colocalisation = length(which(pairedDetailed1v2$Distance<((pairedDetailed1v2$sizech01+pairedDetailed1v2$sizech02)/2)))/length(pairedDetailed1v2[,1]), comp="1v2")
    combColoc = coloc1v2spot
    return(list(combColoc, combColocStep,pairedDetailed1v2))
  }
  else {
    spotPos1v2 = subset(spotPosAll, subset = str_sub(rownames(spotPosAll), -4) %!in% "ch03")
    spotPos1v3 = subset(spotPosAll, subset = str_sub(rownames(spotPosAll), -4) %!in% "ch02")
    spotPos2v3 = subset(spotPosAll, subset = str_sub(rownames(spotPosAll), -4) %!in% "ch01")
    #Calculate the distance between spots
    dist1v2 = as.matrix(dist(spotPos1v2))
    dist1v3 = as.matrix(dist(spotPos1v3))
    dist2v3 = as.matrix(dist(spotPos2v3))
    paired1v2 = getClosestSpots(dm = dist1v2, comp = "1v2")
    paired1v3 = getClosestSpots(dm = dist1v3, comp = "1v3")
    paired2v3 = getClosestSpots(dm = dist2v3, comp = "2v3")
    pairedDetailed1v2 = getDetails(pairedSpots = paired1v2, comp = "1v2", distDim=distDim)
    pairedDetailed1v3 = getDetails(pairedSpots = paired1v3, comp = "1v3",distDim=distDim)
    pairedDetailed2v3 = getDetails(pairedSpots = paired2v3, comp = "2v3",distDim=distDim)
    #Generate stepwise co-localisation data
    coloc1v2 = colocCalcStep(pairedDetailed = pairedDetailed1v2, step=step)
    coloc1v3 = colocCalcStep(pairedDetailed = pairedDetailed1v3, step=step)
    coloc2v3 = colocCalcStep(pairedDetailed = pairedDetailed2v3, step=step)
    combColocStep = rbind(coloc1v2, coloc1v3, coloc2v3)
    #Generate per-spot co-localisation data, which takes into account spot width etc.
    coloc1v2spot = data.frame(colocalisation = length(which(pairedDetailed1v2$Distance<((pairedDetailed1v2$sizech01+pairedDetailed1v2$sizech02)/2)))/length(pairedDetailed1v2[,1]), comp="1v2")
    coloc1v3spot = data.frame(colocalisation = length(which(pairedDetailed1v3$Distance<((pairedDetailed1v3$sizech01+pairedDetailed1v3$sizech03)/2)))/length(pairedDetailed1v2[,1]), comp="1v3")
    coloc2v3spot = data.frame(colocalisation = length(which(pairedDetailed2v3$Distance<((pairedDetailed2v3$sizech02+pairedDetailed2v3$sizech03)/2)))/length(pairedDetailed1v2[,1]), comp="2v3")
    combColoc = rbind(coloc1v2spot, coloc1v3spot, coloc2v3spot)
    return(list(combColoc, combColocStep,pairedDetailed1v2, pairedDetailed1v3, pairedDetailed2v3))
  }
  rm(ch01,ch02,ch03)
  #Get the closest spot in the opposing channel - these spots are then 'paired' #Need if statements so if ch03 doesn't exist, it doesnt try and do it
} #Compare spots in up to 3 channels for co-localization


compSpotsByCell <- function(ch01=NULL, ch02=NULL, ch03=NULL, pixelsize = 160, distDim = 2, step=5){ ######Compares spots from channels in pairwise comparisons e.g. 1v2, 1v3, 2v3 etc. V slow unless data has been filtered to remove 'low quality' spots.
  #ch01 etc is the output from FISHQuant that has been scraped ending in scraped detailed - it can be filtered beforehand using the filtFISH functions
  #pixelsize is only used if the correct pixelsize was not set in matlab
  #distDim is the number of dimensions used to calculate distance, 2 = xy, 3 =xyz
  #step is the size of the step between distances used to calculate colocalisation, default = 5, may be more appropriate to use "max/100".
  ch01<<-ch01
  ch02<<-ch02
  ch03<<-ch03
  if(!is.null(ch01)){
    ch01$Pos_X = ch01$Pos_X/(160/pixelsize)
    ch01$Pos_Y = ch01$Pos_Y/(160/pixelsize)
    ch01$Pos_Z = ch01$Pos_Z/(160/pixelsize)
  }
  if(!is.null(ch02)){
    ch02$Pos_X = ch02$Pos_X/(160/pixelsize)
    ch02$Pos_Y = ch02$Pos_Y/(160/pixelsize)
    ch02$Pos_Z = ch02$Pos_Z/(160/pixelsize)
  }
  if(!is.null(ch03)){
    ch03$Pos_X = ch03$Pos_X/(160/pixelsize)
    ch03$Pos_Y = ch03$Pos_Y/(160/pixelsize)
    ch03$Pos_Z = ch03$Pos_Z/(160/pixelsize)
  }
  if(is.null(ch01)){
    pairedDetailed2v3 = data.frame()
    for(i in 1:max(ch02$Cell)){
      print(paste("Cell",i,sep=" "))
      spotPosAll = getSpotsXYZ(ch01=NULL, ch02[ch02$Cell ==i,], ch03[ch02$Cell==i,], distDim)
      spotPosAll = na.omit(spotPosAll)
      dist2v3 = as.matrix(dist(spotPosAll))
      paired2v3 = getClosestSpots(dm = dist2v3, comp = "2v3")
      pairedDetailedCell2v3 = getDetails(pairedSpots = paired2v3, comp = "2v3",distDim=distDim)
      rbind(pairedDetailed2v3, pairedDetailedCell2v3)
    }
    coloc2v3 = colocCalcStep(pairedDetailed = pairedDetailed2v3, step=step)
    combColocStep = coloc2v3
    #coloc2v3 = data.frame(colocalisation = length(which(pairedDetailed2v3$Distance<((pairedDetailed2v3$sizech02+pairedDetailed2v3$sizech03)/2))), comp="2v3")
    coloc2v3spot = data.frame(colocalisation = length(which(pairedDetailed2v3$Distance<((pairedDetailed2v3$sizech02+pairedDetailed2v3$sizech03)/2)))/length(pairedDetailed2v3[,1]), comp="2v3")
    combColoc = coloc2v3spot
    return(list(combColoc, combColocStep,pairedDetailed2v3))
  }
  if(is.null(ch02)){
    pairedDetailed1v3 = data.frame()
    for(i in 1:max(ch03$Cell)){
      print(paste("Cell",i,sep=" "))
      spotPosAll = getSpotsXYZ(ch01=ch01[ch01$Cell==i,], NULL, ch03[ch02$Cell==i,], distDim)
      spotPosAll = na.omit(spotPosAll)
      dist1v3 = as.matrix(dist(spotPosAll))
      paired1v3 = getClosestSpots(dm = dist1v3, comp = "1v3")
      pairedDetailedCell1v3 = getDetails(pairedSpots = paired1v3, comp = "1v3",distDim=distDim)
      rbind(pairedDetailed1v3, pairedDetailedCell1v3)
    }
    coloc1v3 = colocCalcStep(pairedDetailed = pairedDetailed1v3, step=step)
    combColocStep = coloc1v3
    #coloc1v3 = data.frame(colocalisation = length(which(pairedDetailed1v3$Distance<((pairedDetailed1v3$sizech01+pairedDetailed1v3$sizech03)/2))), comp="1v3")
    coloc1v3spot = data.frame(colocalisation = length(which(pairedDetailed1v3$Distance<((pairedDetailed1v3$sizech01+pairedDetailed1v3$sizech03)/2)))/length(pairedDetailed1v3[,1]), comp="1v3")
    combColoc = coloc1v3spot
    return(list(combColoc, combColocStep,pairedDetailed1v3))
  }
  if(is.null(ch03)){
    pairedDetailed1v2 = data.frame()
    for(i in 1:max(ch01$Cell)){
      print(paste("Cell",i,sep=" "))
      spotPosAll = getSpotsXYZ(ch01=ch01[ch01$Cell==i,], ch02[ch02$Cell ==i,], NULL, distDim)
      spotPosAll = na.omit(spotPosAll)
      dist1v2 = as.matrix(dist(spotPosAll))
      paired1v2 = getClosestSpots(dm = dist1v2, comp = "1v2")
      pairedDetailedCell1v2 = getDetails(pairedSpots = paired1v2, comp = "1v2",distDim=distDim)
      rbind(pairedDetailed1v2, pairedDetailedCell1v2)
    }
    coloc1v2 = colocCalcStep(pairedDetailed = pairedDetailed1v2, step=step)
    combColocStep = coloc1v2
    #coloc1v2 = data.frame(colocalisation = length(which(pairedDetailed1v2$Distance<((pairedDetailed1v2$sizech01+pairedDetailed1v2$sizech02)/2))), comp="1v2")
    coloc1v2spot = data.frame(colocalisation = length(which(pairedDetailed1v2$Distance<((pairedDetailed1v2$sizech01+pairedDetailed1v2$sizech02)/2)))/length(pairedDetailed1v2[,1]), comp="1v2")
    combColoc = coloc1v2spot
    return(list(combColoc, combColocStep,pairedDetailed1v2))
  }
  else {
    pairedDetailed1v2 = data.frame()
    pairedDetailed1v3 = data.frame()
    pairedDetailed2v3 = data.frame()
    pb = txtProgressBar(min=0,max=max(ch01$Cell),style = 3)
    print("Pairing Spots...")
    for(i in 1:max(ch01$Cell)){
      setTxtProgressBar(pb,i)
      spotPosAll = getSpotsXYZ(ch01=ch01[ch01$Cell==i,], ch02[ch02$Cell ==i,], ch03[ch03$Cell ==i,], distDim)
      if(length(ch01[ch01$Cell==i,1]) != 0 && length(ch02[ch02$Cell==i,1]) !=0){
        spotPos1v2 = subset(spotPosAll, subset = str_sub(rownames(spotPosAll), -4) %!in% "ch03")
        dist1v2 = as.matrix(dist(spotPos1v2))
        paired1v2 = getClosestSpots(dm = dist1v2, comp = "1v2")
        pairedDetailedCell1v2 = getDetails(pairedSpots = paired1v2, comp = "1v2",distDim=distDim)
        rbind(pairedDetailed1v2, pairedDetailedCell1v2)
        print(paired1v2)
      }
      if(length(ch01[ch01$Cell==i,1]) != 0 && length(ch03[ch03$Cell==i,1]) !=0){
        spotPos1v3 = subset(spotPosAll, subset = str_sub(rownames(spotPosAll), -4) %!in% "ch02")
        dist1v3 = as.matrix(dist(spotPos1v3))
        paired1v3 = getClosestSpots(dm = dist1v3, comp = "1v3")
        pairedDetailedCell1v3 = getDetails(pairedSpots = paired1v3, comp = "1v3",distDim=distDim)
        rbind(pairedDetailed1v3, pairedDetailedCell1v3)
      }
      if(length(ch03[ch03$Cell==i,1]) != 0 && length(ch02[ch02$Cell==i,1]) !=0){
        spotPos2v3 = subset(spotPosAll, subset = str_sub(rownames(spotPosAll), -4) %!in% "ch01")
        dist2v3 = as.matrix(dist(spotPos2v3))
        paired2v3 = getClosestSpots(dm = dist2v3, comp = "2v3")
        pairedDetailedCell2v3 = getDetails(pairedSpots = paired2v3, comp = "2v3",distDim=distDim)
        rbind(pairedDetailed2v3, pairedDetailedCell2v3)
      }
    }
    print("Spots Paired")
    #Generate stepwise co-localisation data
    print("Step Coloc Calculations...")
    pairedDetailed1v2<<-pairedDetailed1v2
    coloc1v2 = colocCalcStep(pairedDetailed = pairedDetailed1v2, step=step)
    coloc1v3 = colocCalcStep(pairedDetailed = pairedDetailed1v3, step=step)
    coloc2v3 = colocCalcStep(pairedDetailed = pairedDetailed2v3, step=step)
    combColocStep = rbind(coloc1v2, coloc1v3, coloc2v3)
    #Generate per-spot co-localisation data, which takes into account spot width etc.
    print("Per Channel Spot Colocalisation...")
    coloc1v2spot = data.frame(colocalisation = length(which(pairedDetailed1v2$Distance<((pairedDetailed1v2$sizech01+pairedDetailed1v2$sizech02)/2)))/length(pairedDetailed1v2[,1]), comp="1v2")
    coloc1v3spot = data.frame(colocalisation = length(which(pairedDetailed1v3$Distance<((pairedDetailed1v3$sizech01+pairedDetailed1v3$sizech03)/2)))/length(pairedDetailed1v2[,1]), comp="1v3")
    coloc2v3spot = data.frame(colocalisation = length(which(pairedDetailed2v3$Distance<((pairedDetailed2v3$sizech02+pairedDetailed2v3$sizech03)/2)))/length(pairedDetailed1v2[,1]), comp="2v3")
    combColoc = rbind(coloc1v2spot, coloc1v3spot, coloc2v3spot)
    return(list(combColoc, combColocStep,pairedDetailed1v2, pairedDetailed1v3, pairedDetailed2v3))
  }
  rm(ch01,ch02,ch03)
  print("Done")
  #Get the closest spot in the opposing channel - these spots are then 'paired' #Need if statements so if ch03 doesn't exist, it doesnt try and do it
} #Compare spots in up to 3 channels for co-localization

compSpotsRange <- function(ch01=NULL, ch02=NULL, ch03=NULL, pixelsize=160, distDim=2, step=5, thresholdRange=NULL, filtVar = "SC_det_norm"){
  #filtVar defines the column to filter based upon - defaults to SC_det_norm.
  #Performing for 3 channels will take a very long time unless strict filtering is applied. Do not run without filtering!
  #thresholdRange can be a list of cut off values to test the effect of different cut offs for each channel.
  #It should be as follows in order of thresholds for ch01, ch02, ch03 thresholdRange=list(c(0,0.2,0.5), c(0,0.1,0.8), c(0,0.4,0.9))
  if(is.null(thresholdRange) | !is.list(thresholdRange)){
    print("Must supply threshold range as a list of cut off values for filtering")
  }
  if(length(c(ch01[1,1], ch02[1,1],ch03[1,1]))!=length(thresholdRange)){
    print("Error, number of lists for filtering does not equal number of channels being analysed\nif you do not want to apply filtering, simply make the filtering NA")
  }
  ########IN PROGRESS
  filteringMatrix =expand.grid(thresholdRange[[1]],thresholdRange[[2]], thresholdRange[[3]])
  #remove filtering for ch0x if list = NA
  filteringMatrix[sapply(filteringMatrix, function(filteringMatrix) !any(is.na(filteringMatrix)))]
  if(!("ch01" %in% paste("ch0", (substr(names(filteringMatrix),4,5)),sep=""))){
    filtch01=ch01
  }
  if(!("ch02" %in% paste("ch0", (substr(names(filteringMatrix),4,5)),sep=""))){
    filtch02=ch02
  }
  if(!("ch03" %in% paste("ch0", (substr(names(filteringMatrix),4,5)),sep=""))){
    filtch03=ch03
  }
  listComps = list()
  listColoc = list()
  for(i in 1:length(filteringMatrix[,1])){
    argnames <- match.call()
    ch01<<-get(paste(as.list(argnames)$ch01))
    ch02<<-get(paste(as.list(argnames)$ch02))
    ch03<<-get(paste(as.list(argnames)$ch03))
    print(sprintf("Cycle %iof%i", i, length(filteringMatrix[,1])))
    for(z in 1:length(filteringMatrix)){
      assign(paste("filt",paste("ch0",substr(names(filteringMatrix),4,5)[z],sep=""),sep=""),get(paste("ch0",substr(names(filteringMatrix),4,5)[z],sep=""))[get((paste("ch0",substr(names(filteringMatrix),4,5)[z],sep="")))[filtVar]>filteringMatrix[i,z], ])
    }
    compColoc=compSpots(ch01=filtch01, ch02=filtch02,ch03=filtch03, step=step, distDim = distDim, pixelsize = pixelsize)
    listComps[[i]]=compColoc[[2]]
    listComps[[i]]$FiltID = paste(names(filteringMatrix),filteringMatrix[i,], collapse='', sep="_")
    listColoc[[i]]=compColoc[[1]]
    listColoc[[i]]$FiltID = paste(names(filteringMatrix),filteringMatrix[i,], collapse='', sep="_")
  }
  dfColoc=do.call(rbind, listColoc)
  dfComps=do.call(rbind, listComps)
  return(list(dfColoc, dfComps, pairedDetailed1v2, pairedDetailed1v3, pairedDetailed2v3))
} #Iterates through comp spots for varying SC_det_norm thresholds. Simple way to view impact of filtering on spot pairing.

compSpotsRangeByCell <- function(ch01=NULL, ch02=NULL, ch03=NULL, pixelsize=160, distDim=2, step=5, thresholdRange=NULL, filtVar = "SC_det_norm"){
  #filtVar defines the column to filter based upon - defaults to SC_det_norm.
  #Performing for 3 channels will take a very long time unless strict filtering is applied. Do not run without filtering!
  #thresholdRange can be a list of cut off values to test the effect of different cut offs for each channel.
  #It should be as follows in order of thresholds for ch01, ch02, ch03 thresholdRange=list(c(0,0.2,0.5), c(0,0.1,0.8), c(0,0.4,0.9))
  if(is.null(thresholdRange) | !is.list(thresholdRange)){
    print("Must supply threshold range as a list of cut off values for filtering")
  }
  if(length(c(ch01[1,1], ch02[1,1],ch03[1,1]))!=length(thresholdRange)){
    print("Error, number of lists for filtering does not equal number of channels being analysed\nif you do not want to apply filtering, simply make the filtering NA")
  }
  ########IN PROGRESS
  filteringMatrix =expand.grid(thresholdRange[[1]],thresholdRange[[2]], thresholdRange[[3]])
  #remove filtering for ch0x if list = NA
  filteringMatrix[sapply(filteringMatrix, function(filteringMatrix) !any(is.na(filteringMatrix)))]
  if(!("ch01" %in% paste("ch0", (substr(names(filteringMatrix),4,5)),sep=""))){
    filtch01=ch01
  }
  if(!("ch02" %in% paste("ch0", (substr(names(filteringMatrix),4,5)),sep=""))){
    filtch02=ch02
  }
  if(!("ch03" %in% paste("ch0", (substr(names(filteringMatrix),4,5)),sep=""))){
    filtch03=ch03
  }
  listComps = list()
  listColoc = list()
  for(i in 1:length(filteringMatrix[,1])){
    argnames <- match.call()
    ch01<<-get(paste(as.list(argnames)$ch01))
    ch02<<-get(paste(as.list(argnames)$ch02))
    ch03<<-get(paste(as.list(argnames)$ch03))
    print(sprintf("Cycle %iof%i", i, length(filteringMatrix[,1])))
    for(z in 1:length(filteringMatrix)){
      assign(paste("filt",paste("ch0",substr(names(filteringMatrix),4,5)[z],sep=""),sep=""),get(paste("ch0",substr(names(filteringMatrix),4,5)[z],sep=""))[get((paste("ch0",substr(names(filteringMatrix),4,5)[z],sep="")))[filtVar]>filteringMatrix[i,z], ])
    }
    compColoc=compSpotsByCell(ch01=filtch01, ch02=filtch02,ch03=filtch03, step=step, distDim = distDim, pixelsize = pixelsize)
    listComps[[i]]=compColoc[[2]]
    listComps[[i]]$FiltID = paste(names(filteringMatrix),filteringMatrix[i,], collapse='', sep="_")
    listColoc[[i]]=compColoc[[1]]
    listColoc[[i]]$FiltID = paste(names(filteringMatrix),filteringMatrix[i,], collapse='', sep="_")
  }
  dfColoc=do.call(rbind, listColoc)
  dfComps=do.call(rbind, listComps)
  return(list(dfColoc, dfComps, pairedDetailed1v2, pairedDetailed1v3, pairedDetailed2v3))
} #Iterates through comp spots for varying SC_det_norm thresholds. Simple way to view impact of filtering on spot pairing.

stepPlot<-function(comparedSpots, colour=NULL, collated=FALSE, shareAxes = TRUE, size=1){
  #####comparedSpots is the 2nd output from compSpots function e.g. combColoc[[2]]
  #####colour can be a list with multiple colours e.g. colours = list(c("red","blue"), c("green","yellow")). Then each list can be applied to a specific comparison
  #####this might be particularly useful for instances where 3 pairwise comparisons have been performed. In this instance, the step plot will be a gradient between these values.
  #####for colour to be a list of gradients, the data must be output from compSpotsRange. IT WILL NOT WORK WITH compSpots output.
  #####collated refers to the output style of the plot. All pairwise comparisons can be plotted on the same graph, alternatively, if collated=FALSE,
  #####each pairwise comparison will have its own plot.
  #####shareAxes dictates whether these plots can share axes, which may aid visualisation of differential co-localisation. defualt = TRUE.
  #####Size sets the thickness of the line. default =1
  colocPropPlot = ggplot(comparedSpots, aes(distance, colocalisation), group = comp)+ geom_line(aes(col=comp),size=size) + theme_Publication_Gridlines()+labs(x="Distance (nm)", y = "Colocalisation (%)")
  if(collated==FALSE & shareAxes==TRUE){
    colocPropPlot=colocPropPlot+facet_grid(comp~.)
  }
  if(collated==FALSE & shareAxes==FALSE){
    colocPropPlot=colocPropPlot+facet_wrap(comp~.,scales = "free",ncol = 1, strip.position = "right")
  }
  if(collated==TRUE & shareAxes==TRUE){
    shareAxes=FALSE
    print("Cannot share axes on collated graph, setting shareAxes to FALSE")
  }
  if(collated==FALSE & shareAxes==FALSE){
    colocPropPlot=colocPropPlot
  }
  if(is.list(colour)){
    grobs = list()
    for(i in 1:length(colour)){
      colname = sprintf("cols%s",i)
      grobname = sprintf("grob%s",i)
      assign(colname, colorRampPalette(colour[[i]]))
      assign(grobname, ggplot(subset(comparedSpots,comparedSpots$comp==unique(comparedSpots$comp)[i]), aes(distance, colocalisation), group = comp)+geom_line(aes(col=FiltID),size=size)+facet_grid(comp~.)
             +scale_color_manual(values = get(colname)(length(unique(comparedSpots$FiltID)))) + labs(x="Distance (nm)", y = "Colocalisation (%)")+theme_Publication_Gridlines_CB()+guides(col=guide_legend(nrow=2)))
      grobs[[i]]= get(grobname)
      colocPropPlot = grid.arrange(grobs=grobs)
      #######NEED TO MAKE IT SO THAT CAN HAVE COLOUR AS LIST, BUT COLLATED - and also one that shares axes - i think just add a line with scale_x_continuous - min and max to all of them if shareAxes = TRUE!
    }
  }
  if(is.list(colour)&!is.null(colour)){
    if(length(unique(comparedSpots$comp))!=length(colour)){
      print(paste("Number of colours supplied != number of comparisons, number of comparisons is n=", length(unique(comparedSpots$comp)), sep=""))
    }
    else {
      colocPropPlot = colocPropPlot + scale_color_manual(values = colour)
    }
  }
  return(colocPropPlot)
} #Plot of spot co-localization %age over different distances. E.g. spots considered to co-loc if centroids are within 10nm, 20nm,....,1000nm etc.

colocCalcStep<- function(pairedDetailed, step=5,maxColoc=NULL){ 
  #####Calculate co-localisation over a sliding window.
  #####Function will return proportion of co-localising spots based on a range of distances between 0nm and maxColoc nm.
  #####If left blank, maxColoc will default to the max distance between spots (giving 100% co-localisation)
  ##### step is the step size, e.g. increment co-localisation distance by 5nm, 10nm etc. Defaults to 5nm. - if set to "max/100" then step will be max distance/100, providing 100 steps.
  if(step=="max/100"){
    step=(max(pairedDetailed$Distance)/100)
    print(step)
    print(mround(max(pairedDetailed$Distance)))
  }
  if(is.null(maxColoc)){
    coloc = lapply(seq(from=0, to = mround(max(pairedDetailed$Distance),base=step), by = step), function(i, distData=pairedDetailed$Distance) length(distData[distData<i])/length(distData))
    dfColoc = (as.data.frame(do.call(rbind, coloc)))
    return(data.frame(colocalisation = dfColoc$V1, distance = seq(from=0, to = mround(max(pairedDetailed$Distance),base=step), by = step),comp=rep(sprintf("%sv%s",substr(names(pairedDetailed[1]),5,nchar(names(pairedDetailed))),  substr(names(pairedDetailed[35]),5,nchar(names(pairedDetailed)))),length(dfColoc$V1))))
  }
  else{
    coloc = lapply(seq(from=0, to = mround(maxColoc,base=step), by = step), function(i, distData=pairedDetailed$Distance) length(distData[distData<i])/length(distData))
    dfColoc = (as.data.frame(do.call(rbind, coloc)))
    return(data.frame(colocalisation = dfColoc$V1, distance = seq(from=0, to = mround(maxColoc,base=step),by=step),comp=rep(sprintf("%sv%s",substr(names(pairedDetailed[1]),5,nchar(names(pairedDetailed))),  substr(names(pairedDetailed[35]),5,nchar(names(pairedDetailed)))),length(dfColoc$V1))))
  }
}

cellColocalised <- function(compSpotsDF, leniency=0){
  compSpotsDF[[3]]$ColocalisedStrict = 0
  compSpotsDF[[3]]$ColocalisedStrict[which((compSpotsDF[[3]]$sizech01+compSpotsDF[[3]]$sizech02)/2>compSpotsDF[[3]]$Distance)] =1
  compSpotsDF[[3]]$ColocalisedLenient = 0
  compSpotsDF[[3]]$ColocalisedLenient[which(((compSpotsDF[[3]]$sizech01+compSpotsDF[[3]]$sizech02)/2)+leniency>compSpotsDF[[3]]$Distance)] =1
  
  compSpotsDF[[3]]$Colocalised2Strict = 0
  compSpotsDF[[3]]$Colocalised2Strict[which(((compSpotsDF[[3]]$sizech01/2)>compSpotsDF[[3]]$Distance)|((compSpotsDF[[3]]$sizech02/2)>compSpotsDF[[3]]$Distance))] =1
  compSpotsDF[[3]]$Colocalised2StrictLeniency = 0
  compSpotsDF[[3]]$Colocalised2StrictLeniency[which(((compSpotsDF[[3]]$sizech01/2)+leniency>compSpotsDF[[3]]$Distance)|((compSpotsDF[[3]]$sizech02/2)+leniency>compSpotsDF[[3]]$Distance))] =1
  compSpotsDF[[3]]$Colocalised2xStrict = 0
  compSpotsDF[[3]]$Colocalised2xStrict[which(((compSpotsDF[[3]]$sizech01/2)>compSpotsDF[[3]]$Distance)&((compSpotsDF[[3]]$sizech02/2)>compSpotsDF[[3]]$Distance))] =1
  compSpotsDF[[3]]$Colocalised2xStrictLeniency = 0
  compSpotsDF[[3]]$Colocalised2xStrictLeniency[which(((compSpotsDF[[3]]$sizech01/2)+leniency>compSpotsDF[[3]]$Distance)&((compSpotsDF[[3]]$sizech02/2)+leniency>compSpotsDF[[3]]$Distance))] =1
  
  proportions=data.frame()
  for(i in 1:length(unique(compSpotsDF[[3]]$Cellch01))){
    df=data.frame(cell=i, 
                  propStrict =length(which(compSpotsDF[[3]]$ColocalisedStrict==1&compSpotsDF[[3]]$Cellch01==i))/length(which(compSpotsDF[[3]]$Cellch01==i)), 
                  propLenient=length(which(compSpotsDF[[3]]$ColocalisedLenient==1&compSpotsDF[[3]]$Cellch01==i))/length(which(compSpotsDF[[3]]$Cellch01==i)),
                  prop2Strict = length(which(compSpotsDF[[3]]$Colocalised2Strict==1&compSpotsDF[[3]]$Cellch01==i))/length(which(compSpotsDF[[3]]$Cellch01==i)),
                  prop2StrictLeniency = length(which(compSpotsDF[[3]]$Colocalised2StrictLeniency==1&compSpotsDF[[3]]$Cellch01==i))/length(which(compSpotsDF[[3]]$Cellch01==i)),
                  prop2xStrict = length(which(compSpotsDF[[3]]$Colocalised2xStrict==1&compSpotsDF[[3]]$Cellch01==i))/length(which(compSpotsDF[[3]]$Cellch01==i)),
                  prop2xStrictLeniency = length(which(compSpotsDF[[3]]$Colocalised2xStrictLeniency==1&compSpotsDF[[3]]$Cellch01==i))/length(which(compSpotsDF[[3]]$Cellch01==i)),
                  numFoci = length(which(compSpotsDF[[3]]$Cellch01==i)))
    proportions=rbind(proportions,df)
  }
  return(proportions)
}

IDFierMerger <- function(listOfDFs, listOfIDs=NULL) {
  print(listOfIDs)
  if(length(listOfIDs)<1){
    listOfIDs = seq(1,length(listOfDFs),1)
  }
  print(listOfIDs)
  for(i in 1:length(listOfDFs)){
    listOfDFs[i][[1]]$IDFier = listOfIDs[i]
  }
  return(rbindlist(listOfDFs))
}