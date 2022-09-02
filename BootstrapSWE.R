# generate paired datasets for SWE from partial duration analysis
require(dssrip)
require(lubridate)
require(ggplot2)
# require(hydroutils) # github.com/eheisman/hydroutils
require(reshape2)
require(plyr)
require(rJava)

nEvents = 250
primaryBasin <- "MNTREGION"



sweFile <- opendss("./SWE_UA_basinAve.dss")
records <- pathsToDataFrame(getPaths(sweFile, "F=UA-POUDRE", searchFunction=pathByPartsWildcard))
subbasins <- records$LOCATION


pathTemplate <- "//%s/SWE//1DAY/UA-POUDRE/"

getRecordForBasin <- function(basin){
  dataset = fortify.zoo(getTSC(sweFile, sprintf(pathTemplate, basin),fullTSC=TRUE))
  dataset = data.frame(year=year(dataset$Index), month=month(dataset$Index), day=day(dataset$Index), yday=yday(dataset$Index), 
                       swe=round(dataset$SWE, 1))
  dataset = subset(dataset, month %in% c(3:6)) # optionally | (month == 7 & day < 15))
}

durationBasin = getRecordForBasin(primaryBasin)

ggplot(durationBasin) + geom_line(aes(x=yday, y=swe, group=year))

SWE_probs = quantile(durationBasin$swe, probs=c(1:nEvents)/(nEvents+1), type=3) 
# grab quantiles by weibull position
# use type 1 to get inverse of emperical distribution, ensure we always have a result to match in record. (e.g. no averaging)
#subset(durationBasin, swe == SWE_probs[1] )

allBasins = ldply(subbasins, function(sb){
  dataset=getRecordForBasin(sb)
  dataset$location = sb
  return(dataset)
})

mergedBasins = dcast(allBasins, year + yday ~ location, value.var="swe")

# make this repeatable per watershed
set.seed(ncol(mergedBasins))

sweSamples = ldply(SWE_probs, function(swe) {
  quantile_subset = subset(mergedBasins, MNTREGION==swe)
  index = sample(1:nrow(quantile_subset),1)
  quantile_subset[index,]
})

# resample to get a random order
resampleIndex = sample(1:nrow(sweSamples), replace=F, size=nrow(sweSamples))
outSWE = sweSamples[resampleIndex,]
plot(outSWE$MNTREGION) # these should look adequately randomized
plot(outSWE$S_CASCADECK)
plot(outSWE$S_SHEEPCK)

outDSS = opendss("./parameter_pairs_SWE.dss")
writeSample <- function(basin, values){
  outPDC = .jnew("hec.io.PairedDataContainer")
  #outPDC$setNumberCurves(as.integer(1))
  #outPDC$setNumberOrdinates(length(values))
  outPDC$setValues(.jarray(as.numeric(1:length(values))),
                   .jarray(t(as.matrix(values, ncol=1)), dispatch=TRUE))
  outPDC$fullName = sprintf("//%s/INITIAL SWE///SAMPLED/", basin)
  outDSS$put(outPDC)
  #return(outPDC)
}

l_ply(subbasins, function(sb){
  writeSample(sb, outSWE[,sb])
})

write.csv(outSWE, "output_SWE_pairs.csv")

