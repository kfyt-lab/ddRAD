library(data.table)
library(foreach)
library(doSNOW)
library(sjstats)
library(doParallel)
library(genetics)

rawDir ="/media/hrivnak/SSD/Hrivnak/PiceaRAD20" # Raw data directory.
setwd("/media/hrivnak/SSD/Hrivnak/PiceaRAD20/Base")

# ----- Filter catalog.calls (to be small enough to handle) -----
iFile = file(paste0(rawDir, "/stacks/catalog.calls"), "r") # STACKS output. Nucleotide summaries.
oFile = file("catalog.calls.rC2", "w")
app   = F
while (T) {
	line = readLines(iFile, n=1)
	if (length(line) == 0) {
		break }
	line = strsplit(line, "\t")[[1]]
	if ((substr(line[1], 1, 1) == "#") | (line[6] != ".")) {
		cat(paste(line[1:min(length(line), 8)], collapse="\t"), "\n", sep="", file=oFile, append=app)
		app=T } }
close(iFile)
close(oFile)

# ----- Filter the STRUCTURE data -----

# - Load data -
structure    = as.data.frame(fread(paste0(rawDir, "/stacks/populations.structure"), sep="\t", header=F)) # STACKS output. Input file for STRUCTURE, to be a general source of genotypes.
report       = as.data.frame(fread(paste0(rawDir, "/stacks/Sequencing_report.txt"), sep="\t", header=T)) # STACKS output. Sample info, number of reads etc.
callCat      = as.data.frame(fread("catalog.calls.rC2", sep="\t", header=F)) # Nucleotide summaries.

populations  = as.data.frame(fread("StructurePops.csv", sep="\t", header=T)) # Population numbers for samples.
sampleCheck  = rep(NA, nrow(populations))
for (i in 1:length(sampleCheck)) {
	sampleCheck[i] = grep(populations$Code1[i], structure[i+1, 1]) }
sampleCheck
all(sampleCheck==1)

# - Remove markers with Q<13 (p<95%) -
strMarkers  = c(structure[1, 3:ncol(structure)], recursive=T)
strMarkersQ = rep(-1, length(strMarkers))
for (i in 1:length(strMarkers)) {
	marker = as.integer(strsplit(strMarkers[i], "_")[[1]])
	strMarkersQ[i] = callCat$V6[(callCat$V1 == marker[1]) & (callCat$V2 == marker[2]+1)] }
strMarkersQ = as.numeric(strMarkersQ)
mask = strMarkers[strMarkersQ >= 13]
summary(strMarkersQ)
hist(strMarkersQ)
1-sum(mask)/length(mask)
structureGD1 = as.data.frame(structure)
structureGD1 = structureGD1[, c(T, T, strMarkers %in% mask)]

# - Remove markers with >20% NAs -
markerNA = colMeans(structureGD1 == "0")
markerNA = markerNA[3:length(markerNA)]
mask = markerNA <= .2
summary(mask)
1-sum(mask)/length(mask)

structureGD1 = structureGD1[, c(T, T, mask)]
cbind(structureGD1[,1], c(NA, populations$Code1))
structureGD1[, 2] = c(NA, populations$Pop)
fwrite(structureGD1, "populations.structureGD1", sep="\t", col.names=F)

# - Remove markers with <5% minor allelle frequency -
majorF = function(x) {
	t1 = table(x)
	t1 = sort(t1[c("1","2","3","4")], decreasing=T)
	return(t1[1]/sum(t1)) }
frequencies = apply(structureGD1[2:nrow(structureGD1), 3:ncol(structureGD1)], 2, majorF)
summary(frequencies)
mask = frequencies < .95
mean(mask)
summary(mask)

structureGD2 = structureGD1[, c(T, T, mask)]
fwrite(structureGD2, "populations.structureGD2", sep="\t", col.names=F)

# - Remove markers with +-2SD read depth -
strMarkers   = c(structureGD2[1, 3:ncol(structureGD2)], recursive=T)
strMarkersDP = rep(-1, length(strMarkers))
for (i in 1:length(strMarkers)) {
	marker = as.integer(strsplit(strMarkers[i], "_")[[1]])
	info   = callCat$V8[(callCat$V1 == marker[1]) & (callCat$V2 == marker[2]+1)]
	info   = strsplit(info, ";")[[1]]
	DP     = info[grep("DP=", info)]
	DP     = gsub("DP=", "", DP)
	strMarkersDP[i] = as.numeric(DP) }
summary(strMarkersDP)
strMarkersDPnorm = (strMarkersDP - mean(strMarkersDP)) / sd(strMarkersDP)
summary(strMarkersDPnorm)
mask = (strMarkersDPnorm >= -2) & (strMarkersDPnorm <= 2)
summary(mask)
summary(strMarkersDPnorm[mask])

structureGD3 = structureGD2[, c(T, T, mask)]
fwrite(structureGD3, "populations.structureGD3", sep="\t", col.names=F)

# - Check NA frequency in samples -
getNAs = function(x) {
	return(mean(x == "0")) }
NAs = apply(structureGD3[2:nrow(structureGD3), 3:ncol(structureGD3)], 1, getNAs)
hist(NAs)
summary(NAs)

# - Check sample number of reads -
rMean = mean(report$DemultiplexedReads)
rSD   = sd(report$DemultiplexedReads)
hist((report$DemultiplexedReads - rMean) / rSD)
summary((report$DemultiplexedReads - rMean) / rSD)

# ----- Filter correlations from the STRUCTURE data -----

structureGD4 = structureGD3
nCores = 12
registerDoParallel(cores=nCores)
strMarkers = c(structureGD4[1, 3:ncol(structureGD4)], recursive=T)
strNRow    = nrow(structureGD4)

toMajorMinor = function(x) {
	t1 = table(x)
	t1 = sort(t1[c("1","2","3","4")], decreasing=T)
	x[x == "0"] = NA
	x[x == names(t1)[1]] = "ma"
	x[x == names(t1)[2]] = "mi"
	if (length(table(x)) > 2) { stop("Error: length(table(x)) > 2") }
	x[x == "ma"] = 1
	x[x == "mi"] = 0
	return(as.numeric(x)) }

goods     = rep(NA, length(strMarkers))
pb        = txtProgressBar(min=1, initial=1, max=length(strMarkers), style=3)
goods[1]  = strMarkers[1]
nGoods    = 1
start     = Sys.time()
for (i in 2:length(strMarkers)) {
	good = T
	data1 = toMajorMinor(c(structureGD4[2:strNRow, c(F, F, strMarkers == strMarkers[i])], recursive=T))
	data1 = data1[c(T,F)] + data1[c(F,T)]
	for (j0 in floor(nGoods/nCores):0) {
		js = j0 * nCores + nCores:1
		js = js[js <= nGoods]
		if (length(js > 0)) {
			correlations = foreach(j=js, .inorder=T, .combine=c) %dopar% {
				data2 = toMajorMinor(c(structureGD4[2:strNRow, c(F, F, strMarkers == goods[j])], recursive=T))
				data2 = data2[c(T,F)] + data2[c(F,T)]
				mean(data1 == data2, na.rm=T) }
			if (any(correlations >= .8)) {
				good = F
				break } } }
	if (good) {
		nGoods = nGoods + 1
		goods[nGoods] = strMarkers[i] } 
	setTxtProgressBar(pb, i) }
stop = Sys.time()
stop - start
length(goods[!is.na(goods)])

structureGD4 = structureGD4[, c(T, T, strMarkers %in% goods)]
fwrite(structureGD4, "populations.structureGD4", sep="\t", col.names=F)

# ----- Filter correlations from the STRUCTURE data (using LD) -----

toSlashes = function(x) {
	x = paste(x[c(T,F)], x[c(F,T)], sep="/")
	x[x == "N/N"] = NA
	return(genotype(x)) }

structureGD4 = as.data.frame(fread("populations.structureGD3", sep="\t", header=F))
structureGD4[structureGD4 == "0"] = "N"
structureGD4[structureGD4 == "1"] = "A"
structureGD4[structureGD4 == "2"] = "C"
structureGD4[structureGD4 == "3"] = "G"
structureGD4[structureGD4 == "4"] = "T"

nCores = 12
registerDoParallel(cores=nCores)
strMarkers = c(structureGD4[1, 3:ncol(structureGD4)], recursive=T)
strNRow    = nrow(structureGD4)

goods     = rep(NA, length(strMarkers))
goods[1]  = strMarkers[1]
nGoods    = 1
pb        = txtProgressBar(min=1, initial=1, max=length(strMarkers), style=3)
start     = Sys.time()
#for (i in 2:100) {
for (i in 2:length(strMarkers)) {
	good = T
	data1 = toSlashes(c(structureGD4[2:strNRow, c(F, F, strMarkers == strMarkers[i])], recursive=T))
	for (j0 in floor(nGoods/nCores):0) {
		js = j0 * nCores + nCores:1
		js = js[js <= nGoods]
		if (length(js > 0)) {
			correlations = foreach(j=js, .inorder=T, .combine=c) %dopar% {
				data2 = toSlashes(c(structureGD4[2:strNRow, c(F, F, strMarkers == goods[j])], recursive=T))
				LD(data1, data2)$D }
			if (any(abs(correlations) > .2)) {
				good = F
				break } } }
	if (good) {
		nGoods = nGoods + 1
		goods[nGoods] = strMarkers[i] } 
	setTxtProgressBar(pb, i) }
stop = Sys.time()
stop - start
length(goods[!is.na(goods)])

structureGD4 = as.data.frame(fread("populations.structureGD3", sep="\t", header=F))
structureGD4 = structureGD4[, c(T, T, strMarkers %in% goods)]
fwrite(structureGD4, "populations.structureGD4B", sep="\t", col.names=F)

# --- Remove differentiated individuals ---
structure = as.data.frame(fread("populations.structureGD3", sep="\t", header=F))

toRemove = c("15-P1470-G02", "54-P1462-F07", "47-P1474-G06", "78-P1466-F10", "88-P1528-H11", "29-P1443-E04", "69-P1449-E09", 
             "20-P1429-D03", "95-P1494-G12", "79-P1491-G10", "6-P1454-F01", "13-P1441-E02", "80-P1527-H10", "96-P1529-H12", 
             "68-P1436-D09", "53-P1447-E07", "61-P1448-E08")
toRemove[!(toRemove %in% structure$V1)] # Check if names are correct
length(toRemove) == 17 # Check if all individuals will be removed

structurens = structure[!(structure$V1 %in% toRemove),]
fwrite(structurens, "populations.structureGD3ns", sep="\t", col.names=F)
