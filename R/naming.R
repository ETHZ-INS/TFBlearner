# Feature matrix columns are named in the following way.

# <featureType>_<featureName> depending of the feature followed by (<motifSuffix>, <i>, <normedSuffix>)

# feature-types: ---------------------------------------------------------------

siteFeat <- "siteFeat"
contextFeat <- "contextFeat"
panContextFeat <- "panContextFeat"
tfFeat <- "tfFeat"
contextTfFeat <- "contextTfFeat"

# base-experiment-names --------------------------------------------------------
motifExp <- "Motifs"
atacExp <- "ATAC"
atacPromExp <- "ATAC.Promoters"
chIPExp <- "ChIP"
assocExp <- "Activity.Association"
actExp <- "chromVAR.Activity"

# various often used affixes ---------------------------------------------------

motifAffix <- "Motif"
pearsonAffix <- "Pearson"
cohenAffix <- "CohenKappa"
promoterAffix <- "promoter"
coMotifAffix <- "co"
exMotifAffix <- "ex"
overlapsAffix <- "overlaps"
insertsAffix <- "inserts"
normedAffix <- "normed"
gcNormedAffix <- "normedGC"
normedMaxAffix <- "normedMaxATAC"
predPrefix <- "pred"

# base-assay-names -------------------------------------------------------------
matchAssay <- "matchScores"
devAssay <- "deviations"
normDevAssay <- "normedDev"
peakAssay <- "peaks"
assocAssay <- paste0(tolower(pearsonAffix), "Activity")

# for fragment overlap & insert counts
# => see totalOverlapsName/overlapTypeFeatNames and totalInsertsName/insertTypeFeatNames

# col-/rowData names -----------------------------------------------------------

tfNameCol <- "tfName"
motifNameCol <- "motifName"
preSelMotifCol <- paste0("presel", motifAffix)
preSelActCol <- paste0("preselAct", motifAffix)
isTestCol <- "isTesting"
isTrainCol <- "isTraining"
featTypeCol <- "featureType"
tfCofactorsCol <- "tfCofactors"
topVarSitesCol <- "topVarSite"
chromVarExpCol <- "chromVarExpecations"
chromVarBgCol <- "chromVarBackgrounds"
maxScoreCol <- "maxScore"

# return names -----------------------------------------------------------------

retScoresName <- "motifInsertCounts"
reProfileName <- "insertProfiles"

# metaData names ---------------------------------------------------------------

# motif-suffixes
tfCofactorPrefix <- "tfCofactor"
tfCofactorMotifPrefix <- paste0(tfCofactorPrefix, motifAffix)
assocActPrefix <- paste0("assocation",motifAffix,"Activity")
tfMotifPrefix <- paste0("tf", motifAffix)
priorMotifPrefix <- paste0("prior", motifAffix) # motifs chosen based on matching names
selMotifPrefix <- paste0("selected", motifAffix)
ctcfMotifPrefix <- paste0("ctcf", motifAffix)

# feature-names ----------------------------------------------------------------
# the same logic applies in code, if more than the affices or feature names
# need to be changed, e.g. the order of the affices,
# then it would have to be also changed at the respective places in the code

# siteFeats
consScoreFeatName <- "conservationScores"
cpgDensFeatName <- "cpgDensity"
gcContFeatName <- "gcContent"
widthFeatName <- "width"

# tfFeats
promFeatName <- paste0(promoterAffix, c(pearsonAffix, cohenAffix))
motifFeatName <- paste0(tolower(motifAffix), "Match")
actAssocFeatName <- assocActPrefix
coBindFeatName <- "coBind"
coCountFeatName <- "coMotifCount"
cScoreFeatName <- "cScore"
ctcfFeatName <- "ctcfBind"
patContextFeatName <- "patternContext"
patTfFeatName <- "patternTf"

# panContextFeats
atacVarFeatName <- "atacVariance"
maxAtacFeatName <- "atacMax"

# contextFeats
mdsDimFeatName <- "mdsContext"
totalOverlapsFeatName <- paste("total", overlapsAffix, sep=".")
totalInsertsFeatName <-  paste("total", insertsAffix, sep=".")
typeNames <- c("nucleosomefree", "mononucleosome",
               "dinucleosome", "multinucleosome")
overlapTypeFeatNames <- paste(typeNames, overlapsAffix, sep=".")
insertTypeFeatNames <- paste(typeNames, insertsAffix, sep=".")

# contextTfFeats
insertFeatName <- "inserts"
wInsertsFeatName <- "weightedInserts"
devFeatName <- "chi2DevProfile"
chromVarFeatName <- "chromVAR.Activity"
labelName <- "label"

# feature-matrix column names --------------------------------------------------
# some feature matrix column names which get used over and over


# labelBinCol

# column-labels
labelColName <- paste(contextTfFeat, labelName, sep="_")
countColName <- paste(contextFeat, totalOverlapsFeatName, normedAffix, sep="_")
motifFeatColName <- paste(tfFeat, motifFeatName, tfMotifPrefix, 1, sep="_")
maxATACColName <- paste(panContextFeat, maxAtacFeatName, sep="_")
gcContentColName <- paste(siteFeat, gcContFeatName)
cScoreColName <- paste(tfFeat, cScoreFeatName, sep="_")
binLabelName <- "labelBin"

# model component names --------------------------------------------------------

modelTopWeightName <- "top_weighted_pos"
modelMedWeightName <- "med_weighted_pos"
modelAllWeigthName <- "all_weighted_pos"
modelAllName <- "all_pos"
modelStackedSuffix <- "stacked"

# stacking strategy names
stackModelLast <- "last"
stackModelWeightLast <- "wLast"
stackModelWeightMean <- "wMean"
stackModelBoostTree <- "boostTree"

stackingStratEntry <- "stacking_strategy"


