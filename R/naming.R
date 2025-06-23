# Feature matrix columns are named in the following way.

# <featureType>_<featureName> depending of the feature followed by (<motifSuffix>, <i>, <normedSuffix>)

# feature-types: ---------------------------------------------------------------

SITEFEAT <- "siteFeat"
CONTEXTFEAT <- "contextFeat"
PANCONTEXTFEAT <- "panContextFeat"
TFFEAT <- "tfFeat"
CONTEXTTFFEAT <- "contextTfFeat"

# base-experiment-names --------------------------------------------------------
MOTIFEXP <- "Motifs"
ATACEXP <- "ATAC"
ATACPROMEXP <- "ATAC.Promoters"
CHIPEXP <- "ChIP"
ASSOCEXP <- "Activity.Association"
ACTEXP <- "chromVAR.Activity"

# various often used affixes ---------------------------------------------------

MOTIFAFFIX <- "Motif"
PEARSONAFFIX <- "Pearson"
COHENAFFIX <- "CohenKappa"
PROMOTERAFFIX <- "promoter"
COMOTIFAFFIX <- "co"
EXMOTIFAFFIX <- "ex"
OVERLAPSAFFIX <- "overlaps"
INSERTSAFFIX <- "inserts"
NORMEDAFFIX <- "normed"
GCNORMEDAFFIX <- "normedGC"
NORMEDMAXAFFIX <- "normedMaxATAC"
PREDPREFIX <- "pred"

# base-assay-names -------------------------------------------------------------
MATCHASSAY <- "matchScores"
DEVASSAY <- "deviations"
NORMDEVASSAY <- "normedDev"
PEAKASSAY <- "peaks"
ASSOCASSAY <- paste0(tolower(PEARSONAFFIX), "Activity")

# for fragment overlap & insert counts
# => see totalOverlapsName/OVERLAPTYPEFEATNAMES and totalInsertsName/INSERTTYPEFEATNAMES

# col-/rowData names -----------------------------------------------------------

TFNAMECOL <- "tfName"
MOTIFNAMECOL <- "motifName"
PRESELMOTIFCOL <- paste0("presel", MOTIFAFFIX)
PRESELACTCOL <- paste0("preselAct", MOTIFAFFIX)
ISTESTCOL <- "isTesting"
ISTRAINCOL <- "isTraining"
FEATTYPECOL <- "featureType"
TFCOFACTORSCOL <- "tfCofactors"
TOPVARSITESCOL <- "topVarSite"
CHROMVAREXPCOL <- "chromVarExpecations"
CHROMVARBGCOL <- "chromVarBackgrounds"
MAXSCORECOL <- "maxScore"
SCORECOL <- "score"
COOCCURRENCECOL <- "cooccurrences"
MDSSUBROWCOL <- "mdsSubRows"
BASEDIRCOL <- "baseDir"

# return names -----------------------------------------------------------------

RETSCORESNAME <- "motifInsertCounts"
REPROFILENAME <- "insertProfiles"

# metaData names ---------------------------------------------------------------

MDSDIMSTATSENTRY <- "MDSRes"

# motif-suffixes ---------------------------------------------------------------
TFCOFACTORPREFIX <- "tfCofactor"
TFCOFACTORMOTIFPREFIX <- paste0(TFCOFACTORPREFIX, MOTIFAFFIX)
ASSOCACTPREFIX <- paste0("assocation",MOTIFAFFIX,"Activity")
TFMOTIFPREFIX <- paste0("tf", MOTIFAFFIX)
PRIORMOTIFPREFIX <- paste0("prior", MOTIFAFFIX) # motifs chosen based on matching names
SELMOTIFPREFIX <- paste0("selected", MOTIFAFFIX)
CTCFMOTIFPREFIX <- paste0("ctcf", MOTIFAFFIX)

# feature-names ----------------------------------------------------------------
# the same logic applies in code, if more than the affices or feature names
# need to be changed, e.g. the order of the affices,
# then it would have to be also changed at the respective places in the code

# siteFeats
CONSSCOREFEATNAME <- "conservationScores"
CPGDENSEATNAME <- "cpgDensity"
GCCONTFEATNAME <- "gcContent"
WIDTHFEATNAME <- "width"

# tfFeats
PROMFEATNAME <- paste0(PROMOTERAFFIX, c(PEARSONAFFIX, COHENAFFIX))
MOTIFFEATNAME <- paste0(tolower(MOTIFAFFIX), "Match")
ACTASSOCFEATNAME <- ASSOCACTPREFIX
COBINDFEATNAME <- "coBind"
COCOUNTFEATNAME <- "coMotifCount"
CSCOREFEATNAME <- "cScore"
CTCFFEATNAME <- "ctcfBind"
PATCONTEXTFEATNAME <- "patternContext"
PATTFFEATNAME <- "patternTf"

# panContextFeats
ATACVARFEATNAME <- "atacVariance"
MAXATACFEATNAME <- "atacMax"

# contextFeats
MDSDIMFEATNAME <- "mdsContext"
TOTALOVERLAPSFEATNAME <- paste("total", OVERLAPSAFFIX, sep=".")
TOTALINSERTSFEATNAME <-  paste("total", INSERTSAFFIX, sep=".")
TYPENAMES <- c("nucleosomefree", "mononucleosome",
               "dinucleosome", "multinucleosome")
OVERLAPTYPEFEATNAMES <- paste(TYPENAMES, OVERLAPSAFFIX, sep=".")
INSERTTYPEFEATNAMES <- paste(TYPENAMES, INSERTSAFFIX, sep=".")

# contextTfFeats
INSERTFEATNAME <- "inserts"
WINSERTSFEATNAME <- "weightedInserts"
DEVFEATNAME <- "chi2DevProfile"
CHROMVARFEATNAME <- "chromVAR.Activity"
LABELNAME <- "label"

# feature-matrix column names --------------------------------------------------
# some feature matrix column names which get used over and over


# labelBinCol

# column-labels
LABELCOLNAME <- paste(CONTEXTTFFEAT, LABELNAME, sep="_")
COUNTCOLNAME <- paste(CONTEXTFEAT, TOTALOVERLAPSFEATNAME, NORMEDAFFIX, sep="_")
MOTIFFEATCOLNAME <- paste(TFFEAT, MOTIFFEATNAME, TFMOTIFPREFIX, 1, sep="_")
MAXATACCOLNAME <- paste(PANCONTEXTFEAT, MAXATACFEATNAME, sep="_")
GCCONTENTCOLNAME <- paste(SITEFEAT, GCCONTFEATNAME)
CSCORECOLNAME <- paste(TFFEAT, CSCOREFEATNAME, sep="_")
BINLABELNAME <- "labelBin"

# model component names --------------------------------------------------------

MODELTOPWEIGHTNAME <- "top_weighted_pos"
MODELMEDWEIGHTNAME <- "med_weighted_pos"
MODELALLWEIGHTNAME <- "all_weighted_pos"
MODELALLNAME <- "all_pos"
MODELNAMES <- c(MODELTOPWEIGHTNAME, MODELMEDWEIGHTNAME,
                MODELALLWEIGHTNAME, MODELALLNAME)
MODELSTACKEDSUFFIX <- "stacked"

# stacking strategy names
STACKMODELLAST <- "last"
STACKMODELWEIGHTLAST <- "wLast"
STACKMODELWEIGHTMEAN <- "wMean"
STACKMODELBOOSTTREE <- "boostTree"

STACKINGSTRATENTRY <- "stacking_strategy"

PACKAGEVERSION <- "package_version"
SPARSETHR <- "sparse_thr"
