# In order to have one point to easily change the naming of assays & experiments

# general column-names & prefix ------------------------------------------------

tfNameCol <- "tf_name"

# cofactor / motif
tfCofactorsCol <- "tf_cofactors"
preSelMotifCol <- "preselected_motifs"

isTestCol <- "is_testing"
isTrainCol <- "is_training"
featTypeCol <- "feature_type"

# prefices
promoterPrefix <- "promoter"
tfCofactorPrefix <- "tfCofactor"
tfMotifPrefix <- "tf_motif"
assocMotifPrefix <- "associated_motif"

# ATAC-experiment --------------------------------------------------------------

# experiment-name:
atacExp <- "ATAC"
atacPromeExp <- "ATAC_promoters"

# assay-names / prefices:
# typeNames are defined via argument of .getFragType
totalOverlapsName <- "total_overlaps"
normTotalOverlapsName <- "norm_total_overlaps"
totalInsertsName <- "total_inserts"

typeOverlapSuffix <- "overlaps"
typeInsertsSuffix <- "inserts"

# ChIP-experiment --------------------------------------------------------------

# experiment-name:
chIPExp <- "ChIP"

# assay-names:
peakAssayName <- "peaks"

# Motif-experiment -------------------------------------------------------------

# experiment-name:
motifExp <- "Motifs"

# assay-names:
matchAssayName <- "match_scores"

# coldata-names:
maxScoreCol <- "max_score"
motifNameCol <- "motif"

# will be removed
matchRangesExp <- "match_ranges"

# Activity-experiment ----------------------------------------------------------

actExp <- "Activity"

# siteFeat-experiment ----------------------------------------------------------

# experiment-name:
siteFeat <- "siteFeat"

# assay-names:
consScoreFeatName <- "conservation_scores"
cpgDensFeatName <- "cpg_density"
gcContFeatName <- "gc_content"
widthFeatName <- "width"

# coldata-names:
chromVarExpName <- "ChromVAR_expectations"
chromVarBgName <- "ChromVAR_background_peaks"
topVarSitesName <- "top_var_sites"

# tfFeat-experiment ------------------------------------------------------------

# experiment-name:
tfFeat <- "tfFeat"

# assay-names / prefices:
cScoreFeatName <- "c_score"
ctcfSigFeatName <- "CTCF_signal"
motifMatchesFeatName <- "motif_matches" # associated motif-matching scores

tfPatternPrefix <- "pattern_tf"
contextPatternPrefix <- "pattern_context"
assocCohenPrefix <- "Cohen_Kappa"
assocPearsonPrefix <- "Pearson"

cofactorBindingPrefix <- "cofactor_binding"
coTfMotifPrefix <- "n_cooccuring_motifs_tf"
coCofactorMotifPrefix <- "n_cooccuring_motifs_tfCofactor"

# contextTfFeat-experiment -----------------------------------------------------

# experiment-name

contextTfFeat <- "contextTfFeat"

# assay-names /prefices:
insertFeatName <- "insert_counts"
wInsertsFeatName <- "weighted_insert_counts"
devFeatName <- "chi2_dev_profile"
mdsDimFeatName <- "MDS_Context"
atacVarFeatName <- "ATAC_Variance"
maxAtacFeatName <- "Max_ATAC_Signal"

chromVarAssocSuffix <-  "ChromVAR_ATAC"
chromVarScoreSuffix <- "ChromVAR_score"
labelFeatName <- "label"

maxAtacScaledPrefix <- "maxATAC_scaled"

# getInsertionProfiles return names --------------------------------------------

retScoresName <- "motifScores"
reProfileName <- "profile"


# Training ---------------------------------------------------------------------

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

# Prediction -------------------------------------------------------------------

binLabelName <- "label_bin"
