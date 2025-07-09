# TFBlearner

 [![R-CMD-check](https://github.com/ETHZ-INS/TFBlearner/actions/workflows/r_cmd_check.yaml/badge.svg)](https://github.com/ETHZ-INS/TFBlearner/actions/workflows/r_cmd_check.yaml)
 [![codecov](https://codecov.io/gh/ETHZ-INS/TFBlearner/graph/badge.svg)](https://codecov.io/gh/ETHZ-INS/TFBlearner)
 
Package for computing features and training of TF-specific models for TF-binding predictions based on ATAC-seq data.
<br />
<br />

![*Package functionality overview*](./schemes/overview_small.png)   

<br />
<br />
     

 The package can be installed with: 

``` r
 devtools::install_github("https://github.com/ETHZ-INS/TFBlearner")
```

## Input data

As inputs for training a model motif-matches, ATAC-seq data, and ChIP-seq peaks (used as labels) are required.

``` r
library(TFBlearner)
library(BSgenome.Hsapiens.UCSC.hg38)
library(phastCons100way.UCSC.hg38)

# load example data
# genomic coordinates of sites
data("example_coords")

# atac-seq data
exampleATAC <- list(A549=system.file("extdata", "example_atac_A549.bed", package = "TFBlearner"),
                    K562=system.file("extdata", "example_atac_K562.bed", package = "TFBlearner"))

# chip-seq data
exampleChIP <- list(K562_CTCF=system.file("extdata", "example_chIP_K562_ctcf.tsv", package = "TFBlearner"),
                    A549_CTCF=system.file("extdata", "example_chIP_A549_ctcf.tsv", package = "TFBlearner"),
                    K562_JUN=system.file("extdata", "example_chIP_K562_jun.tsv", package = "TFBlearner"))
```

In a first step motif-matches can be obtained and prepared in the required format by calling `prepMotifs()`.

``` r
data("example_pwms")
exampleMotif <- prepMotifs(example_coords, example_pwms,
                           genome=BSgenome.Hsapiens.UCSC.hg38,
                           outputFolder=".")
```

The arguments can be passed as named lists, with names corresponding to the cellular contexts of the data to `prepData()` which outputs a custom [MultiAssayExperiment](https://www.bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html) object with the different modalities in the `experiments()`.
Cellular contexts used for testing / validation can already be defined at this point with the argument `testSet` and the ChIP-seq peaks of these will be excluded from all downstream feature construction steps.

``` r
mae <- prepData(example_coords,
                motifData=exampleMotif,
                atacData=exampleATAC,
                chIPData=exampleChIP,
                testSet="A549")
```

## Feature construction

The features compiled can be split in fives categories:            
1. site-specific       
2. context-specific          
3. pan-context        
4. TF-specific   
5. cellular context and TF-specific       

An overview of the features can be obtained with: 

``` r
listFeatures()
```

Features are added to the provided MultiAssayExperiment objects as experiments. 

``` r
# site-specific features
mae <- siteFeatures(mae, 
                    phast=phastCons100way.UCSC.hg38,
                    genome=BSgenome.Hsapiens.UCSC.hg38)
                    
# context-specific & pan-context features
mae <- panContextFeatures(mae, features=c("Max_ATAC_Signal"))

# TF-specific features (JUN used here for examplary purposes as a cofactor of CTCF)
mae <- tfFeatures(mae, tfName="CTCF", tfCofactors="JUN")

# cellular context & TF-specific features
mae <- contextTfFeatures(mae, tfName="CTCF", addLabels=TRUE)
```

## Training TF-specific models

Training is split up in two steps, first training of single gradient-boosted tree models using the [lightgbm library](https://lightgbm.readthedocs.io/en/stable/R/reference/) on different  strata of the data (different stringency on whats included as a ChIP-seq peak), second training of a stacking model combining the predictions of the single models.
Both steps are performed when calling `trainTfModel()`.

``` r
# get feature matrix for training
fmTrain <- getFeatureMatrix(mae, tfName="CTCF", whichCol="OnlyTrain")

# train models
mods <- trainTfModel(fm=fmTrain, tfName="CTCF", evalRounds=5, stackingStrat="wMean")
```

## Predicting binding

In a last step binding predictions can be obtained with `predictTfBinding()`.

``` r
# get feature matrix on which to predict
fmVal <- getFeatureMatrix(mae, tfName="CTCF", whichCol="Col", colSel="A549")

# predict bindings
predsValBagged <- predictTfBinding(mods, fmVal)
```

