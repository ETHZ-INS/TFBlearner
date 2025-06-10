test_that("It works :)", {
  expect_no_error(listFeatures())
})

test_that("All feature name options in args of functions are listed :)", {
  siteFeats <- eval(formals(TFBlearner::siteFeatures)$features)
  tfFeats <- eval(formals(TFBlearner::tfFeatures)$features)
  conFeats <- eval(formals(TFBlearner::panContextFeatures)$features)
  conTfFeats <- eval(formals(TFBlearner::contextTfFeatures)$features)
  argFeats <- c(siteFeats, tfFeats, conFeats, conTfFeats)

  lDt <- listFeatures()
  expect_in(argFeats, lDt$feature_name)
})
