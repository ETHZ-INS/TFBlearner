test_that("Saving & loading Hdf5 mae", {

  outDir <- tempdir()
  # save model
  maeFilePath <- file.path(outDir, "maeTestHdf5.rds")
  expect_no_error(saveMae(maeTestHdf5, filepath=maeFilePath))
  expect_true(file.exists(maeFilePath))

  # load model
  maeLoaded <- NULL
  expect_no_error(maeLoaded <- loadMae(maeFilePath))
  expect_equal(names(experiments(maeTestHdf5)),
               names(experiments(maeLoaded)))
  file.remove(maeFilePath)
})

test_that("Change the paths of the HDF5 files", {
  outDir <- tempdir()
  file.copy("./test_data/ATAC_mapped.h5", to=file.path(outDir,
                                                       "ATAC_mapped.h5"))
  file.copy("./test_data/ChIP_mapped.h5", to=file.path(outDir,
                                                       "ChIP_mapped.h5"))
  file.copy("./test_data/Motif_mapped.h5", to=file.path(outDir,
                                                       "Motif_mapped.h5"))
  maeRebased <- NULL
  expect_no_error(maeRebased <- rebaseMaeH5paths(maeTestHdf5, newBase=outDir))
  assayMat <- assays(experiments(maeRebased)[[atacExp]])[[totalOverlapsFeatName]]
  expect_equal(dirname(assayMat@seed@seed@filepath), outDir)

  file.remove(file.path(outDir, "ATAC_mapped.h5"))
  file.remove(file.path(outDir, "ChIP_mapped.h5"))
  file.remove(file.path(outDir, "Motif_mapped.h5"))
})
