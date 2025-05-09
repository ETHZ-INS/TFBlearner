test_that("Shifting checks", {
  #TODO: move this test to the corresponding file
  assayTableSimple1$strand <- c("+", "-", "+", "-")
  atacFragShifted <- .processData(assayTableSimple1, shift=TRUE)
  expect_equal(atacFragShifted$start, assayTableSimple1$start+c(-4,-4,-4,-4))
  expect_equal(atacFragShifted$end, assayTableSimple1$end+c(-4,-4,-4,-4))
})
