testthat::test_that(
  "GitHub Action download to runner tempdir works", {
  tdir <- tempdir()
  downloadlink <- "https://geodacenter.github.io/data-and-lab/data/natregimes.zip"
  target_filename <- "natregimes.zip"
  target_fullname <- file.path(tdir, target_filename)
  download.file(downloadlink, target_fullname, "libcurl")
  testthat::expect_true(file.exists(target_fullname))
  }
)