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


testthat::test_that("distJSD", {
  withr::local_package("sf")
  withr::local_options(list(sf_use_s2 = FALSE))

  ncpath <- system.file("gpkg/nc.gpkg", package = "sf")
  nc <- sf::read_sf(ncpath)
  ncn <- nc[, 9:14]
  ncn <- sf::st_drop_geometry(ncn)
  ncn <- as.matrix(ncn) + 1L

  jsd_ncn <- distJSD(ncn)

  testthat::expect_equal(jsd_ncn[1, 1], 0L, tolerance = 1e-8)
  testthat::expect_equal(jsd_ncn[2, 1], 218.7963, tolerance = 7e-5)
})


testthat::test_that("distJSD2", {
  withr::local_package("sf")
  withr::local_options(list(sf_use_s2 = FALSE))

  ncpath <- system.file("gpkg/nc.gpkg", package = "sf")
  nc <- sf::read_sf(ncpath)
  ncn <- nc[, 9:14]
  ncn <- sf::st_drop_geometry(ncn)
  ncn <- as.matrix(ncn) + 1L

  jsd_ncn <- distJSD2(ncn[1:3, ], ncn[60:67, ])

  testthat::expect_equal(dim(jsd_ncn), c(3, 8))
  testthat::expect_equal(jsd_ncn[1, 5], 3415.587, tolerance = 1e-3)
  testthat::expect_equal(jsd_ncn[3, 7], 687.8441, tolerance = 1e-4)

})

testthat::test_that("speedmat-legacy test suites", {
  withr::local_package("sf")
  withr::local_options(list(sf_use_s2 = FALSE))

  ncpath <- system.file("gpkg/nc.gpkg", package = "sf")
  nc <- sf::read_sf(ncpath)
  ncn <- nc[, 9:14]
  ncn <- sf::st_drop_geometry(ncn)
  ncn <- as.matrix(ncn) + 1L

  # error cases
  testthat::expect_error(
    speedmat_legacy(nc, "WHATEVER")
  )
  testthat::expect_error(
    speedmat_legacy(
      sf = nc,
      mode = "CDE"
    )
  )
  testthat::expect_error(
    speedmat_legacy(
      sf = nc,
      mode = "CDE"
    )
  )
})


testthat::test_that("speedmat test suites", {
  withr::local_package("sf")
  withr::local_options(list(sf_use_s2 = FALSE))

  ncpath <- system.file("gpkg/nc.gpkg", package = "sf")
  nc <- sf::read_sf(ncpath)
  ncn <- nc[, 9:14]
  ncc <- sf::st_centroid(nc)
  ncc <- cbind(ncc, data.frame(sf::st_coordinates(ncc)))
  ncn <- sf::st_drop_geometry(ncn)
  ncn <- as.matrix(ncn) + 1L

  # error cases
  testthat::expect_error(
    speedmat_legacy(
      sf = ncc,
      formula =
      reformulate(colnames(ncc)[1:5], colnames(ncc)[6]),
      mode = "CDE"
    )
  )

})