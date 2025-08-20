
# Run these once from the repo root -----------------------------

# 0) Tools you’ll need on this machine
install.packages(c(
  "usethis","devtools","roxygen2","testthat","Rcpp",
  "signal","zoo","pkgdown","covr","lintr",
  "tibble","dplyr","ggplot2","tidyr"
))

# 1) Turn the cloned repo into a proper R package
usethis::create_package(".")   # creates DESCRIPTION, NAMESPACE scaffolding, R/, man/, etc.
usethis::use_rstudio()         # (optional) creates an .Rproj

# 2) Initialize git metadata if needed (you’re using GitHub Desktop already)
usethis::use_git()

# 3) Connect to GitHub if not already connected by Desktop (skip if it is)
# usethis::use_github(protocol = "https", private = FALSE)

# 4) Licensing + DESCRIPTION clean up
usethis::use_mit_license("Victor Rwandarwacu")
usethis::use_tidy_description()

# 5) Rcpp infrastructure
usethis::use_rcpp()
Rcpp::compileAttributes()
devtools::document()           # generates NAMESPACE + man/ from roxygen
