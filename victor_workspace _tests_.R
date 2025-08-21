# Synthetic
n <- 5000; fs <- 200
t <- seq(0, 25, length.out = n)
abp <- 90 + 25*pmax(0, sin(2*pi*1.2*t)) + rnorm(n,0,0.5)
ecg <- sin(2*pi*1.2*t) + rnorm(n,0,0.1)
ppg <- pmax(0, sin(2*pi*1.2*t)) + 0.1*rnorm(n)

head(extractFeaturesEnhanced(abp, ppg, ecg, fs))

#test on csv file  (work directory R folder )

untreated <- process_condition("untreated_signals.csv")
treated   <- process_condition("treated_signals.csv")

summary(untreated)
summary(treated)

devtools::document()
devtools::load_all()
untreated <- process_condition(
  system.file("extdata", "untreated_signals.csv", package = "abpwaveforms")
)
treated <- process_condition(
  system.file("extdata", "treated_signals.csv", package = "abpwaveforms")
)
summary(untreated)
summary(untreated)




#create test
usethis::use_testthat(3)               # sets up tests/ + testthat edition 3
usethis::use_test("extractFeaturesEnhanced")  # creates tests/testthat/test-extractFeaturesEnhanced.R

#full check

devtools::document()
devtools::test()
devtools::check()
