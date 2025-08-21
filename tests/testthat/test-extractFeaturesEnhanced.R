test_that("extractFeaturesEnhanced runs and returns expected columns", {
  n <- 3000; fs <- 200
  t <- seq(0, 15, length.out = n)
  abp <- 90 + 20*pmax(0, sin(2*pi*1.1*t))
  ecg <- sin(2*pi*1.1*t)
  ppg <- pmax(0, sin(2*pi*1.1*t))

  out <- extractFeaturesEnhanced(abp, ppg, ecg, fs)

  expect_true(is.data.frame(out))
  expect_true(all(c("pulse_pressure","AIx","RWTT","PAT","HR","ppg_width_half_height") %in% names(out)))
  expect_false(any(is.infinite(out$AIx), na.rm = TRUE))

  med_hr <- stats::median(out$HR, na.rm = TRUE)
  expect_true(is.finite(med_hr) && med_hr > 30 && med_hr < 220)
})

