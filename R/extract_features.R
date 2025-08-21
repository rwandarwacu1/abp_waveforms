#' Enhanced ABP/PPG/ECG feature extraction (200 Hz pipeline)
#' @param abp Numeric vector (ABP, mmHg) at ~1000 Hz; resampled to 200 Hz internally.
#' @param ppg Numeric vector (PPG).
#' @param ecg Numeric vector (ECG).
#' @param fs Numeric. Target sampling rate in Hz (default 200).
#' @return data.frame of per-beat features.
#' @export
extractFeaturesEnhanced <- function(abp, ppg, ecg, fs = 200) {
  res <- BP_resample(abp, origFs = 1000)
  signal200 <- stats::na.omit(res$newWaveform)

  dd <- doubleDerive(signal200)
  zoi <- rep(FALSE, length(signal200)); zoi[seq(100, length(zoi), by = 200)] <- TRUE

  foot <- FixIndex(getFootIndex(dd$waveformDDPlus, zoi))
  peak <- getPeakIndex(signal200, foot)
  dic  <- getDicroticNotchAndPeak(signal200, foot, peak, fs)
  notch <- dic$notch; dicrotic <- dic$dicrotic

  n <- min(length(foot), length(peak), length(notch), length(dicrotic)) - 1L
  foot <- foot[1:(n+1L)]; peak <- peak[1:n]; notch <- notch[1:n]; dicrotic <- dicrotic[1:n]

  pulse_pressure <- signal200[peak] - signal200[foot[1:n]]
  AIx <- (signal200[peak] - signal200[notch]) / pulse_pressure * 100
  RWTT <- (notch - foot[1:n]) / fs
  rise_time <- (peak - foot[1:n]) / fs
  pulse_duration <- (foot[2:(n+1L)] - foot[1:n]) / fs
  HR <- 60 / pulse_duration

  ecg_norm <- (ecg - mean(ecg, na.rm = TRUE)) / stats::sd(ecg, na.rm = TRUE)
  bf <- signal::butter(2, c(5, 15)/(fs/2), type = "pass")
  ecg_filt <- signal::filtfilt(bf, ecg_norm)
  ecg_diff <- diff(ecg_filt); ecg_sq <- ecg_diff^2
  win <- max(1L, round(0.150 * fs))
  ecg_mwi <- stats::filter(ecg_sq, rep(1/win, win), sides = 1)
  thr <- mean(ecg_mwi, na.rm = TRUE) + 0.5 * stats::sd(ecg_mwi, na.rm = TRUE)
  r_peaks <- which(ecg_mwi > thr); r_peaks <- r_peaks[c(TRUE, diff(r_peaks) > (0.25 * fs))]

  PAT <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    prev_r <- suppressWarnings(max(r_peaks[r_peaks < foot[i]]))
    if (is.finite(prev_r)) PAT[i] <- (foot[i] - prev_r) / fs
  }

  valid_range <- ppg >= 0 & ppg <= 5
  dppg <- c(0, diff(ppg)); valid_grad <- abs(dppg) < 0.5
  ppg[!(valid_range & valid_grad)] <- NA
  ppg <- zoo::na.approx(ppg, na.rm = FALSE)

  ppg_rise_time <- (peak - foot[1:n]) / fs
  ppg_sys <- ppg[peak]; ppg_dia <- ppg[foot[1:n]]
  s_d_ratio <- ppg_sys / ppg_dia
  dicrotic_timing <- (dicrotic - foot[1:n]) / fs
  dicrotic_amplitude <- ppg[dicrotic] - ppg[foot[1:n]]

  half_height <- (ppg[foot[1:n]] + ppg[peak]) / 2
  pw_half <- vapply(seq_len(n), function(i) {
    if (is.na(foot[i]) || is.na(peak[i]) || foot[i] >= peak[i]) return(NA_real_)
    seg_idx <- foot[i]:peak[i]; seg <- ppg[seg_idx]
    if (length(seg) < 3 || all(is.na(seg))) return(NA_real_)
    t_seg <- seg_idx / fs; hh <- half_height[i]
    above <- which(!is.na(seg) & seg >= hh)
    if (length(above) >= 2) t_seg[max(above)] - t_seg[min(above)] else {
      out <- tryCatch(approx(x = seg, y = t_seg, xout = hh)$y, error = function(e) NA_real_)
      if (length(out) == 2 && all(!is.na(out))) diff(range(out)) else NA_real_
    }
  }, numeric(1))

  max_dVdt <- vapply(seq_len(n), function(i) {
    seg <- ppg[foot[i]:peak[i]]
    if (length(seg) > 1) max(diff(seg), na.rm = TRUE) * fs else NA_real_
  }, numeric(1))

  data.frame(
    pulse_pressure, rise_time, AIx, RWTT, PAT, max_dVdt, HR, pulse_duration,
    ppg_amplitude = ppg_sys - ppg_dia,
    ppg_rise_time, s_d_ratio, dicrotic_timing, dicrotic_amplitude,
    ppg_width_half_height = pw_half
  )
}
