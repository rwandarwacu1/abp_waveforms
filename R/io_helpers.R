#' Process a CSV with ABP/ECG/PPG columns
#' @param file CSV path containing columns ABP, ECG, PPG
#' @param fs Target Hz (default 200)
#' @return Data frame of features
#' @export
process_condition <- function(file, fs = 200) {
  dat <- utils::read.csv(file)
  extractFeaturesEnhanced(abp = dat$ABP, ppg = dat$PPG, ecg = dat$ECG, fs = fs)
}
