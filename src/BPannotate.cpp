#include <Rcpp.h>
#include <limits>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
List BP_resample(NumericVector waveform, double origFs) {
  double Fs = 200; // target resample frequency
  int n = waveform.size();

  double duration = n / origFs;
  int newLen = round(Fs * duration);

  NumericVector oldx(n);
  NumericVector newx(newLen);

  for(int i=0; i<n; i++){
    oldx[i] = i / origFs;
  }
  for(int i=0; i<newLen; i++){
    newx[i] = i / Fs;
  }

  NumericVector newWaveform(newLen);

  // linear interpolation
  for(int i=0; i<newLen; i++){
    double xi = newx[i];
    // find the left index
    int left = floor(xi * origFs);
    if(left < 0) left = 0;
    if(left >= n-1) left = n-2;
    double t0 = oldx[left];
    double t1 = oldx[left+1];
    double y0 = waveform[left];
    double y1 = waveform[left+1];
    double frac = (xi - t0)/(t1 - t0);
    newWaveform[i] = y0 + frac*(y1 - y0);
  }

  return List::create(
    Named("newWaveform") = newWaveform,
    Named("newx") = newx,
    Named("oldx") = oldx
  );
}

// [[Rcpp::export]]
NumericVector BP_lowpass(NumericVector x) {
  Rcpp::Environment signal = Rcpp::Environment::namespace_env("signal");
  Rcpp::Function butter  = signal["butter"];
  Rcpp::Function filtfilt = signal["filtfilt"];
  List bf = butter(4, 20.0/(200/2.0)); // 4th order Butterworth
  NumericVector filtered = filtfilt(bf, x);
  return filtered;
}



// [[Rcpp::export]]
  List doubleDerive(NumericVector waveform) {
    int n = waveform.size();

    NumericVector waveformD(n, NA_REAL);
    for (int i = 0; i < n - 1; i++) {
      waveformD[i] = waveform[i + 1] - waveform[i];
    }

    NumericVector waveformDD(n - 1, NA_REAL);
    for (int i = 0; i < n - 2; i++) {
      waveformDD[i] = waveformD[i + 1] - waveformD[i];
    }

    // apply lowpass to second derivative
    Rcpp::Environment signal = Rcpp::Environment::namespace_env("signal");
    Rcpp::Function butter  = signal["butter"];
    Rcpp::Function filtfilt = signal["filtfilt"];
    int fs = 200;
    List bf = butter(4, 20.0 / (fs / 2.0)); // 4th order, 20 Hz cutoff

    // Remove NAs manually
    std::vector<double> tmp;
    tmp.reserve(waveformDD.size());
    for (int i = 0; i < waveformDD.size(); ++i) {
      if (!NumericVector::is_na(waveformDD[i])) tmp.push_back(waveformDD[i]);
    }

    NumericVector waveformDD_filt_padded(n, NA_REAL); // pad to n (like your original)
    if (!tmp.empty()) {
      NumericVector validDD(tmp.begin(), tmp.end());
      NumericVector waveformDD_filt = filtfilt(bf, validDD);
      int m = waveformDD_filt.size();
      for (int i = 0; i < m && i < n; ++i) {
        waveformDD_filt_padded[i] = waveformDD_filt[i];
      }
    }

    // waveformDDPlus logic
    NumericVector waveformDDPlus(n, NA_REAL);
    for (int i = 0; i < n - 1; i++) {
      if ((waveformD[i] > 0) && (waveformDD_filt_padded[i] > 0)) {
        waveformDDPlus[i] = std::pow(waveformDD_filt_padded[i], 2);
      } else {
        waveformDDPlus[i] = 0.0;
      }
    }

    return List::create(
      Named("waveformD")     = waveformD,
      Named("waveformDD")    = waveformDD_filt_padded,
      Named("waveformDDPlus")= waveformDDPlus
    );
  }




// [[Rcpp::export]]
IntegerVector getFootIndex(NumericVector waveformDDPlus, LogicalVector zoneOfInterest) {
  int n = zoneOfInterest.size();

  IntegerVector zoneWall(n - 1, 0);
  for (int i = 0; i < n - 1; i++) {
    zoneWall[i] = zoneOfInterest[i + 1] - zoneOfInterest[i];
  }

  std::vector<int> BP_start, BP_stop;
  for (int i = 0; i < zoneWall.size(); i++) {
    if (zoneWall[i] == 1)  BP_start.push_back(i + 1); // 1-based
    if (zoneWall[i] == -1) BP_stop.push_back(i + 1);  // 1-based
  }

  // drop leading stop before first start
  while (!BP_stop.empty() && !BP_start.empty() && BP_stop.front() < BP_start.front()) {
    BP_stop.erase(BP_stop.begin());
  }

  int nfeet = std::min(BP_start.size(), BP_stop.size());
  IntegerVector footIndex(nfeet);

  for (int i = 0; i < nfeet; i++) {
    double maxval = -std::numeric_limits<double>::infinity();
    int maxidx = BP_start[i];
    for (int j = BP_start[i]; j < BP_stop[i]; j++) {
      // j is 1-based; convert to 0-based when indexing the C++ vector
      int jj = j - 1;
      if (jj >= 0 && jj < waveformDDPlus.size()) {
        if (waveformDDPlus[jj] > maxval) {
          maxval = waveformDDPlus[jj];
          maxidx = j; // keep 1-based for return
        }
      }
    }
    footIndex[i] = maxidx; // 1-based
  }
  return footIndex;
}


// [[Rcpp::export]]
IntegerVector FixIndex(IntegerVector footIndex, double fs = 200.0,
                       double minIntervalSec = 0.3,
                       double maxIntervalSec = 2.0) {
  if (footIndex.size() < 2) return footIndex;

  std::vector<int> fixed;
  fixed.push_back(footIndex[0]);

  for (int i = 1; i < footIndex.size(); ++i) {
    int prev = fixed.back();
    int curr = footIndex[i];
    double interval = (curr - prev) / fs;

    if (interval >= minIntervalSec && interval <= maxIntervalSec) {
      fixed.push_back(curr);
    }
    // else discard point as it's either too close or too far
  }

  return wrap(fixed);
}

// [[Rcpp::export]]
IntegerVector getPeakIndex(NumericVector waveform, IntegerVector footIndex) {
  int n = footIndex.size();
  IntegerVector peakIndex(n - 1);

  for (int i = 0; i < n - 1; ++i) {
    int start = footIndex[i]     - 1; // 0-based
    int end   = footIndex[i + 1] - 1; // 0-based

    if (end > start && end < waveform.size()) {
      double maxVal = waveform[start];
      int maxIdx = start;
      for (int j = start + 1; j < end; ++j) {
        if (waveform[j] > maxVal) { maxVal = waveform[j]; maxIdx = j; }
      }
      peakIndex[i] = maxIdx + 1;  // return 1-based
    } else {
      peakIndex[i] = NA_INTEGER;
    }
  }
  return peakIndex;
}


// [[Rcpp::export]]
DataFrame computeSUT_PA(NumericVector waveform, IntegerVector foot, IntegerVector peak, double fs) {
  int n = foot.size();
  NumericVector sut(n), pa(n);

  for (int i = 0; i < n; i++) {
    int f = foot[i] - 1;  // 0-based
    int p = peak[i] - 1;  // 0-based
    if (f >= 0 && p >= 0 && p < waveform.size() && f < waveform.size() && p > f) {
      sut[i] = (p - f) / fs;
      pa[i]  = waveform[p] - waveform[f];
    } else {
      sut[i] = NA_REAL; pa[i] = NA_REAL;
    }
  }

  return DataFrame::create(
    Named("foot") = foot,
    Named("peak") = peak,
    Named("SUT")  = sut,
    Named("PA")   = pa
  );
}


// [[Rcpp::export]]
DataFrame getDicroticNotchAndPeak(NumericVector signal, IntegerVector foot, IntegerVector peak, double fs = 200.0) {
  int n = std::min((int)foot.size() - 1, (int)peak.size());
  IntegerVector notch(n), dicrotic(n);

  for (int i = 0; i < n; ++i) {
    int i_start = peak[i]     - 1;         // 0-based
    int i_end   = foot[i + 1] - 1;         // 0-based

    if (i_end <= i_start || i_end >= signal.size()) {
      notch[i] = NA_INTEGER; dicrotic[i] = NA_INTEGER; continue;
    }

    double y1 = signal[i_start], y2 = signal[i_end];
    double slope = (y2 - y1) / (i_end - i_start);
    double intercept = y1 - slope * i_start;

    int len = i_end - i_start + 1;
    NumericVector residual(len);
    for (int j = 0; j < len; ++j) {
      int idx = i_start + j;
      double baseline = slope * idx + intercept;
      residual[j] = signal[idx] - baseline;
    }

    int rel_notch = which_min(residual);
    int notch_idx = i_start + rel_notch;   // 0-based
    notch[i] = notch_idx + 1;              // return 1-based

    int win_len = std::min((int)(0.25 * fs), (int)signal.size() - notch_idx - 1);
    if (win_len > 1) {
      double max_val = signal[notch_idx + 1];
      int max_idx = notch_idx + 1;
      for (int k = 1; k < win_len; ++k) {
        int idx = notch_idx + 1 + k;
        if (signal[idx] > max_val) { max_val = signal[idx]; max_idx = idx; }
      }
      dicrotic[i] = max_idx + 1;          // 1-based
    } else {
      dicrotic[i] = NA_INTEGER;
    }
  }

  return DataFrame::create(
    Named("notch")   = notch,
    Named("dicrotic")= dicrotic
  );
}
