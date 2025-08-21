# abpwaveforms <a href="https://github.com/rwandarwacu1/abp_waveforms"><img src="https://img.shields.io/badge/github-abp__waveforms-blue?logo=github" alt="GitHub"></a>

# abpwaveforms <a href="https://github.com/rwandarwacu1/abpwaveforms"><img src="https://img.shields.io/badge/github-abpwaveforms-blue?logo=github" alt="GitHub"></a>

**Arterial Blood Pressure (ABP), ECG, and PPG waveform feature extraction in R.**

`abpwaveforms` is an R package that provides reproducible, Rcpp-accelerated tools for cardiovascular waveform analysis.  
Inspired by MATLAB **BPannotate**, it extracts key hemodynamic features from ABP, ECG, and PPG signals.

---

## ✨ Features

- 📉 **Fiducials**: foot (onset), systolic peak, dicrotic notch/peak  
- ⚡ **Features**: Pulse pressure (PP), Augmentation Index (AIx), Reflected Wave Transit Time (RWTT),
  Pulse Arrival Time (PAT), Heart Rate (HR), max dV/dt, PPG amplitude/rise time/half-width,
  dicrotic timing & amplitude  
- 🚀 **Fast Rcpp backend** for resampling, derivatives, and peak/notch detection

---

## 📦 Installation

```r
install.packages("devtools")
devtools::install_github("rwandarwacu1/abpwaveforms")
