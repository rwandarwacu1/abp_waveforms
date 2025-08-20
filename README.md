# abpwaveforms <a href="https://github.com/rwandarwacu1/abp_waveforms"><img src="https://img.shields.io/badge/github-abp__waveforms-blue?logo=github" alt="GitHub"></a>

**Arterial Blood Pressure (ABP), ECG, and PPG waveform feature extraction in R.**

`abpwaveforms` is an R package that provides reproducible, Rcpp-accelerated tools for cardiovascular waveform analysis.  
It is inspired by **BPannotate (MATLAB)** and enables researchers to extract key hemodynamic features from ABP, ECG, and PPG signals.

---

## ✨ Features

- 📉 **Fiducial detection**: foot (onset), systolic peak, dicrotic notch/peak  
- ⚡ **Enhanced feature extraction**:
  - Pulse pressure (PP)  
  - Augmentation index (AIx)  
  - Reflected wave transit time (RWTT)  
  - Pulse arrival time (PAT)  
  - Heart rate (HR)  
  - Max dV/dt  
  - PPG amplitude, rise time, half-width, dicrotic timing/amplitude    
- 📊 **Visualization utilities**: waveform plots, feature trends, summary tables  
- 🚀 **Fast Rcpp backend** for resampling, derivatives, and peak/notch detection

---

## 📦 Installation

```r
# install directly from GitHub
install.packages("devtools")
devtools::install_github("rwandarwacu1/abp_waveforms")
