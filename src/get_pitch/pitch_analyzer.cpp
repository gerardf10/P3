/// @file

#include <iostream>
#include <math.h>
#include "pitch_analyzer.h"

using namespace std;

/// Name space of UPC
namespace upc {
  void PitchAnalyzer::autocorrelation(const vector<float> &x, vector<float> &r) const {
    unsigned int N = x.size();
    // Para cada retardo l (lag)
    for (unsigned int l = 0; l < r.size(); ++l) {
        float sum = 0.0F;
        // r_{xx}[m] = \frac{1}{N} \sum_0^{N-m} x[n] x[n+m]
        // Se recorre la señal hasta el índice N - l para evitar acceder fuera de rango
        for (unsigned int n = 0; n < N - l; ++n) {
            sum += x[n] * x[n + l];
        }
        // Autocorrelación sesgada: se divide por el número total de muestras N
        r[l] = sum / static_cast<float>(N);
    }

    // Se ajusta r[0] para evitar problemas posteriores (como división por cero o logaritmos)
    if (r[0] == 0.0F)
        r[0] = 1e-10F;
  }


  void PitchAnalyzer::set_window(Window win_type) {
    if (frameLen == 0)
      return;

    window.resize(frameLen);

    switch (win_type) {
    case HAMMING:
      for (unsigned int n = 0; n < frameLen; ++n) {
        window[n] = 0.54f - 0.46f * cos(2.0f * 3.14159 * n / (frameLen - 1));
      }
      break;
    case RECT:
    default:
      window.assign(frameLen, 1);
    }
  }

  void PitchAnalyzer::set_f0_range(float min_F0, float max_F0) {
    npitch_min = (unsigned int) samplingFreq/max_F0;
    if (npitch_min < 2)
      npitch_min = 2;  // samplingFreq/2

    npitch_max = 1 + (unsigned int) samplingFreq/min_F0;

    //frameLen should include at least 2*T0
    if (npitch_max > frameLen/2)
      npitch_max = frameLen/2;
  }

  bool PitchAnalyzer::unvoiced(float pot, float norm_r1, float norm_rpeak) const {
    // Further refined thresholds based on error analysis
    const float threshold_power = -52.0f;    // Lowered to reduce voiced-as-unvoiced errors
    const float threshold_r1 = 0.45f;        // Reduced to better capture voiced frames
    const float threshold_rpeak = 0.42f;     // Reduced to catch more voiced frames
    
    // Special case for strong periodicity - clearly voiced
    if (norm_rpeak > 0.6f && pot > -48.0f) {
      return false;  // Very strong peak and decent power - definitely voiced
    }
    
    // Additional special case for borderline frames with good periodicity
    if (norm_rpeak > 0.52f && norm_r1 > 0.5f && pot > -55.0f) {
      return false;  // Good periodicity indicators even with lower power
    }
    
    // Keep the original OR structure but with refined thresholds
    if ((pot < threshold_power) ||
        (norm_r1 < threshold_r1) ||
        (norm_rpeak < threshold_rpeak))
    {
      // Expanded exception for borderline cases with good power and decent periodicity
      if ((pot > -42.0f && norm_rpeak > 0.36f && norm_r1 > 0.42f) ||
          (pot > -46.0f && norm_rpeak > 0.45f)) {
        return false;  // More cases classified as voiced
      }
      
      return true;  // Most frames will be classified as unvoiced
    } else {
      return false;  // Voiced frame
    }
  }

  
float PitchAnalyzer::compute_pitch(vector<float> & x) const {
  if (x.size() != frameLen)
      return -1.0F;

  // Apply window to input frame
  for (unsigned int i = 0; i < x.size(); ++i)
      x[i] *= window[i];

  vector<float> r(npitch_max);
  autocorrelation(x, r);

  // Improved frequency range definitions
  const unsigned int maxPitch = 400;
  const unsigned int minPitch = 80;
  unsigned int lag_min = samplingFreq / maxPitch;
  unsigned int lag_max = samplingFreq / minPitch;

  if (lag_max >= r.size())
      lag_max = r.size() - 1;

  // Enhanced peak finding logic with parabolic interpolation
  unsigned int best_lag = lag_min;
  float best_corr = r[lag_min];
  
  // First pass: find the coarse maximum peak
  for (unsigned int l = lag_min + 1; l <= lag_max; ++l) {
      if (r[l] > best_corr) {
          best_corr = r[l];
          best_lag = l;
      }
  }
  
  // Check if we have a reliable peak by comparing to neighboring values

  float pot = 10 * log10(r[0]);


  // Enhanced reliability assessment
  bool reliable_peak = (best_lag > lag_min && best_lag < lag_max) && 
  (r[best_lag] > r[best_lag-1] && r[best_lag] > r[best_lag+1]);

  // Additional check for strong harmonics (helps with voiced detection)
  if (!reliable_peak && best_lag > lag_min*2) {
  // Check if we have a harmonic peak at half the frequency (double the lag)
  unsigned int harmonic_lag = best_lag / 2;
    if (harmonic_lag >= lag_min && 
      r[harmonic_lag] > 0.38f * r[0] &&  
      r[harmonic_lag] > r[harmonic_lag-1] && 
      r[harmonic_lag] > r[harmonic_lag+1]) {
      reliable_peak = true;
      best_lag = harmonic_lag;
    }
  }
// Add this code before returning the pitch in compute_pitch():

// Apply parabolic interpolation for sub-sample pitch accuracy
float delta = 0.0f;
if (reliable_peak && best_lag > 0 && best_lag < r.size()-1) {
    float y1 = r[best_lag-1];
    float y2 = r[best_lag];
    float y3 = r[best_lag+1];
    delta = 0.5f * (y1 - y3) / (y1 - 2.0f*y2 + y3 + 1e-10f);
    // Limit delta to reasonable range
    if (delta < -0.5f) delta = -0.5f;
    if (delta > 0.5f) delta = 0.5f;
}

// Use the corrected pitch value
if (!reliable_peak || unvoiced(pot, r[1] / r[0], r[best_lag] / r[0]))
    return 0;
else
    return static_cast<float>(samplingFreq) / static_cast<float>(best_lag + delta);
}
}