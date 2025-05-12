/// @file

#include <iostream>
#include <math.h>
#include "pitch_analyzer.h"

using namespace std;
void makeTukeyWindow(vector<float>& w, float alpha) {
    int N = w.size();
    int edge = int(alpha * (N-1) / 2);
    for (int n = 0; n < N; ++n) {
        if (n < edge)
            w[n] = 0.5f * (1 + cos( 3.14159 * (2.0f*n/alpha/(N-1) - 1) ));
        else if (n <= N-1-edge)
            w[n] = 1.0f;
        else
            w[n] = 0.5f * (1 + cos( 3.14159 * (2.0f*(n-(N-1- edge))/alpha/(N-1) + 1) ));
    }
}

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
    case TUKEY:
      // Tukey window with alpha = 0.5
      makeTukeyWindow(window, 0.2f);
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

    // 1) Ventaneo
    for (unsigned int i = 0; i < x.size(); ++i)
        x[i] *= window[i];

    // 2) Autocorrelación
    vector<float> r(npitch_max);
    autocorrelation(x, r);

    // 3) Definir rango de búsqueda
    const unsigned int maxPitch = 400;
    const unsigned int minPitch = 80;
    unsigned int lag_min = samplingFreq / maxPitch;
    unsigned int lag_max = samplingFreq / minPitch;
    if (lag_max >= r.size()) lag_max = r.size() - 1;

    // 4) Búsqueda del pico grueso
    unsigned int best_lag = lag_min;
    float best_corr = r[lag_min];
    for (unsigned int l = lag_min + 1; l <= lag_max; ++l) {
        if (r[l] > best_corr) {
            best_corr = r[l];
            best_lag = l;
        }
    }

    // 5) Medida de potencia y normalizaciones
    float pot = 10 * log10(r[0]);
    float norm_r1   = r[1]        / r[0];
    float norm_rpeak = r[best_lag] / r[0];

    // 6) Validar fiabilidad del pico
    bool reliable_peak = (best_lag > lag_min && best_lag < lag_max)
        && (r[best_lag] > r[best_lag-1] && r[best_lag] > r[best_lag+1]);

    // <<< NUEVO >>> 7) Interpolación parabólica para sub-muestra
    float delta = 0.0f;
    if (reliable_peak && best_lag > 0 && best_lag < r.size()-1) {
        float y1 = r[best_lag - 1];
        float y2 = r[best_lag];
        float y3 = r[best_lag + 1];
        // fórmula Δ = 0.5*(y1 - y3)/(y1 - 2*y2 + y3)
        delta = 0.5f * (y1 - y3) / ( (y1 - 2.0f*y2 + y3) + 1e-10f );
        // acotar Δ para evitar saltos excesivos
        delta = std::max(std::min(delta, 0.5f), -0.5f);
    }

    // 8) Decidir voiced/unvoiced y devolver f0
    if (!reliable_peak || unvoiced(pot, norm_r1, norm_rpeak)) {
        return 0.0f;
    } else {
        // uso best_lag + delta para frecuencia sub-muestral
        return samplingFreq / (best_lag + delta);
    }
}





}
