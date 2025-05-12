PAV - P3: estimación de pitch
=============================

Guillem Moreno Garcia i Gerard Ferret Colomer
---------------------------------------------

Esta práctica se distribuye a través del repositorio GitHub [Práctica 3](https://github.com/albino-pav/P3).
Siga las instrucciones de la [Práctica 2](https://github.com/albino-pav/P2) para realizar un `fork` de la
misma y distribuir copias locales (*clones*) del mismo a los distintos integrantes del grupo de prácticas.

Recuerde realizar el *pull request* al repositorio original una vez completada la práctica.

Ejercicios básicos
------------------

- Complete el código de los ficheros necesarios para realizar la estimación de pitch usando el programa
  `get_pitch`.

   * Complete el cálculo de la autocorrelación e inserte a continuación el código correspondiente.
	```cpp
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
 	```

   * Inserte una gŕafica donde, en un *subplot*, se vea con claridad la señal temporal de un segmento de
     unos 30 ms de un fonema sonoro y su periodo de pitch; y, en otro *subplot*, se vea con claridad la
	 autocorrelación de la señal y la posición del primer máximo secundario.

	 NOTA: es más que probable que tenga que usar Python, Octave/MATLAB u otro programa semejante para
	 hacerlo. Se valorará la utilización de la biblioteca matplotlib de Python.

		Subplot generat amb el codi "subplot.py" ubicat a la carpeta "scripts".

		![Subplot](/subplot.png)

   * Determine el mejor candidato para el periodo de pitch localizando el primer máximo secundario de la
     autocorrelación. Inserte a continuación el código correspondiente.

```cpp
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
```

   * Implemente la regla de decisión sonoro o sordo e inserte el código correspondiente.

```cpp
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
```

   * Puede serle útil seguir las instrucciones contenidas en el documento adjunto `código.pdf`.

- Una vez completados los puntos anteriores, dispondrá de una primera versión del estimador de pitch. El 
  resto del trabajo consiste, básicamente, en obtener las mejores prestaciones posibles con él.

  * Utilice el programa `wavesurfer` para analizar las condiciones apropiadas para determinar si un
    segmento es sonoro o sordo. 
	
	  - Inserte una gráfica con la estimación de pitch incorporada a `wavesurfer` y, junto a ella, los 
	    principales candidatos para determinar la sonoridad de la voz: el nivel de potencia de la señal
		(r[0]), la autocorrelación normalizada de uno (r1norm = r[1] / r[0]) y el valor de la
		autocorrelación en su máximo secundario (rmaxnorm = r[lag] / r[0]).

		Puede considerar, también, la conveniencia de usar la tasa de cruces por cero.

	    Recuerde configurar los paneles de datos para que el desplazamiento de ventana sea el adecuado, que
		en esta práctica es de 15 ms.

		Subplot generat amb el codi "candidats.py" ubicat a la carpeta "scripts".

		![Candidats](/candidats.png)

      - Use el estimador de pitch implementado en el programa `wavesurfer` en una señal de prueba y compare
	    su resultado con el obtenido por la mejor versión de su propio sistema.  Inserte una gráfica
		ilustrativa del resultado de ambos estimadores.
     
		Aunque puede usar el propio Wavesurfer para obtener la representación, se valorará
	 	el uso de alternativas de mayor calidad (particularmente Python).

		Subplot generat amb el codi "comparacio.py" ubicat a la carpeta "scripts".

		![Comparacio](/comparacio.png)
  
  * Optimice los parámetros de su sistema de estimación de pitch e inserte una tabla con las tasas de error
    y el *score* TOTAL proporcionados por `pitch_evaluate` en la evaluación de la base de datos 
	`pitch_db/train`..

	Fent servir el fitxer "run_get_pitch.sh" que trobem a la carpeta "scripts", principalment, podem constatar que d'acord amb la decisió que pren "unvoiced" el més important són els llindars de l'autocorrelació normalitzada d'1 (r1norm) i del valor de l'autocorrelació al seu màxim secundari (rmaxnorm). D'aquesta manera aconseguim una puntuació màxima total de 92.3%.

Taula:
| Error type                  | Number of errors      | %      |
|----------------------------|-----------------------|--------|
| Unvoiced frames as voiced  | 235/7045              | 3.05   |
| Voiced frames as unvoiced  | 268/4155              | 7.10   |
| Gross voiced errors (+20%) | 41/3887               | 0.67   |
| MSE of fine errors         |                       | 2.71   |
| **TOTAL**                  |                       | **92.39** |

Captura de pantalla:

![Score](/score.png)


Ejercicios de ampliación
------------------------

- Usando la librería `docopt_cpp`, modifique el fichero `get_pitch.cpp` para incorporar los parámetros del
  estimador a los argumentos de la línea de comandos.
  
  Esta técnica le resultará especialmente útil para optimizar los parámetros del estimador. Recuerde que
  una parte importante de la evaluación recaerá en el resultado obtenido en la estimación de pitch en la
  base de datos.

  * Inserte un *pantallazo* en el que se vea el mensaje de ayuda del programa y un ejemplo de utilización
    con los argumentos añadidos.

    ![get_pitch1](/get_pitch1.png)

    ![get_pitch2](/get_pitch2.png)

- Implemente las técnicas que considere oportunas para optimizar las prestaciones del sistema de estimación
  de pitch.

  Entre las posibles mejoras, puede escoger una o más de las siguientes:

  * Técnicas de preprocesado: filtrado paso bajo, diezmado, *center clipping*, etc.
  * Técnicas de postprocesado: filtro de mediana, *dynamic time warping*, etc.
  * Métodos alternativos a la autocorrelación: procesado cepstral, *average magnitude difference function*
    (AMDF), etc.
  * Optimización **demostrable** de los parámetros que gobiernan el estimador, en concreto, de los que
    gobiernan la decisión sonoro/sordo.
  * Cualquier otra técnica que se le pueda ocurrir o encuentre en la literatura.

  Encontrará más información acerca de estas técnicas en las [Transparencias del Curso](https://atenea.upc.edu/pluginfile.php/2908770/mod_resource/content/3/2b_PS%20Techniques.pdf)
  y en [Spoken Language Processing](https://discovery.upc.edu/iii/encore/record/C__Rb1233593?lang=cat).
  También encontrará más información en los anexos del enunciado de esta práctica.

  Incluya, a continuación, una explicación de las técnicas incorporadas al estimador. Se valorará la
  inclusión de gráficas, tablas, código o cualquier otra cosa que ayude a comprender el trabajo realizado.

  También se valorará la realización de un estudio de los parámetros involucrados. Por ejemplo, si se opta
  por implementar el filtro de mediana, se valorará el análisis de los resultados obtenidos en función de
  la longitud del filtro.

  Després del que es va comentar a la sessió de laboratori i després que algunes d'aquestes tècniques s'hagin vist a classe, així com haver consultat la documentació que se'ns posava a disposició i haver navegat per Internet, hem implementat el center clipping, el filtre pasbaix i el filtre de mediana, entre d'altres millores que es detallen a continuació.

- Center clipping:

```cpp
// Función auxiliar para aplicar central clipping a la señal
// Enhanced central clipping with adaptive thresholding
// Improved central clipping with gentler thresholding
void centralClipping(vector<float>& x, float threshold_factor = 0.3f) {
    // Process in frames for better local adaptation
    const int frame_size = 320; // ~20ms at 16kHz
    
    for (size_t frame_start = 0; frame_start < x.size(); frame_start += frame_size) {
        size_t frame_end = min(frame_start + frame_size, x.size());
        
        // Find peak amplitude in this frame
        float peak = 0.0f;
        for (size_t i = frame_start; i < frame_end; ++i) {
            if (fabs(x[i]) > peak) peak = fabs(x[i]);
        }
        
        // Set adaptive threshold as a percentage of the frame peak
        float threshold = peak * threshold_factor;
        
        // Apply central clipping with preserved amplitudes
        for (size_t i = frame_start; i < frame_end; ++i) {
            if (fabs(x[i]) < threshold) {
                x[i] = 0.0f;
            } else if (x[i] > 0) {
                x[i] = x[i] - threshold; // Preserve amplitude information
            } else {
                x[i] = x[i] + threshold; // Preserve amplitude information
            }
        }
    }
}
```


- Filtre de mediana:

```cpp
// Enhanced median filter with voiced continuity protection
vector<float> medianFilter(const vector<float>& f0, int window_size = 5) {
    vector<float> filtered(f0.size());
    int half = window_size / 2;
    
    // First pass: standard median filtering
    for (size_t i = 0; i < f0.size(); ++i) {
        vector<float> window;
        for (int j = -half; j <= half; ++j) {
            int idx = i + j;
            if (idx < 0) 
                window.push_back(f0.front());
            else if (idx >= static_cast<int>(f0.size()))
                window.push_back(f0.back());
            else
                window.push_back(f0[idx]);
        }
        // Filter out zeros before finding median (for better pitch continuity)
        vector<float> non_zero;
        for (float val : window) {
            if (val > 0) non_zero.push_back(val);
        }
        
        if (non_zero.empty()) {
            filtered[i] = 0.0f; // No voiced frames in window
        } else {
            // Find median of non-zero values
            sort(non_zero.begin(), non_zero.end());
            filtered[i] = non_zero[non_zero.size()/2];
            
            // If original was unvoiced but surrounded by similar voiced frames, use median
            if (f0[i] == 0.0f) {
                // Count voiced frames in window
                int voiced_count = 0;
                for (float val : window) {
                    if (val > 0) voiced_count++;
                }
                
                // If less than half are voiced, keep as unvoiced
                if (voiced_count <= window_size/2) {
                    filtered[i] = 0.0f;
                }
            }
        }
    }
    
    return filtered;
}
```


- Fix isolated unvoiced segments:
  
```cpp
// New function: Fix isolated unvoiced frames in voiced regions
// Enhanced isolated frame correction with relaxed criteria
vector<float> fixIsolatedFrames(const vector<float>& f0, int context_size = 3) {
    vector<float> fixed = f0;
    
    for (size_t i = context_size; i < f0.size() - context_size; ++i) {
        // Check for unvoiced frame that might be incorrectly classified
        if (f0[i] == 0.0f) {
            // Count surrounding voiced frames
            int voiced_count = 0;
            float sum_pitch = 0.0f;
            
            for (int j = -context_size; j <= context_size; ++j) {
                if (j == 0) continue; // Skip current frame
                
                if (f0[i+j] > 0.0f) {
                    voiced_count++;
                    sum_pitch += f0[i+j];
                }
            }
            
            // More relaxed condition: only 2/3 of surrounding frames need to be voiced
            float threshold = 2.0f * context_size / 3.0f;
            if (voiced_count >= threshold) {
                fixed[i] = sum_pitch / voiced_count;
            }
        }
    }
    
    return fixed;
}
```


- Fix brief unvoiced segments:
  
```cpp
// New function to fix brief unvoiced segments in voiced regions
vector<float> fixBriefUnvoicedSegments(const vector<float>& f0, int max_length = 3) {
    vector<float> fixed = f0;
    
    for (size_t i = max_length; i < f0.size() - max_length; ++i) {
        // Check if we're at the start of an unvoiced segment
        if (f0[i] == 0.0f && f0[i-1] > 0.0f) {
            // Find length of unvoiced segment
            int unvoiced_length = 0;
            while (i + unvoiced_length < f0.size() && f0[i+unvoiced_length] == 0.0f) {
                unvoiced_length++;
            }
            
            // Only fix short segments that are followed by voiced frames
            if (unvoiced_length <= max_length && i+unvoiced_length < f0.size() && f0[i+unvoiced_length] > 0.0f) {
                // Linear interpolation between surrounding voiced frames
                float start_pitch = f0[i-1];
                float end_pitch = f0[i+unvoiced_length];
                
                for (int j = 0; j < unvoiced_length; j++) {
                    float t = static_cast<float>(j+1) / (unvoiced_length+1);
                    fixed[i+j] = start_pitch * (1-t) + end_pitch * t;
                }
            }
            
            // Skip the unvoiced segment we just processed
            i += unvoiced_length;
        }
    }
    
    return fixed;
}
```


- Filtre pasbaix:
  
 ```cpp
// Low-pass filter implementation
void lowPassFilter(vector<float>& x, int sampleRate, float cutoffFreq) {
    // Simple first-order IIR low-pass filter
    // y[n] = alpha * x[n] + (1-alpha) * y[n-1]
    
    // Calculate alpha based on cutoff frequency
    float dt = 1.0f / static_cast<float>(sampleRate);
    float RC = 1.0f / (2.0f * 3.14159 * cutoffFreq);
    float alpha = dt / (dt + RC);
    
    // Apply filter
    float y_prev = x[0];
    for (size_t i = 0; i < x.size(); i++) {
        float y = alpha * x[i] + (1.0f - alpha) * y_prev;
        x[i] = y;
        y_prev = y;
    }
    
    cout << "Applied low-pass filter with cutoff: " << cutoffFreq << " Hz\n";
}
```


- Selective voice recovery filter:
  
 ```cpp
// Selective voice recovery filter - focus specifically on likely voiced frames
vector<float> recoverMissedVoicedFrames(const vector<float>& f0) {
    vector<float> fixed = f0;
    
    // First pass: identify regions of stable pitch
    vector<float> pitch_stability(f0.size(), 0.0f);
    for (size_t i = 3; i < f0.size() - 3; ++i) {
        if (f0[i] > 0) {
            float prev_pitch = 0.0f;
            int count = 0;
            float sum_deviation = 0.0f;
            
            // Compute average deviation from surrounding voiced frames
            for (int j = -3; j <= 3; ++j) {
                if (j == 0) continue;
                
                if (i+j >= 0 && i+j < f0.size() && f0[i+j] > 0) {
                    if (prev_pitch > 0) {
                        sum_deviation += fabs(f0[i+j] - prev_pitch) / prev_pitch;
                    }
                    prev_pitch = f0[i+j];
                    count++;
                }
            }
            
            if (count > 1) {
                pitch_stability[i] = 1.0f - (sum_deviation / count);  // Higher value = more stable
            }
        }
    }
    
    // Second pass: recover missed voiced frames in stable regions
    for (size_t i = 3; i < f0.size() - 3; ++i) {
        if (f0[i] == 0.0f) {  // Unvoiced frame
            // Check if surrounded by stable voiced frames
            int voiced_neighbors = 0;
            float avg_stability = 0.0f;
            float sum_pitch = 0.0f;
            
            for (int j = -2; j <= 2; ++j) {
                if (j == 0) continue;
                
                if (i+j >= 0 && i+j < f0.size() && f0[i+j] > 0) {
                    voiced_neighbors++;
                    avg_stability += pitch_stability[i+j];
                    sum_pitch += f0[i+j];
                }
            }
            
            // Only convert frames that are surrounded by stable voiced frames
            if (voiced_neighbors >= 3 && avg_stability / voiced_neighbors > 0.85f) {
                fixed[i] = sum_pitch / voiced_neighbors;
            }
        }
    }
    
    return fixed;
}
```

- Highly selective isolated frame correction:
  
 ```cpp
// Highly selective isolated frame correction
vector<float> fixIsolatedFramesWithPitchConsistency(const vector<float>& f0) {
    vector<float> fixed = f0;
    
    for (size_t i = 2; i < f0.size() - 2; ++i) {
        // Only target isolated unvoiced frames surrounded by voiced frames
        if (f0[i] == 0.0f && 
            f0[i-1] > 0.0f && f0[i-2] > 0.0f && 
            f0[i+1] > 0.0f && f0[i+2] > 0.0f) {
            
            // Check pitch consistency of surrounding frames
            float pitch_diff = fabs(f0[i-1] - f0[i+1]) / max(f0[i-1], f0[i+1]);
            
            // Only interpolate if surrounding pitches are very consistent (within 10%)
            if (pitch_diff < 0.10) {
                // Use average of neighboring frames
                fixed[i] = (f0[i-1] + f0[i+1]) / 2.0f;
            }
        }
    }
    
    return fixed;
}
```

Evaluación *ciega* del estimador
-------------------------------

Antes de realizar el *pull request* debe asegurarse de que su repositorio contiene los ficheros necesarios
para compilar los programas correctamente ejecutando `make release`.

Con los ejecutables construidos de esta manera, los profesores de la asignatura procederán a evaluar el
estimador con la parte de test de la base de datos (desconocida para los alumnos). Una parte importante de
la nota de la práctica recaerá en el resultado de esta evaluación.
