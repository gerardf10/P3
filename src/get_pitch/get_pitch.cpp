/// @file

#include <iostream>
#include <fstream>
#include <string.h>
#include <errno.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <map>
#include "wavfile_mono.h"
#include "pitch_analyzer.h"
#include "docopt.h"

#define FRAME_LEN   0.030 /* 30 ms. */
#define FRAME_SHIFT 0.015 /* 15 ms. */

using namespace std;
using namespace upc;

// Se añaden nuevas opciones: --clip, --lowpass y --median
static const char USAGE[] = R"(
get_pitch - Pitch Estimator

Usage:
    get_pitch [options] <input-wav> <output-txt>
    get_pitch (-h | --help)
    get_pitch --version

Options:
    -h, --help          Show this screen
    --version           Show the version of the project
    --clip              Enable central clipping preprocessing.
    --lowpass=<fc>      Apply low-pass filtering with cutoff frequency <fc> (Hz) [default: 8000]
    --median            Apply median filtering postprocessing.

Arguments:
    input-wav   Wave file with the audio signal
    output-txt  Output file: ASCII file with the result of the estimation:
                    - One line per frame with the estimated f0
                    - If considered unvoiced, f0 must be set to f0 = 0
)";

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

// Enhanced median filter with voiced continuity protection
// Enhanced median filter optimized for pitch continuity
vector<float> medianFilter(const vector<float>& f0, int window_size = 5) {
    vector<float> filtered = f0; // Start with original values
    int half = window_size / 2;
    
    for (size_t i = 0; i < f0.size(); ++i) {
        // Skip processing if this is a reliable voiced frame
        if (f0[i] > 0) {
            // Don't filter reliable frames with good context
            bool has_consistent_context = true;
            int context_voiced = 0;
            
            // Check consistency with surrounding frames
            for (int j = -2; j <= 2; ++j) {
                if (j == 0) continue;
                size_t idx = i + j;
                if (idx < f0.size() && f0[idx] > 0) {
                    context_voiced++;
                    float ratio = max(f0[i], f0[idx]) / min(f0[i], f0[idx]);
                    // Detect octave jumps (ratio close to 2) or large deviations
                    if (ratio > 1.8f && ratio < 2.2f) {
                        has_consistent_context = false; // Possible octave error
                    }
                    else if (ratio > 1.3f) {
                        has_consistent_context = false; // Large discontinuity
                    }
                }
            }
            
            // Skip filtering if frame is in a consistent voiced region
            if (has_consistent_context && context_voiced >= 2) {
                continue; // Keep original value
            }
        }
        
        // Collect values from window
        vector<float> window;
        for (int j = -half; j <= half; ++j) {
            int idx = i + j;
            if (idx >= 0 && idx < static_cast<int>(f0.size()))
                window.push_back(f0[idx]);
        }
        
        // Separate voiced and unvoiced values
        vector<float> voiced;
        int unvoiced_count = 0;
        
        for (float val : window) {
            if (val > 0) {
                voiced.push_back(val);
            } else {
                unvoiced_count++;
            }
        }
        
        // Apply different strategies based on frame context
        if (f0[i] == 0) {
            // Current frame is unvoiced
            if (voiced.size() > window.size() * 0.7f) {
                // Strong evidence of voicing in context
                
                // Find clusters in voiced values to avoid octave errors
                sort(voiced.begin(), voiced.end());
                
                // Find most common pitch range
                float best_pitch = 0;
                int max_cluster = 0;
                
                for (size_t v = 0; v < voiced.size(); v++) {
                    int cluster_size = 1;
                    for (size_t v2 = 0; v2 < voiced.size(); v2++) {
                        if (v == v2) continue;
                        float ratio = max(voiced[v], voiced[v2]) / min(voiced[v], voiced[v2]);
                        if (ratio < 1.2f) {
                            cluster_size++;
                        }
                    }
                    
                    if (cluster_size > max_cluster) {
                        max_cluster = cluster_size;
                        best_pitch = voiced[v];
                    }
                }
                
                if (max_cluster >= 3) {
                    // Found a strong cluster - correct to this value
                    filtered[i] = best_pitch;
                }
            }
        } else {
            // Current frame is voiced
            if (unvoiced_count > window.size() * 0.7f) {
                // Strong evidence of unvoicing in context
                filtered[i] = 0;
            } else if (!voiced.empty()) {
                // Look for possible octave errors
                sort(voiced.begin(), voiced.end());
                
                // Compare current pitch with median in context
                float median_pitch = voiced[voiced.size()/2];
                float ratio = max(f0[i], median_pitch) / min(f0[i], median_pitch);
                
                // Detect and fix potential octave errors
                if (ratio > 1.8f && ratio < 2.2f) {
                    // Potential octave error - correct towards neighborhood
                    filtered[i] = median_pitch;
                }
                else if (ratio > 1.5f) {
                    // Large jump - use weighted correction
                    filtered[i] = 0.7f * f0[i] + 0.3f * median_pitch;
                }
            }
        }
    }
    
    return filtered;
}





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

int main(int argc, const char *argv[]) {
    /// \TODO 
    ///  Modify the program syntax and the call to **docopt()** in order to
    ///  add options and arguments to the program.
    std::map<std::string, docopt::value> args = docopt::docopt(USAGE,
                                                                {argv + 1, argv + argc},
                                                                true,
                                                                "2.0");

    std::string input_wav = args["<input-wav>"].asString();
    std::string output_txt = args["<output-txt>"].asString();

    // Opciones de pre y postprocesado extraídas de los argumentos
    bool apply_clip = args["--clip"].asBool();
    bool apply_median = args["--median"].asBool();
    float lowpassCutoff = std::stof(args["--lowpass"].asString());

    // Read input sound file
    unsigned int rate;
    vector<float> x;
    if (readwav_mono(input_wav, rate, x) != 0) {
        cerr << "Error reading input file " << input_wav << " (" << strerror(errno) << ")\n";
        return -2;
    }

    int n_len = rate * FRAME_LEN;
    int n_shift = rate * FRAME_SHIFT;

    // Define analyzer.
    // Aquí se definen por ejemplo los límites de pitch deseados. 
    // En este ejemplo, se utiliza: min pitch = 50 Hz y max pitch = 500 Hz.
        // Pitch estimation with a Tukey window
        upc::PitchAnalyzer pitchAnalyzer( n_len,
                                        rate,
                                        upc::PitchAnalyzer::HAMMING,
                                        /*min_F0=*/80,
                                        /*max_F0=*/350 );
   // force Tukey estimator to *always* output some F0
      pitchAnalyzer.setSkipUnvoicedTest(true);
        // VAD decision with a Rectangular window
        upc::PitchAnalyzer vadAnalyzer(   n_len,
                                        rate,
                                        upc::PitchAnalyzer::RECT,
                                        /*min_F0=*/80,
                                        /*max_F0=*/350 );

    
    /// \TODO
    /// Preprocess the input signal in order to ease pitch estimation. For instance,
    /// central-clipping or low pass filtering may be used.
    //if(1){
    if(apply_clip) {
        cout << "Applying central clipping preprocessing...\n";
        centralClipping(x);
    }
    // Aquí podrías añadir también un filtrado paso bajo, si lo deseas.
    // Dado que se pasa un cutoff (lowpassCutoff) se podría aplicar un filtro (no implementado aquí)
    // Por ejemplo: lowPassFilter(x, rate, lowpassCutoff);
    if (lowpassCutoff > 0 && lowpassCutoff < rate/2) {
        // Use a lower cutoff to focus on fundamental frequency range
        float effective_cutoff = min(lowpassCutoff, 900.0f); // Cap at 900 Hz for better voice focus
        cout << "Applying low-pass filter with cutoff " << effective_cutoff << " Hz...\n";
        lowPassFilter(x, rate, effective_cutoff);
    }
    // Iterate for each frame and save values in f0 vector
    vector<float> frame(n_len);
    vector<float> f0;
    for (auto it = x.begin(); it + n_len < x.end(); it += n_shift) {
        std::copy(it, it + n_len, frame.begin());

        float pot2, nr1_2, nrp2; // For vadAnalyzer
        unsigned int lag2; bool rp2;
        // Declare dummy variables for correlation values for vadAnalyzer call
        float r_m1_vad, r_0_vad, r_p1_vad; 
        vadAnalyzer.analyze_frame(frame, pot2, nr1_2, nrp2, lag2, rp2, r_m1_vad, r_0_vad, r_p1_vad); // Pass new arguments

        if (vadAnalyzer.unvoiced(pot2, nr1_2, nrp2)) {
            f0.push_back(0.0f);
        } else {
            float pot1, nr1_1, nrp1;
            unsigned int lag1; bool rp1;
            // These are already correctly declared for pitchAnalyzer
            float r_m1, r_0, r_p1; 
            pitchAnalyzer.analyze_frame(frame, pot1, nr1_1, nrp1, lag1, rp1, r_m1, r_0, r_p1);
            f0.push_back(pitchAnalyzer.lag_to_f0(lag1, rp1, pot1, nr1_1, nrp1, r_m1, r_0, r_p1));
        }
    }


    /// \TODO
    /// Postprocess the estimation in order to supress errors. For instance, a median filter
    /// or time-warping may be used.
    
    //if(1){
if (apply_median) {
    
    cout << "Applying enhanced median filtering...\n";
    f0 = medianFilter(f0, 3);

    cout << "Fixing brief unvoiced segments...\n";
    f0 = fixBriefUnvoicedSegments(f0, 2);
    

    cout << "Applying precision voice recovery...\n";
    f0 = fixIsolatedFramesWithPitchConsistency(f0);
    
    cout << "Applying selective voice recovery...\n";
    f0 = recoverMissedVoicedFrames(f0);
    

    

}
    

    // Write f0 contour into the output file
    ofstream os(output_txt);
    if (!os.good()) {
        cerr << "Error reading output file " << output_txt << " (" << strerror(errno) << ")\n";
        return -3;
    }

    // Escribir la estimación: se incluye pitch=0 al inicio y al final de la señal
    os << 0 << '\n'; // pitch at t=0
    for (vector<float>::iterator i = f0.begin(); i != f0.end(); ++i) 
        os << *i << '\n';
    os << 0 << '\n'; // pitch at t=Duration

    return 0;
}