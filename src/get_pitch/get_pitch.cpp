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
    PitchAnalyzer analyzer(n_len, rate, PitchAnalyzer::RECT, 50, 500);
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
    vector<float> f0;
    for (vector<float>::iterator iX = x.begin(); iX + n_len < x.end(); iX += n_shift) {
        float f = analyzer(iX, iX + n_len);
        f0.push_back(f);
    }

    /// \TODO
    /// Postprocess the estimation in order to supress errors. For instance, a median filter
    /// or time-warping may be used.
    
    //if(1){
    if (apply_median) {
    //No aplico la media porque va peor entonces xd
    // Apply selective voice recovery filter - this should help with voiced-as-unvoiced errors
        cout << "Applying selective voice recovery...\n";
        f0 = recoverMissedVoicedFrames(f0);
        // Fix brief unvoiced segments in voiced regions
        cout << "Fixing brief unvoiced segments...\n";
        f0 = fixBriefUnvoicedSegments(f0, 2);
        
        // Add an additional targeted fix for isolated unvoiced frames
        // with very strict criteria to avoid false positives
        cout << "Applying precision voice recovery...\n";
        f0 = fixIsolatedFramesWithPitchConsistency(f0);
        //cout << "Applying enhanced median filtering...\n";
        //f0 = medianFilter(f0, 5);  
        
        // Fix isolated unvoiced frames with increased context and relaxed criteria
        //cout << "Fixing isolated unvoiced frames...\n";
        //f0 = fixIsolatedFrames(f0, 3);
        
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
