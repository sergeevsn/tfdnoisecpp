#include "tfd.h"
#include <fftw3.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <numeric>

namespace tfd {

// Optimized window calculation 
real_vector window_hann(size_t n) {
    real_vector window(n);
    const float factor = 2.0f * M_PI / (n - 1);
    for (size_t i = 0; i < n; ++i) {
        window[i] = 0.5f * (1.0f - std::cos(factor * i));
    }
    return window;
}

// Vectorized reflection padding
real_vector reflect_padding(const real_vector& trace, size_t target_length) {
    real_vector padded = trace;
    while (padded.size() < target_length) {
        real_vector reflected(trace.rbegin(), trace.rend()-1);
        padded.insert(padded.end(), reflected.begin(), reflected.end());
    }
    return padded;
}

// Optimized median calculation using nth_element 
float compute_median(real_vector values) {
    if (values.empty()) return 0.0f;
    auto mid = values.begin() + values.size()/2;
    std::nth_element(values.begin(), mid, values.end());
    return values.size() % 2 ? *mid : (*mid + *std::max_element(values.begin(), mid)) / 2.0f;
}

complex_spectrogram compute_stft(const seismic_data& seismic, size_t n_fft, size_t n_overlap) {
    if (n_fft == 0 || n_overlap >= n_fft)
        throw std::invalid_argument("Invalid STFT parameters");

    complex_spectrogram spectrograms;
    const size_t window_step = n_fft - n_overlap;
    const real_vector window = window_hann(n_fft); // Precomputed once 

    // Reuse FFTW buffers
    std::vector<float> in(n_fft);
    std::vector<fftwf_complex> out(n_fft/2 + 1);
    fftwf_plan plan = fftwf_plan_dft_r2c_1d(n_fft, in.data(), out.data(), FFTW_MEASURE);

    for (const auto& trace : seismic) {
        if (trace.empty()) continue;

        size_t n_frames = (trace.size() - n_overlap + window_step) / window_step;
        size_t target_length = n_frames * window_step + n_fft;
        const real_vector padded = reflect_padding(trace, target_length);

        complex_vector spectrum(n_frames * (n_fft/2 + 1));
        for (size_t frame = 0; frame < n_frames; ++frame) {
            size_t offset = frame * window_step;

            // Apply window using vectorized operations
            std::transform(padded.begin() + offset, 
                          padded.begin() + offset + n_fft,
                          window.begin(),
                          in.begin(),
                          std::multiplies<>());

            fftwf_execute(plan);

            for (size_t i = 0; i < n_fft/2 + 1; ++i) {
                spectrum[frame*(n_fft/2 +1)+i] = {out[i][0], out[i][1]};
            }
        }
        spectrograms.emplace_back(std::move(spectrum));
    }

    fftwf_destroy_plan(plan);
    return spectrograms;
}

seismic_data compute_istft(const complex_spectrogram& spectrograms, 
                          size_t n_fft, size_t n_overlap, 
                          size_t original_length) 
{
    seismic_data reconstructed;
    const size_t window_step = n_fft - n_overlap;
    const real_vector window = window_hann(n_fft);
    const real_vector window_sq = [&]{
        real_vector sq(n_fft);
        std::transform(window.begin(), window.end(), sq.begin(),
                      [](float w){ return w*w; });
        return sq;
    }(); // Precompute squared window

    std::vector<fftwf_complex> in(n_fft/2 + 1);
    std::vector<float> out(n_fft);
    fftwf_plan plan = fftwf_plan_dft_c2r_1d(n_fft, in.data(), out.data(), FFTW_MEASURE);

    for (const auto& spectrum : spectrograms) {
        const size_t n_frames = spectrum.size() / (n_fft/2 + 1);
        const size_t output_size = (n_frames-1)*window_step + n_fft;
        
        real_vector output(output_size, 0.0f);
        real_vector overlap_count(output_size, 0.0f);

        for (size_t frame = 0; frame < n_frames; ++frame) {
            size_t offset = frame * window_step;

            // Copy frequency data
            for (size_t i = 0; i < n_fft/2 +1; ++i) {
                in[i][0] = spectrum[frame*(n_fft/2+1)+i].real();
                in[i][1] = spectrum[frame*(n_fft/2+1)+i].imag();
            }

            fftwf_execute(plan);

            // Overlap-add with precomputed window 
            for (size_t i = 0; i < n_fft; ++i) {
                if (offset + i >= output_size) break;
                output[offset+i] += out[i] * window[i] / n_fft;
                overlap_count[offset+i] += window_sq[i];
            }
        }

        // Normalize using precomputed values [[7]]
        std::transform(output.begin(), output.end(), overlap_count.begin(), output.begin(),
                      [](float val, float cnt) { return cnt > 1e-6f ? val/cnt : val; });

        // Resize to original length [[5]]
        if (output.size() > original_length) {
            output.resize(original_length);
        } else if (output.size() < original_length) {
            output.insert(output.end(), original_length - output.size(), output.back());
        }

        reconstructed.emplace_back(std::move(output));
    }

    fftwf_destroy_plan(plan);
    return reconstructed;
}

complex_spectrogram tfd_noise_rejection(
    const complex_spectrogram& stft_data,
    size_t trace_aperture,
    const std::vector<std::pair<float, float>>& threshold_params,
    float sample_rate,
    size_t n_fft,
    size_t n_overlap)
{
    complex_spectrogram filtered = stft_data;
    const size_t n_traces = stft_data.size();
    if (n_traces == 0) return filtered;

    const size_t n_freq = n_fft/2 + 1;
    const size_t window_step = n_fft - n_overlap;

    // Precompute global median [[5]]
    real_vector all_amplitudes;
    all_amplitudes.reserve(n_traces * stft_data[0].size());
    for (const auto& trace : stft_data) {
        for (const auto& c : trace) {
            all_amplitudes.push_back(std::abs(c));
        }
    }
    float global_median = compute_median(all_amplitudes);

    // Precompute phases [[7]]
    std::vector<real_vector> phases(n_traces);
    for (size_t i = 0; i < n_traces; ++i) {
        phases[i].resize(stft_data[i].size());
        std::transform(stft_data[i].begin(), stft_data[i].end(), phases[i].begin(),
                      [](const complex& c) { return std::arg(c); });
    }

    // Precompute multipliers table [[5]]
    auto compute_multiplier = [&](float time) {
        if (threshold_params.size() == 1) return threshold_params[0].second;
        
        auto it = std::lower_bound(threshold_params.begin(), threshold_params.end(), time,
            [](const auto& p, float t) { return p.first < t; });
        size_t idx = std::distance(threshold_params.begin(), it);
        
        if (idx == 0) return threshold_params[0].second;
        if (idx == threshold_params.size()) return threshold_params.back().second;
        
        float t0 = threshold_params[idx-1].first;
        float t1 = threshold_params[idx].first;
        float m0 = threshold_params[idx-1].second;
        float m1 = threshold_params[idx].second;
        return m0 + (time - t0)*(m1 - m0)/(t1 - t0);
    };

    for (size_t i = 0; i < n_traces; ++i) {
        size_t n_frames = stft_data[i].size() / n_freq;
        std::vector<float> frame_times(n_frames);
        for (size_t t = 0; t < n_frames; ++t) {
            frame_times[t] = (t*window_step + (n_fft-1)/2.0f) * sample_rate;
        }

        // Precompute multipliers for all frames [[7]]
        std::vector<float> multipliers(n_frames);
        std::transform(frame_times.begin(), frame_times.end(), multipliers.begin(),
                      compute_multiplier);

        for (size_t t = 0; t < n_frames; ++t) {
            float threshold = global_median * multipliers[t];
            for (size_t f = 0; f < n_freq; ++f) {
                size_t idx = t*n_freq + f;
                float original_amp = std::abs(stft_data[i][idx]);

                if (original_amp > threshold) {
                    // Optimized neighbor window collection [[8]]
                    real_vector window(2*trace_aperture + 1);
                    size_t start = (i > trace_aperture) ? i - trace_aperture : 0;
                    size_t end = std::min(i + trace_aperture + 1, n_traces);
                    
                    for (size_t j = start, k = 0; j < end; ++j, ++k) {
                        window[k] = std::abs(stft_data[j][idx]);
                    }
                    
                    float median_amp = compute_median(window);
                    filtered[i][idx] = std::polar(median_amp, phases[i][idx]);
                }
            }
        }
    }

    return filtered;
}

} // namespace tfd
