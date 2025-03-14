#ifndef TFD_H
#define TFD_H

#include <vector>
#include <complex>
#include <string>
#include <stdexcept>
#include <cstddef>

namespace tfd {

using complex = std::complex<float>;
using real_vector = std::vector<float>;
using complex_vector = std::vector<complex>;
using seismic_data = std::vector<real_vector>;
using complex_spectrogram = std::vector<complex_vector>;

// Оконная функция Ханна
real_vector window_hann(size_t n);

// Чтение сейсмограммы из бинарного файла
seismic_data read_seismic(const std::string& filename, size_t n_traces, size_t n_time_samples);

// Отраженное дополнение для обработки краев сигнала
real_vector reflect_padding(const real_vector& trace, size_t target_length);

// Вычисление медианы для фильтрации
float compute_median(real_vector values);

// Прямое STFT с использованием FFTW3
complex_spectrogram compute_stft(const seismic_data& seismic, size_t n_fft, size_t n_overlap);

// Обратное STFT с использованием FFTW3
seismic_data compute_istft(const complex_spectrogram& spectrograms, size_t n_fft, size_t n_overlap, size_t original_length);

// Функция подавления шума (переименованная из apply_median_filter)
complex_spectrogram tfd_noise_rejection(const complex_spectrogram& stft_data, size_t trace_aperture,
                                        const std::vector<std::pair<float, float>>& threshold_params,
                                        float sample_rate, size_t n_fft, size_t n_overlap);

}

#endif // TFD_H
