#include <iostream>
#include <fstream>
#include <filesystem>
#include <map>
#include <vector>
#include <string>
#include <chrono>
#include <indicators/progress_bar.hpp>
#include <segyio/segy.h>
#include "tfd.h"

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <param_file>" << std::endl;
        return EXIT_FAILURE;
    }

    std::string param_file = argv[1];
    std::ifstream param_stream(param_file);
    if (!param_stream.is_open()) {
        std::cerr << "Error: Could not open parameter file " << param_file << std::endl;
        return EXIT_FAILURE;
    }

    std::map<std::string, std::string> params;
    std::string line;
    while (std::getline(param_stream, line)) {
        std::istringstream iss(line);
        std::string key, value;
        if (std::getline(iss, key, '=') && std::getline(iss, value)) {
            params[key] = value;
            std::cout << "Key: " << key << "; Value: " << value << std::endl;
        }
    }

    // Extract parameters
    size_t n_fft = std::stoul(params["n_fft"]);
    size_t n_overlap = std::stoul(params["n_overlap"]);
    size_t trace_aperture = std::stoul(params["trace_aperture"]);
    std::string input_file = params["input_file"];
    std::string output_file = params["output_file"];
    int id_offset = std::stoi(params["id_offset"]);

    // Parse threshold multipliers
    std::vector<std::pair<float, float>> threshold_multipliers;
    std::string threshold_str = params["threshold_multipliers"];
    threshold_str.erase(std::remove_if(threshold_str.begin(), threshold_str.end(), ::isspace), threshold_str.end());
    std::istringstream iss_cleaned(threshold_str);
    std::string pair_str;
    while (std::getline(iss_cleaned, pair_str, ',')) {
        std::istringstream pair_ss(pair_str);
        std::string first, second;
        if (std::getline(pair_ss, first, '-') && std::getline(pair_ss, second)) {
            threshold_multipliers.emplace_back(std::stof(first), std::stof(second));
        } else {
            std::cerr << "Error: Invalid threshold multipliers format" << std::endl;
            return EXIT_FAILURE;
        }
    }

    // Copy input to output file
    std::cout << "Copying input file to output file..." << std::endl;
    try {
        std::filesystem::copy(input_file, output_file, std::filesystem::copy_options::overwrite_existing);
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    // Open files
    segy_file* file = segy_open(input_file.c_str(), "r");
    segy_file* out_file = segy_open(output_file.c_str(), "r+");
    if (!file || !out_file) {
        std::cerr << "Failed to open files for processing." << std::endl;
        return EXIT_FAILURE;
    }

    // Read binary header
    char binheader[SEGY_BINARY_HEADER_SIZE];
    if (segy_binheader(file, binheader) != 0) {
        std::cerr << "Failed to read binary header." << std::endl;
        segy_close(file);
        return EXIT_FAILURE;
    }

    int format = segy_format(binheader);
    std::cout << "Data format code: " << format << std::endl;

    // Get trace parameters
    int trace_count;
    long trace0 = segy_trace0(binheader);
    int samples = segy_samples(binheader);
    int trace_bsize = segy_trace_bsize(samples);
    float sample_interval = 0.0f;
    segy_sample_interval(file, 0.0f, &sample_interval);
    sample_interval /= 1e6; // Convert microseconds to seconds

    if (segy_traces(file, &trace_count, trace0, trace_bsize) != 0) {
        std::cerr << "Failed to get number of traces." << std::endl;
        segy_close(file);
        return EXIT_FAILURE;
    }

    // Build trace map
    std::map<int, std::vector<size_t>> seismogram_traces;
    char trace_header[SEGY_TRACE_HEADER_SIZE];
    for (size_t i = 0; i < trace_count; ++i) {
        if (segy_traceheader(file, i, trace_header, trace0, trace_bsize) != 0) {
            std::cerr << "Failed to read trace header for trace " << i << std::endl;
            continue;
        }

        int32_t value;
        if (segy_get_field(trace_header, id_offset, &value) != 0) {
            std::cerr << "Failed to get field from trace header for trace " << i << std::endl;
            continue;
        }
        seismogram_traces[value].push_back(i);
    }

    // Initialize progress tracking 
    auto start_time = std::chrono::steady_clock::now();
    indicators::ProgressBar progress_bar{
        indicators::option::BarWidth{50},
        indicators::option::Start{"["},
        indicators::option::Fill{"="},
        indicators::option::Lead{">"},
        indicators::option::Remainder{" "},
        indicators::option::End{"]"},
        indicators::option::ForegroundColor{indicators::Color::green},
        indicators::option::ShowElapsedTime{true},
        indicators::option::ShowRemainingTime{true}
    };

    progress_bar.set_progress(0);
    progress_bar.set_option(indicators::option::MaxProgress{
        static_cast<size_t>(seismogram_traces.size())});

    size_t processed_records = 0;
    for (const auto& [record_id, trace_indices] : seismogram_traces) {
        size_t n_traces = trace_indices.size();
        tfd::seismic_data data(n_traces, tfd::real_vector(samples));

        for (size_t i = 0; i < n_traces; ++i) {
            std::vector<float> buffer(samples);
            if (segy_readtrace(file, trace_indices[i], buffer.data(), trace0, trace_bsize) != 0) {
                std::cerr << "Read error for trace " << trace_indices[i] << std::endl;
                continue;
            }

            // Convert to native format (IEEE)
            segy_to_native(format, samples, buffer.data());
            data[i].assign(buffer.begin(), buffer.end());
        }

        // Process data
        auto spectrograms = tfd::compute_stft(data, n_fft, n_overlap);
        auto filtered = tfd::tfd_noise_rejection(spectrograms, trace_aperture, threshold_multipliers, sample_interval, n_fft, n_overlap);
        auto processed_data = tfd::compute_istft(filtered, n_fft, n_overlap, samples);

        // Convert back to SEG-Y format and write
        for (size_t i = 0; i < n_traces; ++i) {
            std::vector<float> temp_buffer = processed_data[i];
            segy_from_native(format, samples, temp_buffer.data());
            
            if (segy_writetrace(out_file, trace_indices[i], temp_buffer.data(), trace0, trace_bsize) != 0) {
                std::cerr << "Write error for trace " << trace_indices[i] << std::endl;
            }
        }

        // Update progress bar [[9]][[5]]
        processed_records++;
        progress_bar.set_progress(processed_records);
    }

    // Finalize progress bar and show total duration [[8]]
    progress_bar.mark_as_completed();
    auto end_time = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(
        end_time - start_time).count();
    
    std::cout << "\nProcessing completed in " << duration 
              << " seconds. Results saved to " << output_file << std::endl;

    segy_close(file);
    segy_close(out_file);
    return 0;
}
