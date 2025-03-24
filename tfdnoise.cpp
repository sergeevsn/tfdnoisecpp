#include <boost/mpi.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <indicators/progress_bar.hpp>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <map>
#include <vector>
#include <string>
#include <chrono>
#include <thread>
#include <segyio/segy.h>
#include "tfd.h"
#include "SimpleIni.h"

namespace mpi = boost::mpi;

int main(int argc, char* argv[]) {
    // Initialize MPI environment
    mpi::environment env(argc, argv);
    mpi::communicator world;

    // Check command line arguments
    if (argc != 2) {
        if (world.rank() == 0) {
            std::cerr << "Usage: " << argv[0] << " <param_file>" << std::endl;
        }
        return EXIT_FAILURE;
    }

    // Load and parse configuration file
    CSimpleIniA ini;
    ini.SetUnicode();
    std::string param_file = argv[1];

    if (world.rank() == 0) {
        SI_Error rc = ini.LoadFile(param_file.c_str());
        if (rc < 0) {
            std::cerr << "Error: Could not open parameter file " << param_file << std::endl;
            return EXIT_FAILURE;
        }

        // Display loaded parameters
        CSimpleIniA::TNamesDepend keys;
        ini.GetAllKeys("", keys);
        std::cout << "Parameters from file " << param_file << ":" << std::endl;
        for (const auto& key : keys) {
            const char* value = ini.GetValue("", key.pItem);
            std::cout << key.pItem << " = " << value << std::endl;
        }
    }

    // Get parameters and broadcast to all processes
    size_t n_fft = ini.GetLongValue("", "n_fft", 0);
    size_t trace_aperture = ini.GetLongValue("", "trace_aperture", 0);
    std::string input_file = ini.GetValue("", "input_file", "");
    std::string output_file = ini.GetValue("", "output_file", "");
    int id_offset = ini.GetLongValue("", "gather_byte", 0);

    mpi::broadcast(world, n_fft, 0);
    mpi::broadcast(world, trace_aperture, 0);
    mpi::broadcast(world, input_file, 0);
    mpi::broadcast(world, output_file, 0);
    mpi::broadcast(world, id_offset, 0);

    // Validate required parameters
    if (n_fft == 0 || trace_aperture == 0 || input_file.empty() || output_file.empty()) {
        if (world.rank() == 0) {
            std::cerr << "Error: Missing required parameters" << std::endl;
        }
        return EXIT_FAILURE;
    }

    // Parse threshold multipliers
    std::vector<std::pair<float, float>> threshold_multipliers;
    if (world.rank() == 0) {
        const char* threshold_str = ini.GetValue("", "threshold_multipliers", "");
        std::string cleaned_str(threshold_str);
        cleaned_str.erase(std::remove_if(cleaned_str.begin(), cleaned_str.end(), ::isspace), 
                         cleaned_str.end());
        
        std::istringstream iss_cleaned(cleaned_str);
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
    }
    mpi::broadcast(world, threshold_multipliers, 0);

    // Copy input file to output (only rank 0)
    if (world.rank() == 0) {
        std::cout << "Copying input file to output file..." << std::endl;
        try {
            std::filesystem::copy(input_file, output_file, 
                                 std::filesystem::copy_options::overwrite_existing);
        } catch (const std::filesystem::filesystem_error& e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return EXIT_FAILURE;
        }
    }

    world.barrier();

    // Open SEG-Y files
    segy_file* file = segy_open(input_file.c_str(), "r");
    segy_file* out_file = segy_open(output_file.c_str(), "r+");
    if (!file || !out_file) {
        std::cerr << "Process " << world.rank() << ": Failed to open files for processing." << std::endl;
        return EXIT_FAILURE;
    }

    // Read binary header
    char binheader[SEGY_BINARY_HEADER_SIZE];
    if (segy_binheader(file, binheader) != 0) {
        std::cerr << "Process " << world.rank() << ": Failed to read binary header." << std::endl;
        segy_close(file);
        return EXIT_FAILURE;
    }

    // Get file format and display (rank 0 only)
    int format = segy_format(binheader);
    if (world.rank() == 0) {
        std::cout << "Data format code: " << format << std::endl;
    }

    // Get trace information
    int trace_count;
    long trace0 = segy_trace0(binheader);
    int samples = segy_samples(binheader);
    int trace_bsize = segy_trace_bsize(samples);
    float sample_interval = 0.0f;
    segy_sample_interval(file, 0.0f, &sample_interval);
    sample_interval /= 1e6;

    if (segy_traces(file, &trace_count, trace0, trace_bsize) != 0) {
        std::cerr << "Process " << world.rank() << ": Failed to get number of traces." << std::endl;
        segy_close(file);
        return EXIT_FAILURE;
    }

    // Read trace headers and group by record ID (rank 0 only)
    std::map<int, std::vector<size_t>> seismogram_traces;
    if (world.rank() == 0) {
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
    }
    mpi::broadcast(world, seismogram_traces, 0);

    // Prepare record IDs for distribution
    std::vector<int> record_ids;
    for (const auto& [id, traces] : seismogram_traces) {
        record_ids.push_back(id);
    }

    // Calculate work distribution
    size_t total_records = record_ids.size();
    size_t records_per_process = total_records / world.size();
    size_t remainder = total_records % world.size();
    size_t start_idx = world.rank() * records_per_process + std::min(static_cast<size_t>(world.rank()), remainder);
    size_t end_idx = start_idx + records_per_process + (world.rank() < remainder ? 1 : 0);
    if (end_idx > total_records) end_idx = total_records;

    auto start_time = std::chrono::steady_clock::now();

    // Initialize progress bar (rank 0 only)
    std::unique_ptr<indicators::ProgressBar> progress_bar;
    if (world.rank() == 0) {
        progress_bar = std::make_unique<indicators::ProgressBar>(
            indicators::option::BarWidth{50},
            indicators::option::Start{"["},
            indicators::option::Fill{"="},
            indicators::option::Lead{">"},
            indicators::option::Remainder{" "},
            indicators::option::End{"]"},
            indicators::option::PostfixText{"0/" + std::to_string(total_records) + " Elapsed: 0s Remaining: 0s"},
            indicators::option::ForegroundColor{indicators::Color::green},
            indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}},
            indicators::option::MaxProgress{total_records}
        );
        std::cout << "Total records: " << total_records 
                  << ", Records per process: " << records_per_process 
                  << ", Remainder: " << remainder << std::endl;
    }

    // Initialize progress tracking
    size_t local_processed_records = 0;
    size_t total_processed = 0;
    std::vector<mpi::request> send_requests;
    size_t expected_messages = 0;

    // Helper function to update progress bar with elapsed and remaining time
    auto update_progress_bar = [&](size_t current, size_t total, auto start) {
        if (progress_bar) {
            auto now = std::chrono::steady_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start).count();
            double progress_ratio = static_cast<double>(current) / total;
            double remaining = (progress_ratio > 0) ? (elapsed / progress_ratio) * (1 - progress_ratio) : 0;
            std::string postfix = std::to_string(current) + "/" + std::to_string(total) +
                                 " Elapsed: " + std::to_string(elapsed) + "s" +
                                 " Remaining: " + std::to_string(static_cast<int>(remaining)) + "s";
            progress_bar->set_option(indicators::option::PostfixText{postfix});
            progress_bar->set_progress(current);
        }
    };

    if (world.rank() == 0) {
        for (int rank = 1; rank < world.size(); ++rank) {
            size_t records_per_rank = records_per_process + (rank < remainder ? 1 : 0);
            expected_messages += records_per_rank;
        }
    }

    // Initialize non-blocking receives (rank 0 only)
    std::vector<mpi::request> recv_requests;
    std::vector<int> updates(world.size(), 0);
    if (world.rank() == 0) {
        for (int rank = 1; rank < world.size(); ++rank) {
            size_t records_per_rank = records_per_process + (rank < remainder ? 1 : 0);
            for (size_t i = 0; i < records_per_rank; ++i) {
                recv_requests.push_back(world.irecv(rank, 1, updates[rank]));
            }
        }
    }

    // Main processing loop
    for (size_t idx = start_idx; idx < end_idx; ++idx) {
        int record_id = record_ids[idx];
        const auto& trace_indices = seismogram_traces[record_id];
        size_t n_traces = trace_indices.size();

        // Read trace data
        tfd::seismic_data data(n_traces, tfd::real_vector(samples));
        for (size_t i = 0; i < n_traces; ++i) {
            std::vector<float> buffer(samples);
            if (segy_readtrace(file, trace_indices[i], buffer.data(), trace0, trace_bsize) != 0) {
                std::cerr << "Process " << world.rank() << ": Read error for trace " 
                          << trace_indices[i] << std::endl;
                continue;
            }
            segy_to_native(format, samples, buffer.data());
            data[i].assign(buffer.begin(), buffer.end());
        }

        // Process data
        auto n_overlap = n_fft / 2;
        auto spectrograms = tfd::compute_stft(data, n_fft, n_overlap);
        auto filtered = tfd::tfd_noise_rejection(
            spectrograms, 
            trace_aperture, 
            threshold_multipliers, 
            sample_interval, 
            n_fft, 
            n_overlap
        );
        auto processed_data = tfd::compute_istft(filtered, n_fft, n_overlap, samples);

        // Write processed data
        for (size_t i = 0; i < n_traces; ++i) {
            std::vector<float> temp_buffer = processed_data[i];
            segy_from_native(format, samples, temp_buffer.data());
            if (segy_writetrace(out_file, trace_indices[i], temp_buffer.data(), trace0, trace_bsize) != 0) {
                std::cerr << "Process " << world.rank() << ": Write error for trace " 
                          << trace_indices[i] << std::endl;
            }
        }

        // Update progress
        local_processed_records++;

        // Non-root processes send progress updates
        if (world.rank() != 0) {
            send_requests.push_back(world.isend(0, 1, static_cast<int>(local_processed_records)));
        }

        // Root process updates the progress bar
        if (world.rank() == 0) {
            total_processed++;
            update_progress_bar(total_processed, total_records, start_time);
            
            // Check for incoming messages from other processes
            for (auto it = recv_requests.begin(); it != recv_requests.end(); ) {
                if (it->test()) {
                    total_processed++;
                    update_progress_bar(total_processed, total_records, start_time);
                    it = recv_requests.erase(it);
                } else {
                    ++it;
                }
            }
        }
    }

    // Wait for all sends to complete (non-root processes)
    if (world.rank() != 0) {
        mpi::wait_all(send_requests.begin(), send_requests.end());
    }

    // Root process handles remaining messages and updates progress bar
    if (world.rank() == 0) {
        size_t received_messages = total_processed - (end_idx - start_idx);
        while (received_messages < expected_messages) {
            for (auto it = recv_requests.begin(); it != recv_requests.end(); ) {
                if (it->test()) {
                    total_processed++;
                    received_messages++;
                    update_progress_bar(total_processed, total_records, start_time);
                    it = recv_requests.erase(it);
                } else {
                    ++it;
                }
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }

        if (total_processed != total_records) {
            std::cerr << "Error: Total processed (" << total_processed 
                      << ") does not match total records (" << total_records << ")" << std::endl;
        }
    }

    world.barrier();

    // Final output (rank 0 only)
    if (world.rank() == 0) {
        auto end_time = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
        std::cout << "\nProcessing completed in " << duration 
                  << " seconds. Results saved to " << output_file << std::endl;
    }

    // Cleanup
    segy_close(file);
    segy_close(out_file);
    return 0;
}
