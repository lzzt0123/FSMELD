#include "funcs.h"

#include <chrono>
#include <exception>
#include <iostream>
#include <string>

namespace {

const char* mode_name(int type) {
    switch (type) {
        case 1: return "LP/IP TLB";
        case 2: return "Exact Top-K threshold";
        case 3: return "FSMELD search";
        case 4: return "FSMELD raw search";
        case 6: return "Dynamic FSMELD Top-K threshold";
        case 8: return "Exact elastic-distance scan";
        default: return "Unknown";
    }
}

bool supported_mode(int type) {
    return type == 1 || type == 2 || type == 3 || type == 4 || type == 6 || type == 8;
}

void print_usage(const char* executable) {
    std::cout
        << "FSMELD elastic subsequence search\n\n"
        << "Usage:\n"
        << "  " << executable
        << " BI_T Q type m epsilon metric w arg [top_k] [target_index]"
           " [exclude_start] [exclude_end]\n\n"
        << "Modes:\n"
        << "  1  Compute LP/IP TLB for all valid target windows\n"
        << "  2  Find the exact Top-K elastic-distance threshold\n"
        << "  3  Run FSMELD (LP/IP -> BGLB -> exact distance)\n"
        << "  4  Run FSMELD raw (BGLB -> exact distance)\n"
        << "  6  Run FSMELD while dynamically shrinking the Top-K threshold\n"
        << "  8  Scan exact elastic distance without lower bounds\n\n"
        << "Metrics:\n"
        << "  dtw, erp, edr, lcss, msm\n\n"
        << "Input formats:\n"
        << "  BI_T  Raw little-endian binary doubles without a header\n"
        << "  Q     Whitespace-separated text doubles\n";
}

}  

int main(int argc, char* argv[]) {
    if (argc == 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")) {
        print_usage(argv[0]);
        return 0;
    }
    if (argc < 9) {
        print_usage(argv[0]);
        return 1;
    }

    try {
        const std::string target_path = argv[1];
        const std::string query_path = argv[2];
        const int type = std::stoi(argv[3]);
        const int m = std::stoi(argv[4]);
        const double epsilon = std::stod(argv[5]);
        std::string metric = argv[6];
        const int w = std::stoi(argv[7]);
        const double arg = std::stod(argv[8]);
        const int top_k = argc >= 10 ? std::stoi(argv[9]) : 10;
        const double target_index = argc >= 11 ? std::stod(argv[10]) : -1.0;
        const long long exclude_start = argc >= 12 ? std::stoll(argv[11]) : -1;
        const long long exclude_end = argc >= 13 ? std::stoll(argv[12]) : -1;

        if (!supported_mode(type)) {
            std::cerr << "Unsupported mode: " << type << std::endl;
            return 1;
        }
        if (metric_from_tag(metric.c_str()) == METRIC_UNKNOWN) {
            std::cerr << "Unsupported metric: " << metric << std::endl;
            return 1;
        }

        std::cout << "Mode: " << type << " (" << mode_name(type) << ")\n"
                  << "Target: " << target_path << "\n"
                  << "Query: " << query_path << "\n"
                  << "m=" << m << ", epsilon=" << epsilon
                  << ", metric=" << metric << ", w=" << w
                  << ", arg=" << arg << ", top_k=" << top_k << std::endl;

        const auto start = std::chrono::high_resolution_clock::now();
        char* metric_tag = metric.data();

        switch (type) {
            case 1:
                FSMELD_get_tlb(target_path.c_str(), query_path.c_str(), m, epsilon,
                               metric_tag, w, arg, target_index, exclude_start, exclude_end);
                break;
            case 2:
                FSMELD_find_threshold(target_path.c_str(), query_path.c_str(), m, epsilon,
                                      metric_tag, w, arg, top_k, exclude_start, exclude_end);
                break;
            case 3:
                FSMELD(target_path.c_str(), query_path.c_str(), m, epsilon,
                       metric_tag, w, arg, exclude_start, exclude_end);
                break;
            case 4:
                FSMELD_raw(target_path.c_str(), query_path.c_str(), m, epsilon,
                           metric_tag, w, arg, exclude_start, exclude_end);
                break;
            case 6:
                FSMELD_epsilon(target_path.c_str(), query_path.c_str(), m, epsilon,
                               metric_tag, w, arg, top_k, exclude_start, exclude_end);
                break;
            case 8:
                FSMELD_exact(target_path.c_str(), query_path.c_str(), m, epsilon,
                             metric_tag, w, arg, top_k, exclude_start, exclude_end);
                break;
        }

        const auto end = std::chrono::high_resolution_clock::now();
        const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << "Time taken: " << duration.count() / 1'000'000.0 << " s" << std::endl;
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Invalid arguments: " << error.what() << std::endl;
        return 1;
    }
}
