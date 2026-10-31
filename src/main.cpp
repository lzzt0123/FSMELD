#include "funcs.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
void usage() {
    std::cout << "FSMELD: exact z-normalized subsequence search\n\n"
        "Usage:\n"
        "  fsmeld search TARGET.bin QUERY.txt m cutoff metric w arg [exclude_start exclude_end]\n"
        "  fsmeld cutoff TARGET.bin QUERY.txt m k      metric w arg [exclude_start exclude_end]\n"
        "  fsmeld exact  TARGET.bin QUERY.txt m cutoff metric w arg [exclude_start exclude_end]\n\n"
        "Metrics: dtw, erp, edr, lcss, msm\n"
        "arg: 0 for DTW; gap g for ERP; match epsilon for EDR/LCSS; cost c for MSM.\n"
        "TARGET.bin contains native-endian IEEE-754 doubles; QUERY.txt contains text numbers.\n"
        "Query and candidate windows are z-normalized. Matches use zero-based indices.\n"
        "Optional exclusion is a half-open interval [start,end) of candidate starts.\n"
        "No exclusion is applied by default (self-matches are included).\n";
}
long long integer(const char* text, const char* name) {
    size_t end = 0;
    const auto value = std::stoll(text, &end);
    if (text[end] != '\0') throw std::invalid_argument(std::string("Invalid ") + name);
    return value;
}
double number(const char* text, const char* name) {
    size_t end = 0;
    const double value = std::stod(text, &end);
    if (text[end] != '\0' || !std::isfinite(value))
        throw std::invalid_argument(std::string("Invalid finite ") + name);
    return value;
}
}

int main(int argc, char** argv) {
    if (argc == 1 || (argc == 2 && (std::string(argv[1]) == "--help" ||
                                  std::string(argv[1]) == "-h"))) {
        usage(); return 0;
    }
    try {
        if (argc != 9 && argc != 11) {
            usage(); throw std::invalid_argument("Expected eight arguments, plus optional exclusion bounds.");
        }
        const std::string command = argv[1];
        if (command != "search" && command != "cutoff" && command != "exact")
            throw std::invalid_argument("Unknown command: " + command);
        const auto length = integer(argv[4], "m");
        if (length < 2 || length > std::numeric_limits<int>::max()/8)
            throw std::invalid_argument("m must be at least 2 and fit the FFT indexing range.");
        const int m = static_cast<int>(length);
        const std::string metric = argv[6];
        if (metric != "dtw" && metric != "erp" && metric != "edr" &&
            metric != "lcss" && metric != "msm")
            throw std::invalid_argument("Unsupported metric: " + metric);
        const auto radius = integer(argv[7], "w");
        if (radius < 0 || radius >= m) throw std::invalid_argument("w must satisfy 0 <= w < m.");
        const int w = static_cast<int>(radius);
        const double arg = number(argv[8], "metric parameter");
        if ((metric == "edr" || metric == "lcss" || metric == "msm") && arg < 0)
            throw std::invalid_argument("Matching tolerance and MSM cost must be non-negative.");
        const auto bytes = std::filesystem::file_size(argv[2]);
        if (bytes % sizeof(double) || bytes/sizeof(double) < 2ULL*m-1)
            throw std::invalid_argument("Target must contain at least 2*m-1 complete double values.");
        std::ifstream target(argv[2], std::ios::binary);
        if (!target) throw std::invalid_argument("Cannot open target file.");
        std::ifstream query(argv[3]);
        double sum = 0, squares = 0;
        for (int i = 0; i < m; ++i) {
            double value;
            if (!(query >> value) || !std::isfinite(value))
                throw std::invalid_argument("Query must begin with m finite numeric values.");
            sum += value; squares += value*value;
        }
        const double mean = sum/m;
        const double sigma = std::sqrt(squares/m-mean*mean);
        if (!std::isfinite(sigma) || sigma <= 1e-12)
            throw std::invalid_argument("Query cannot be z-normalized: non-finite or near-zero standard deviation.");
        long long begin = -1, end = -1;
        if (argc == 11) {
            begin = integer(argv[9], "exclude_start"); end = integer(argv[10], "exclude_end");
            if (begin < 0 || end <= begin)
                throw std::invalid_argument("Exclusion requires 0 <= start < end.");
        }
        set_normalization_mode(NormalizationMode::ZScore);
        if (command == "cutoff") {
            const auto k = integer(argv[5], "k");
            if (k < 1 || k > std::numeric_limits<int>::max() ||
                static_cast<unsigned long long>(k) > bytes/sizeof(double)-m+1)
                throw std::invalid_argument("k exceeds the candidate count or supported integer range.");
            FSMELD_find_threshold(argv[2], argv[3], m, std::numeric_limits<double>::max(),
                                  argv[6], w, arg, static_cast<int>(k), begin, end);
        } else {
            const double cutoff = number(argv[5], "cutoff");
            if (cutoff < 0) throw std::invalid_argument("cutoff must be non-negative.");
            if (command == "search")
                FSMELD(argv[2], argv[3], m, cutoff, argv[6], w, arg, begin, end);
            else
                FSMELD_exact(argv[2], argv[3], m, cutoff, argv[6], w, arg, 0, begin, end);
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
