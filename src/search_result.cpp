#include "search_result.h"

#include <chrono>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace {

std::string csv_escape(const std::string& value) {
    if (value.find_first_of(",\"\r\n") == std::string::npos) return value;
    std::string escaped = "\"";
    for (char ch : value) {
        if (ch == '\"') escaped += '\"';
        escaped += ch;
    }
    escaped += '\"';
    return escaped;
}

std::string timestamp_now() {
    const auto now = std::chrono::system_clock::now();
    const std::time_t value = std::chrono::system_clock::to_time_t(now);
    std::tm local_time{};
#ifdef _WIN32
    localtime_s(&local_time, &value);
#else
    localtime_r(&value, &local_time);
#endif
    std::ostringstream stream;
    stream << std::put_time(&local_time, "%Y-%m-%dT%H:%M:%S");
    return stream.str();
}

const char* csv_header() {
    return "dataset,normalization,metric,m,w,xi,method,lb_time,total_time,valid_candidates,"
           "invalid_candidates,pruned_candidates,pruning_rate,survivors,matches,"
           "lbip_time,lbip_pruning_rate,timestamp";
}

std::string csv_row(const ExperimentResult& r) {
    const double pruning_rate = r.valid_candidates == 0
        ? 0.0
        : static_cast<double>(r.pruned_candidates) /
          static_cast<double>(r.valid_candidates);
    std::ostringstream out;
    out << std::setprecision(17)
        << csv_escape(r.dataset) << ',' << csv_escape(r.normalization) << ','
        << csv_escape(r.metric) << ','
        << r.m << ',' << r.w << ',' << r.xi << ',' << csv_escape(r.method) << ','
        << r.lb_time_seconds << ',';
    if (r.has_total_time) out << r.total_time_seconds;
    out << ',' << r.valid_candidates << ',' << r.invalid_candidates << ','
        << r.pruned_candidates << ',' << pruning_rate << ','
        << r.survivors << ',';
    if (r.has_matches) out << r.matches;
    out << ',';
    if (r.has_lbip) out << r.lbip_time_seconds;
    out << ',';
    if (r.has_lbip) out << r.lbip_pruning_rate;
    out << ',' << timestamp_now();
    return out.str();
}

} // namespace

std::string infer_dataset_name(const char* target_path) {
    if (target_path == nullptr || *target_path == '\0') return "unknown";
    std::filesystem::path path(target_path);
    const std::filesystem::path parent = path.parent_path();
    return parent.empty() ? path.stem().string() : parent.filename().string();
}

void emit_experiment_result(const ExperimentResult& result) {
    const std::string row = csv_row(result);
    std::cout << "EXPERIMENT_CSV_HEADER=" << csv_header() << '\n';
    std::cout << "EXPERIMENT_CSV=" << row << std::endl;

    const char* configured_path = std::getenv("FSMELD_RESULTS_CSV");
    const std::filesystem::path csv_path =
        configured_path != nullptr && *configured_path != '\0'
            ? std::filesystem::path(configured_path)
            : std::filesystem::path("experiment_results.csv");

    std::error_code error;
    const bool needs_header = !std::filesystem::exists(csv_path, error) ||
                              std::filesystem::file_size(csv_path, error) == 0;
    std::ofstream output(csv_path, std::ios::app);
    if (!output) {
        std::cerr << "[Warning] Cannot append experiment CSV: "
                  << csv_path.string() << std::endl;
        return;
    }
    if (needs_header) output << csv_header() << '\n';
    output << row << '\n';
    std::cout << "Experiment CSV saved to: " << csv_path.string() << std::endl;
}
