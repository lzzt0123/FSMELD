#ifndef FSMELD_EXPERIMENT_RESULT_H
#define FSMELD_EXPERIMENT_RESULT_H

#include <cstdint>
#include <string>

struct ExperimentResult {
    std::string dataset;
    std::string normalization;
    std::string metric;
    int m = 0;
    int w = 0;
    double xi = 0.0;
    std::string method;
    double lb_time_seconds = 0.0;
    bool has_total_time = true;
    double total_time_seconds = 0.0;
    std::uint64_t valid_candidates = 0;
    std::uint64_t invalid_candidates = 0;
    std::uint64_t pruned_candidates = 0;
    std::uint64_t survivors = 0;
    bool has_matches = true;
    std::uint64_t matches = 0;
    bool has_lbip = false;
    double lbip_time_seconds = 0.0;
    double lbip_pruning_rate = 0.0;
};

std::string infer_dataset_name(const char* target_path);
void emit_experiment_result(const ExperimentResult& result);

#endif
