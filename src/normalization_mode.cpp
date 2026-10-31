#include "funcs.h"

namespace {
thread_local NormalizationMode active_mode = NormalizationMode::ZScore;
}

void set_normalization_mode(NormalizationMode mode) {
    active_mode = mode;
}

NormalizationMode normalization_mode() {
    return active_mode;
}

const char* normalization_mode_name() {
    return active_mode == NormalizationMode::ZScore ? "zscore" : "raw";
}

bool uses_z_normalization() {
    return active_mode == NormalizationMode::ZScore;
}
