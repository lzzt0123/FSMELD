# FSMELD

FSMELD is a C++17 framework for time-series subsequence matching. It supports five elastic measures: DTW, ERP, EDR, LCSS, and MSM.

## Overview

Exact elastic distances are usually evaluated by dynamic programming and can
require quadratic time for every candidate subsequence. This cost becomes the
main bottleneck when a query must be compared with every sliding window of a
long target series.

FSMELD introduces an inner-product lower bound, denoted as `LBIP`. The bound
separates target-side envelope morphology from query-side amplitude penalties
and rewrites lower-bound evaluation as two sliding inner products. These inner
products are evaluated in batches with FFT, reducing the amortized lower-bound
cost to `O(log m)` per candidate, where `m` is the query length.

## Search pipeline

For each target segment, FSMELD performs the following stages:

1. Z-normalize the query and candidate subsequences.
2. Construct upper and lower envelopes under window constraint `w`.
3. Select upper and lower thresholds (`tU` and `tL`) for `LBIP`.
4. Transform the two bound components into sliding inner products and evaluate
   them with FFT.
5. Prune candidates whose `LBIP` value exceeds the active search threshold.
6. Apply the stronger BGLB filter to the remaining candidates.
7. Evaluate the exact elastic distance for final verification.

The implementation processes the target in FFT-friendly segments, avoiding an
index-building stage and supporting direct search over raw time-series data.

## Supported elastic measures

| Metric                          | Tag    | Meaning of `arg`             |
| ------------------------------- | ------ | ---------------------------- |
| Dynamic Time Warping            | `dtw`  | Unused; pass `0`             |
| Edit Distance with Real Penalty | `erp`  | Gap reference value `g`      |
| Edit Distance on Real Sequences | `edr`  | Matching tolerance `epsilon` |
| Longest Common Subsequence      | `lcss` | Matching tolerance `epsilon` |
| Move-Split-Merge                | `msm`  | Split/merge cost `c`         |

## Requirements

- Windows 64-bit
- CMake 3.16 or newer
- A C++17 compiler (the included FFTW import library is tested with MinGW)
- Ninja or another CMake-supported build system

The repository includes the double-precision FFTW import library and runtime
DLL in `lib/`.

## Build

Using Ninja and MinGW:

```powershell
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --target fsmeld
```

The executable is generated as `build/fsmeld.exe`. CMake also copies the FFTW
runtime DLL beside the executable.

## Input files

- `BI_T`: target series stored as raw little-endian binary `double` values.
- `Q`: query series stored as whitespace-separated text `double` values.

## Command line

```text
fsmeld BI_T Q type m epsilon metric w arg [top_k] [target_index] [exclude_start] [exclude_end]
```

| Argument  | Description                                                                  |
| --------- | ---------------------------------------------------------------------------- |
| `BI_T`    | Path to the binary target series                                             |
| `Q`       | Path to the text query series                                                |
| `type`    | Experiment or search mode listed below                                       |
| `m`       | Query and candidate subsequence length                                       |
| `epsilon` | Search threshold `xi`; this is separate from the EDR/LCSS tolerance in `arg` |
| `metric`  | One of `dtw`, `erp`, `edr`, `lcss`, or `msm`                                 |
| `w`       | Temporal constraint/envelope window size                                     |
| `arg`     | Metric-specific parameter from the table above                               |
| `top_k`   | Optional Top-K value; default is `10`                                        |

Optional arguments are positional. To provide `target_index`, first provide a
value for `top_k`.

### Modes

| Type | Operation                                                         |
| ---- | ----------------------------------------------------------------- |
| `1`  | Compute `LBIP / ELD` tightness over all valid target windows      |
| `2`  | Compute the exact Top-K elastic-distance threshold                |
| `3`  | Run the complete `LBIP -> BGLB -> exact distance` FSMELD pipeline |
| `4`  | Run the ablation pipeline `BGLB -> exact distance` without `LBIP` |
| `6`  | Run FSMELD while dynamically shrinking the Top-K threshold        |
| `8`  | Scan exact elastic distances without lower-bound filtering        |

Run `fsmeld --help` to display the same interface from the executable.

### Example

The following command runs complete ERP search with `m = 1024`, `w = 102`,
search threshold `0.5`, and ERP gap value `g = 0`:

```powershell
.\build\fsmeld.exe target_double.bin query.txt 3 1024 0.5 erp 102 0
```

To calculate TLB values for EDR:

```powershell
.\build\fsmeld.exe target_double.bin query.txt 1 1024 1.0 edr 102 0.1
```

## Experimental configuration in the paper

The paper evaluates ECG, TEMP, RANDOMWALK, THCHS, PLA, and PLETH datasets with:

- query lengths `m` in `{128, 256, 512, 1024, 2048, 4096}`;
- window ratios `w/m` in `{0.01, 0.02, 0.05, 0.10}`;
- DTW, ERP, EDR, LCSS, and MSM;

TLB is defined as `LBIP(X, Q) / ELD(X, Q)`, where a value closer to `1` denotes
a tighter lower bound. 

## Source layout

- `FSMELD.cpp`: search, TLB, Top-K, and ablation pipelines.
- `find_t_lp.cpp`: threshold selection for the upper and lower bound components.
- `q_to_h.cpp`, `lb.cpp`, `sip.cpp`: inner-product transformations and FFT-based
  sliding inner products.
- `glb.cpp`: BGLB filters and metric dispatch.
- `ElasticMeasures.cpp`: exact elastic-distance implementations.
- `lower_upper_lemire.cpp`: linear-time envelope construction.
