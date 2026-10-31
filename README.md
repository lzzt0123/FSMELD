# FSMELD

FSMELD is a C++ framework for exact subsequence range search under elastic
distance measures. Given a long time series, a query, and a distance cutoff,
it finds the starting positions of matching subsequences. The query and each
candidate window are z-normalized independently.

FSMELD evaluates the inner-product lower bound `LB_IP` in batches using FFT.
Candidates that survive this filter proceed to GLB filtering and exact
elastic-distance verification. This `LB_IP -> GLB -> ELD` cascade reduces the
number of expensive distance computations while preserving exact verification
of the returned matches.

[中文说明](README.zh-CN.md)

## Features

- `search`: FSMELD search with a fixed distance cutoff.
- `cutoff`: determine the exact distance of the k-th nearest valid candidate.
- `exact`: scan every valid candidate with the elastic-distance verifier.
- DTW, squared-cost ERP, EDR, normalized LCSS distance, and MSM.
- Per-window z-normalization, full blocks and tail windows, and invalid-window filtering.

## Build

Requirements:

- CMake 3.20 or later.
- GCC/MinGW-w64 or Clang with C++17 and GNU variable-length-array extensions.
  The numerical kernels use array extensions that MSVC does not support.
- FFTW3 double-precision headers and library, installed separately.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DFFTW3_ROOT=/path/to/fftw3
cmake --build build --parallel
./build/fsmeld --help
```

If FFTW3 is installed in a system search path, omit `FFTW3_ROOT`. Alternatively,
set `FFTW3_INCLUDE_DIR` to the directory containing `fftw3.h` and
`FFTW3_LIBRARY` to the double-precision library file.

On Windows, run from a MinGW terminal with its compiler tools on `PATH`:

```sh
cmake -S . -B build -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release -DFFTW3_ROOT=C:/fftw3
cmake --build build --parallel
build/fsmeld.exe --help
```

For a shared FFTW installation, put its DLL beside `fsmeld.exe` or make it
available on `PATH`. Compiler runtime dependencies must also be available.

## Input

- `target.bin`: contiguous native-endian IEEE-754 64-bit floating-point values,
  with no header. The extension does not determine the format.
- `query.txt`: whitespace-separated numeric text. The first `m` values are used.
- `m >= 2`; the target must contain at least `2*m-1` values.
- `0 <= w < m`; the cutoff must be finite and non-negative.
- Query and candidate windows are z-normalized independently. Invalid windows
  containing non-finite input, non-finite standard deviation, or standard
  deviation at most `1e-12` are skipped. An invalid query is rejected.

Do not concatenate independent recordings unless windows crossing their
boundaries are intended to be searched.

## Commands

```text
fsmeld search TARGET.bin QUERY.txt m cutoff metric w arg
fsmeld cutoff TARGET.bin QUERY.txt m k      metric w arg
fsmeld exact  TARGET.bin QUERY.txt m cutoff metric w arg
```

| Metric | `arg` | Distance convention |
|---|---|---|
| `dtw` | ignored; use `0` | Square root of accumulated squared cost |
| `erp` | gap value `g` | Square root of accumulated squared match/gap cost |
| `edr` | non-negative match tolerance | Edit count |
| `lcss` | non-negative match tolerance | `1 - LCSS/m` |
| `msm` | non-negative move/split/merge cost `c` | Accumulated cost |

### Example workflow

For a query of length 128, the following example uses ERP with `g=0` and a
temporal constraint of `w=6`:

```sh
# Determine the rank-10 cutoff for ERP.
./build/fsmeld cutoff target.bin query.txt 128 10 erp 6 0

# Use the reported Final Top-K Threshold as the search cutoff.
# Replace 0.75 with the reported cutoff or your chosen search cutoff.
./build/fsmeld search target.bin query.txt 128 0.75 erp 6 0

# Scan all candidates directly at the same cutoff.
./build/fsmeld exact target.bin query.txt 128 0.75 erp 6 0
```

`cutoff` reports the k-th smallest exact distance; tied candidates can make a
range search at that cutoff return more than k matches. It fails when fewer
than k valid candidates exist. `search` and `exact` retain only
computed exact distances `<= cutoff`.

## Output

- Match lines: `A[0] = ...`, `A[1] = ...`, giving candidate starting positions.
- `cutoff` prints `Final Top-K Threshold` and ranked candidate indices.
- Search runs print counts, timings, and a CSV summary, and append to
  `experiment_results.csv` in the working directory.
- Set `FSMELD_RESULTS_CSV` to choose a different CSV destination.
- Set `FSMELD_QUIET_MATCHES=1` to omit individual match lines.

The macro-segment length is `2^ceil(log2(8*m))`, as in the paper. It is
shortened to the target length when the target is smaller than one segment.

## Source layout

- `src/fsmeld.cpp`: normalized search, Top-K cutoff computation, exact scan,
  shared candidate validation, and tail handling.
- `src/find_t_lp.cpp`: heuristic upper/lower threshold selection.
- `src/lb.cpp`, `src/q_to_h.cpp`, `src/lb_normalize.cpp`: masks, penalties,
  and normalized envelopes.
- `src/sip.cpp`: FFT sliding-inner-product primitives.
- `src/glb.cpp`: GLB second-stage filter.
- `src/ElasticMeasures.cpp`: elastic-distance dynamic programming.
- Remaining source files provide envelopes, normalization, reporting, and CLI.

