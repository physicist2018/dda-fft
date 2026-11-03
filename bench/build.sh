#!/bin/sh
# Build the standalone FFT-kernel harness.  It is not part of dda_fft.out:
# `compile` globs *.cpp in the project root, so the harness lives in bench/.
set -e
cd "$(dirname "$0")/.."
g++ -O3 -std=c++17 -I. bench/fft_bench.cpp fft_func.cpp complex.cpp cvector.cpp -o bench/fft_bench
echo "built bench/fft_bench"
