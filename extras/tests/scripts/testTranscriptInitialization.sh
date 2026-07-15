#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
build_dir="$(mktemp -d)"
trap 'rm -rf "${build_dir}"' EXIT

"${CXX:-g++}" \
    -std=c++11 -Wall -Wextra -fsanitize=address,undefined \
    -I"${repo_dir}/source" \
    "${repo_dir}/extras/tests/testTranscriptInitialization.cpp" \
    "${repo_dir}/source/Transcript.cpp" \
    -o "${build_dir}/testTranscriptInitialization"

ASAN_OPTIONS=detect_leaks=0 UBSAN_OPTIONS=halt_on_error=1 \
    "${build_dir}/testTranscriptInitialization"
