#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
build_dir="$(mktemp -d)"
trap 'rm -rf "${build_dir}"' EXIT

"${CXX:-g++}" \
    -std=c++11 -Wall -Wextra -fsanitize=address,undefined \
    -I"${repo_dir}/source" \
    "${repo_dir}/extras/tests/testSuffixComparator.cpp" \
    "${repo_dir}/source/funCompareUintAndSuffixesMemcmp.cpp" \
    -o "${build_dir}/testSuffixComparator"

ASAN_OPTIONS=detect_leaks=0 "${build_dir}/testSuffixComparator"
