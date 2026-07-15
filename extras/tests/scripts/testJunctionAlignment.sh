#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
build_dir="$(mktemp -d)"
trap 'rm -rf "${build_dir}"' EXIT

"${CXX:-g++}" \
    -std=c++11 -Wall -Wextra \
    -I"${repo_dir}/source" \
    "${repo_dir}/extras/tests/testJunctionAlignment.cpp" \
    -o "${build_dir}/testJunctionAlignment"

"${build_dir}/testJunctionAlignment"
