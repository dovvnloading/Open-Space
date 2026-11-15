#!/usr/bin/env bash
set -euo pipefail
build_dir="${1:?Usage: bash build_acceleration.sh /absolute/path/to/build}"
mkdir -p "$build_dir"
build_dir="$(cd "$build_dir" && pwd)"
script_dir="$(cd "$(dirname "$0")" && pwd)"
bash "$script_dir/../orbits/build.sh" "$build_dir"
if ! grep -q 'getenv("FO_ACCEL")' "$build_dir/find_orb/elem_out.cpp"; then
    patch --directory="$build_dir/find_orb" --strip=1 < "$script_dir/acceleration.patch"
fi
make -C "$build_dir/find_orb" -j4 fo PREFIX="$build_dir/prefix"
