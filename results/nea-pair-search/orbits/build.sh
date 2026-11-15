#!/usr/bin/env bash
# Build a forward-only Find_Orb executable and fetch the original DE440 coefficients.
set -euo pipefail
build_dir="${1:?Usage: bash build.sh /path/to/build}"
mkdir -p "$build_dir"
build_dir="$(cd "$build_dir" && pwd)"
script_dir="$(cd "$(dirname "$0")" && pwd)"
prefix="$build_dir/prefix"
mkdir -p "$prefix"/{lib,bin,include}
cd "$build_dir"
fetch_source() {
    local project="$1" revision="$2"
    if [[ ! -d "$project" ]]; then
        local staging
        staging="$(mktemp -d -p "$build_dir")"
        curl --fail --silent --show-error --location "https://codeload.github.com/Bill-Gray/$project/tar.gz/$revision" |
            tar -xz -C "$staging" --strip-components=1
        mv "$staging" "$project"
    fi
}
fetch_source lunar 110ab53d46e3e5f1952318360e7bb01cb4884674
fetch_source jpl_eph a73f25e54d02b99b1c0d9a9d6c61acfbb2fa3a26
fetch_source sat_code 2814ea074e814d0f69b4ec1aab03425ea7a6d342
fetch_source find_orb 1de55cca12e7a585c9c7384acbb80538ea2afc3c
if ! grep -q FO_FORWARD_MODEL find_orb/elem_out.cpp; then
    patch --directory=find_orb --strip=1 < "$script_dir/forward.patch"
fi
cp lunar/*.h "$prefix/include/"
make -C jpl_eph -j4 libjpl.a PREFIX="$prefix"
make -C jpl_eph install PREFIX="$prefix"
make -C lunar -j4 install PREFIX="$prefix"
make -C sat_code -j4 libsatell.a PREFIX="$prefix"
make -C sat_code install_lib PREFIX="$prefix"
make -C find_orb -j4 fo PREFIX="$prefix"
if [[ ! -s de440.bin ]]; then
    curl --fail --silent --show-error --location \
        https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de440/linux_p1550p2650.440 -o de440.partial
    mv de440.partial de440.bin
fi
echo '29915576d0a6555766b99485ac3056ee415e86df4fce282611c31afb329ad062  de440.bin' | sha256sum --check
