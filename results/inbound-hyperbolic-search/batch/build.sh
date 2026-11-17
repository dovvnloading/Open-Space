#!/usr/bin/env bash
set -euo pipefail
work=${1:?supply an empty build directory}
mkdir -p "$work" "$HOME/bin" "$HOME/lib" "$HOME/include"
cd "$work"
git clone -q --depth 1 https://github.com/Bill-Gray/find_orb.git find_orb
git -C find_orb fetch -q --depth 1 origin 1de55cca12e7a585c9c7384acbb80538ea2afc3c
git -C find_orb checkout -q --detach 1de55cca12e7a585c9c7384acbb80538ea2afc3c
git clone -q --depth 1 https://github.com/Bill-Gray/lunar.git lunar
git -C lunar fetch -q --depth 1 origin 110ab53d46e3e5f1952318360e7bb01cb4884674
git -C lunar checkout -q --detach 110ab53d46e3e5f1952318360e7bb01cb4884674
git clone -q --depth 1 https://github.com/Bill-Gray/sat_code.git sat_code
git -C sat_code fetch -q --depth 1 origin 2814ea074e814d0f69b4ec1aab03425ea7a6d342
git -C sat_code checkout -q --detach 2814ea074e814d0f69b4ec1aab03425ea7a6d342
git clone -q --depth 1 https://github.com/Bill-Gray/jpl_eph.git jpl_eph
git -C jpl_eph fetch -q --depth 1 origin a73f25e54d02b99b1c0d9a9d6c61acfbb2fa3a26
git -C jpl_eph checkout -q --detach a73f25e54d02b99b1c0d9a9d6c61acfbb2fa3a26
make -s -C lunar -j2 install
make -s -C jpl_eph -j2 libjpl.a
make -s -C jpl_eph install
make -s -C sat_code -j2 sat_id
make -s -C sat_code install
make -s -C find_orb -j2
make -s -C find_orb install IDIR="$HOME/.find_orb"
ephem="$HOME/.find_orb/linux_p1550p2650.430t"
curl -fsSL --retry 2 --max-time 180 -o "$ephem" https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de430t/linux_p1550p2650.430t
sha256sum "$ephem"
