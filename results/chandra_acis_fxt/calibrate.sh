#!/usr/bin/env bash
set -euo pipefail
: "${CIAO_PREFIX:?CIAO_PREFIX must name a calibration environment}"
# Install only on workers with surviving sources. Keep the host Python separate.
conda create --yes --prefix "$CIAO_PREFIX" --override-channels \
  -c https://cxc.cfa.harvard.edu/conda/ciao -c conda-forge \
  python=3.12 ciao=4.18 caldb_main=4.12.4 numpy=2.3.5
