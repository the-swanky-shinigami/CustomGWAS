#!/usr/bin/env bash
set -euo pipefail

export R_EXECUTABLE="${R_EXECUTABLE:-/usr/bin/Rscript}"
export JUPYTER_TOKEN="${JUPYTER_TOKEN:-customgwas}"

if [[ "$#" -eq 0 ]]; then
  exec jupyter lab --ip=0.0.0.0 --no-browser --allow-root \
    --ServerApp.token="$JUPYTER_TOKEN" \
    --ServerApp.root_dir=/work
fi

exec "$@"
