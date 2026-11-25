#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
work=$(mktemp -d "$PWD/.lean-check.XXXXXX")
trap 'rm -rf "$work"' EXIT
mkdir -p "$work/OpenSpace"
python3 generate.py "$work/OpenSpace/Evidence.lean"
lean -DwarningAsError=true -o "$work/OpenSpace/Foundations.olean" OpenSpace/Foundations.lean
LEAN_PATH="$work" lean -DwarningAsError=true --root="$work" "$work/OpenSpace/Evidence.lean"
echo "Lean checked foundations and the generated exact evidence certificates."
