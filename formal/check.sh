#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
toolchain=$(cat lean-toolchain)
version=${toolchain#leanprover/lean4:v}
case "$(lean --version)" in
  "Lean (version $version,"*) ;;
  *) echo "Lean version must match lean-toolchain ($version)" >&2; exit 1 ;;
esac
work=$(mktemp -d "$PWD/.lean-check.XXXXXX")
trap 'rm -rf "$work"' EXIT
mkdir -p "$work/OpenSpace"
python3 -B -m unittest discover -s . -p 'test_generate.py' -v
python3 -B generate.py "$work/OpenSpace/Evidence.lean"
lean -DwarningAsError=true -o "$work/OpenSpace/Foundations.olean" OpenSpace/Foundations.lean | tee "$work/foundations.log"
export LEAN_PATH="$work"
lean -DwarningAsError=true --root="$work" "$work/OpenSpace/Evidence.lean" | tee "$work/evidence.log"
python3 -B audit_axioms.py "$work/foundations.log" "$work/evidence.log"
python3 -B test_negative.py "$work"
echo "Lean checked foundations and the generated exact evidence certificates."
