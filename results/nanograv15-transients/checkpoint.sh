#!/usr/bin/env bash
# Reuse only unexpired primary-search artifacts from identical science code.
set -euo pipefail
study=results/nanograv15-transients
prior_run=${SCAN_RUN:-}
if [[ ! $prior_run =~ ^[0-9]+$ ]]; then
  exit 0
fi
api="repos/$GITHUB_REPOSITORY/actions/runs/$prior_run"
run_json=$(gh api "$api") || exit 0
workflow=$(jq -r .path <<<"$run_json")
prior_sha=$(jq -r .head_sha <<<"$run_json")
if [[ $workflow != .github/workflows/nanograv15-transients.yml || ! $prior_sha =~ ^[a-f0-9]{40}$ ]]; then
  exit 0
fi
completed=$(gh api "$api/jobs?filter=latest&per_page=100" --jq '.jobs[] | select(.conclusion == "success") | .name') || exit 0
for required in fetch common "reduce (0)" "reduce (1)" "reduce (2)" "reduce (3)" "individual (0)" "individual (1)" "individual (2)" "individual (3)"; do
  grep -Fxq "$required" <<<"$completed" || exit 0
done
git fetch --quiet --no-tags --depth=1 origin "$prior_sha" || exit 0
science_commands() {
  git show "$1:$workflow" | sed -nE '/python "\$STUDY\/(prepare|search|array_search|check_single|check_array)\.py"/p'
}
if [[ $(science_commands "$prior_sha") != "$(science_commands HEAD)" ]]; then
  echo "Science commands differ; computing the primary search."
  exit 0
fi
available=$(gh api --paginate "$api/artifacts?per_page=100" --jq '.artifacts[] | select(.expired == false) | .name') || exit 0
for required in input-receipts common-evidence residuals-{0..3} individual-evidence-{0..3}; do
  if ! grep -Fxq "$required" <<<"$available"; then
    echo "Missing $required; computing the primary search."
    exit 0
  fi
done
if git diff --quiet "$prior_sha" HEAD -- "$study"/{prepare.py,requirements.txt}; then
  echo "ng_run=$prior_run" >> "$GITHUB_OUTPUT"
fi
if ! git diff --quiet "$prior_sha" HEAD -- "$study"/{prepare.py,epta.py,epta_modern.py,check_epta.py,search.py,array_search.py,requirements.txt,check_array.py,check_single.py}; then
  echo "Primary or external analysis code differs; rerunning the affected stages."
  exit 0
fi
echo "Reusing verified primary-search artifacts from $prior_run at $prior_sha"
echo "scan_run=$prior_run" >> "$GITHUB_OUTPUT"
echo "scan_sha=$prior_sha" >> "$GITHUB_OUTPUT"
