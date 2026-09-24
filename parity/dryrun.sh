#!/usr/bin/env bash
# Normalized forced dry run: shell commands and job fields, sorted, without
# timestamps or reasons, so two layouts can be diffed. Args go to snakemake.
set -euo pipefail
err=$(mktemp)
trap 'rm -f "$err"' EXIT
if ! out=$(snakemake -n -p -F --quiet progress "$@" 2>"$err"); then
    cat "$err" >&2
    exit 1
fi
awk '/^Reasons:/ {skip = 1} !skip' <<< "$out" \
    | grep -v -e '^\[' -e 'reason:' -e 'jobid:' -e '^Building' -e 'host:' -e '^Provided' \
    | sed 's/[[:space:]]*$//' \
    | LC_ALL=C sort
