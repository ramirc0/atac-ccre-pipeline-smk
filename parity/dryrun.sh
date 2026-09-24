#!/usr/bin/env bash
# Normalized forced dry run: shell commands and job fields, sorted, without
# timestamps or reasons, so two layouts can be diffed. Args go to snakemake.
set -euo pipefail
snakemake -n -p -F --quiet progress "$@" 2>/dev/null \
    | awk '/^Reasons:/ {exit} 1' \
    | grep -v -e '^\[' -e 'reason:' -e 'jobid:' -e '^Building' -e 'host:' -e '^Provided' \
    | sed 's/[[:space:]]*$//' \
    | LC_ALL=C sort
