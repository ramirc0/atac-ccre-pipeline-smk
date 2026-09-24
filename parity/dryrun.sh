#!/usr/bin/env bash
# Normalized dry run: shell commands and job fields, sorted, without timestamps
# or reason lines, so two layouts can be diffed. Args go to snakemake.
set -euo pipefail
snakemake -n -p --quiet progress "$@" 2>/dev/null \
    | grep -v -e '^\[' -e 'reason:' -e 'jobid:' -e '^Building' -e 'host:' -e '^Provided' \
    | sed 's/[[:space:]]*$//' \
    | LC_ALL=C sort
