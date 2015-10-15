#!/usr/bin/env bash

# Current PID: $$
# Last PID: $!

if [[ -z "$MAX_NPROC" ]]; then
  [[ $(hash nproc 2>/dev/null) ]] && MAX_NPROC=$(nproc) || MAX_NPROC=1
fi

# process queue
function waitForProcess {
  while [[ $(jobs -pr | wc -l) -ge $MAX_NPROC ]]; do
    sleep 1
  done
}
