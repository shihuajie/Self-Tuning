#!/usr/bin/env bash
params="$@"
res=1
while [[ ! $res -eq 0 ]]; do
  # yes n | ...
  ./synth.sh $params <<< $'n\nexit(1);'
  res=$?
done
