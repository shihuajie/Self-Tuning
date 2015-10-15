#!/usr/bin/env bash
DEBUG_PARAMS='-DSAFE_MAT=1'
if [ $# -ge 1 ]; then
  DEBUG_PARAMS="$@ $DEBUG_PARAMS"
fi
./mex_build.sh $DEBUG_PARAMS
