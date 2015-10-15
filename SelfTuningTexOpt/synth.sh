#!/usr/bin/env bash

# input checking
if [[ $# -lt 3 ]]; then
  echo "Usage: ./synth.sh source_file output_dir [param_list]"
  exit 1
fi

INPUT="'$1'"
OUTPUT="'$2'"
shift 2

param_list="$@"
param_count="$#"
if [[ $# -gt 0 ]]; then
  while (( $# )); do
    if [[ -z "$1" ]]; then
      shift
      continue
    #elif [[ $# -lt 2 ]]; then
    #  echo "Current param: $1"
    #  echo "Umbalanced parameter list ($param_count): $param_list"
    #  exit 1
    fi
    # is it a char that is not escaped?
    if [[ "$1" =~ ^[%a-zA-Z] ]]; then
      # escape parameter for matlab as char
      PARAMS="$PARAMS, '$1'"
    else
      # no escaping, it's not a string, or it's already escaped
      PARAMS="$PARAMS, $1"
    fi
    shift
  done
else
  PARAMS=""
fi

# the matlab function to call
[[ -z "$SYNTH_FUNC" ]] && SYNTH_FUNC="synth_func"

# we generate the matlab command
COMMAND="$SYNTH_FUNC($INPUT, $OUTPUT $PARAMS); exit"
echo "$COMMAND"
# we run it! (-nojvm disables exp_figure)
matlab -nosplash -nojvm -nodesktop -logfile remoteAutocode.log -r "$COMMAND"
