#!/usr/bin/env bash

# Call with:
# $ ./while_synth.sh [params ...]
#
# Resume with:
# $ directory=path_to_my_dir ./while_synth.sh
#
# Variables:
#   directory       = resuming from that directory base
#   output_prefix   = prefix for the output path
#   seed            = starting with that seed number
#   suffix          = directory suffix (timestamp by default)

read -r -d '' USAGE << EOF
A script for texture synthesis in an infinite loop.
Usage: `basename $0` [-h] [-l listfile] [-e] [-r resume_dir] [-n seed] [-p output_prefix] [-o output_dir] [-j nb_jobs] synth_params...
Usage: `basename $0` [synth_params...]
Usage: `basename $0` resume_dir

    -h                Shows this help
    -l listfile       Synthesis only data in listfile
    -e                Equivalent for -l essential.dat
    -1                Does only one full synthesis and then exit
    -r resume_dir     Resume a synthesis from a directory
    -n seed           Set the starting RNG seed [0]
    -p output_prefix  Set the output path prefix [Results]
    -o output_dir     Set the output directory [\$output_prefix/pm-sz\$size]
    -j nb_jobs        Set number of simultanious jobs [1]
EOF

# we need arguments
if [[ $# -lt 1 ]]; then
  echo "$USAGE" >&2
  exit 1
fi

# parsing the arguments
while getopts e1l:r:n:p:o:j:h OPT; do
  case $OPT in
  h)  echo "$USAGE"
      exit 0 ;;
  e)  LIST_FILE=essential.dat ;& # pass-through
  l)  if [[ -z "$LIST_FILE" ]]; then
        LIST_FILE="$OPTARG"
      fi
      if [[ ! -f "$LIST_FILE" ]]; then
        echo "List file $LIST_FILE does not exist!" >& 2
        exit 1
      else
        echo "Using file $LIST_FILE as list."
      fi
      ;;
  1)  RUN_ONCE=1 ;;
  r)  RESUME="$OPTARG" ;;
  n)  SEED=$OPTARG ;;
  p)  OUTPUT_PREFIX="$OPTARG" ;;
  o)  BASE_DIR="$OPTARG" ;;
  j)  MAX_NPROC=$OPTARG ;;
  \?) # invalid parameter
      echo "$USAGE" >&2
      exit 1 ;;
  esac
done
shift $((OPTIND-1))

# if we kill this process, children should all die!
trap "kill -- -$BASHPID" EXIT

# default arguments
[[ -z "$MAX_NPROC" ]] && MAX_NPROC=1
[[ -z "$RUN_ONCE" ]] && RUN_ONCE=0

# parse the positional parameters
if [[ $# -gt 0 ]]; then
  # check whether first argument is a resume directory
  if [[ -z "$RESUME" ]] && [[ -d "$1" ]]; then
    RESUME="$1"
    shift
  fi
fi

# should we create a new directory or resume a process?
if [[ -z "$RESUME" ]]; then
  # create a new directory!
  if [[ -z "$BASE_DIR" ]]; then
    [[ -z "$OUTPUT_PREFIX" ]] && OUTPUT_PREFIX=Results
    [[ -z "$SUFFIX" ]] && SUFFIX=$(date +"%s")
    DIR_NAME=synth-$SUFFIX
    BASE_DIR="$OUTPUT_PREFIX"/"$DIR_NAME"
  fi
  echo "New synthesis in $BASE_DIR."
  PARAMS="$@"
else
  echo "Resuming synthesis in $RESUME."
  BASE_DIR="$RESUME"
  DIR_NAME=$(basename "$RESUME")
  # compute the size if not set
  PARAMS=""
  for d in $(ls -d "$RESUME"/*/ --); do
    if [[ -f "$d"/params.dat ]]; then
      PARAMS=$(cat "$d"/params.dat)
    fi
  done
  if [[ -z "$PARAMS" ]] && [[ -f "$RESUME"/params.dat ]]; then
    PARAMS=$(cat "$RESUME"/params.dat)
  fi
  echo "Params: $PARAMS"
fi

# seed
[[ -z "$SEED" ]] && SEED=0
i="$SEED"

# parallelization functions
source func_parallel.sh

# output pid in case we need to get killed and it doesn't work properly
echo "PID: $$"

if [[ -z "$LIST_FILE" ]]; then
  FILE_LIST=$(ls ../../resources/*/*_S.jpg)
else
  FILE_LIST=$(cat "$LIST_FILE")
fi

# record when we start
start_time=$(date +%s)

# loop
RUNNING=1
while [[ "$RUNNING" -eq 1 ]]; do
  # where do we output? in seed_X?
  if [[ "$RUN_ONCE" -eq 1 ]]; then
    OUTPUT_DIR="$BASE_DIR"
    RUNNING=0 # we do not run once more
  else
    OUTPUT_DIR="$BASE_DIR"/seed_$i
  fi
  mkdir -p $OUTPUT_DIR
  echo "$PARAMS" > "$OUTPUT_DIR"/params.dat
  echo "Parameters in $OUTPUT_DIR/params.dat"
  for file in $FILE_LIST; do
    RES_DIR=$(basename $file)
    RES_DIR=${RES_DIR/.jpg/}
    TEX_DIR="$OUTPUT_DIR"/"$RES_DIR"
    if [[ -d "$TEX_DIR" ]]; then
      if [[ -f "$TEX_DIR"/done ]]; then
        echo "Skipping $RES_DIR"
        continue
      else
        echo "Resuming $RES_DIR"
      fi
    else
      echo "Synth $RES_DIR"
    fi
    # do the job
    COMMAND="./checked_synth.sh $file $OUTPUT_DIR rand_seed $i num_threads 1 $PARAMS"
    if [[ "$MAX_NPROC" -le 1 ]]; then
      $COMMAND
    else
      waitForProcess
      mkdir -p "$TEX_DIR"
      $COMMAND 2> "$TEX_DIR"/proc_err.log > "$TEX_DIR"/proc.log &
      queueProcess $!
    fi
  done
  # new seed
  i=$(($i + 1))
done

# wait for all the chidlren
echo "Waiting for the last synthesis"
wait
end_time=$(date +%s)
duration=$(($end_time - $start_time))
echo "Done. ($duration sec)"

# store time result
echo $duration > "$OUTPUT_DIR"/duration.dat
