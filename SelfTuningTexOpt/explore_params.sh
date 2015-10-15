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

# the current command
BASE_LINE="$0 $@"

read -r -d '' USAGE << EOF
A script for the exploration of parameters in texture synthesis
Usage: `basename $0` [-h] [-l list_file|-i file] [-o output_dir] [-j nb_jobs] synth_params...

    -h                Shows this help
    -c synth_comm     Comparison synthesis method
    -l listfile       Synthesis only data in listfile
    -i input_file     Synthesis of the given file only
    -o output_dir     Set the output directory [\$output_prefix/pm-sz\$size]
    -p paramfile      Synthesis using parameters from a file
    -s seq_count      Sequence mode, using the result to seed the next step
    -j nb_jobs        Set number of simultanious jobs [1]
EOF

# we need arguments
if [[ $# -lt 1 ]]; then
  echo "$USAGE" >&2
  exit 1
fi

# default synthesis command
SYNTH_COMMAND="$(dirname $0)/checked_synth.sh"

# parsing the arguments
while getopts l:i:o:p:j:s:c:h OPT; do
  case $OPT in
  h)  echo "$USAGE"
      exit 0 ;;
  l)  if [[ ! -f "$OPTARG" ]]; then
        echo "Invalid option -j, $OPTARG is not a valid file." >& 2
        exit 1
      else
        echo "Using file $OPTARG as list."
        FILE_LIST=$(cat "$OPTARG")
      fi
      ;;
  i)  FILE_LIST="$OPTARG" ;;
  o)  BASE_DIR="$OPTARG" ;;
  p)  PARAM_FILE="$OPTARG" ;;
  j)  MAX_NPROC=$OPTARG ;;
  s)  SEQUENCE_COUNT="$OPTARG" ;;
  c)  SYNTH_COMMAND="$OPTARG" ;;
  \?) # invalid parameter
      echo "$USAGE" >&2
      exit 1 ;;
  esac
done
shift $((OPTIND-1))

# if we kill this process, children should all die!
trap "kill -- -$BASHPID" 2
trap "echo 'Requesting the end.'; LET_ME_DIE=1" 3

# parallelization functions
source func_parallel.sh

# output pid in case we need to get killed and it doesn't work properly
echo "PID: $$"

# the file list
if [[ -z "$FILE_LIST" ]]; then
  FILE_LIST=${1/:/ }
  shift
fi
for file in $FILE_LIST; do
  if [[ ! -f $file ]]; then
    echo "Invalid input file: $file" >&2
    exit 2
  fi
done

# default parameters
[[ -z "$BASE_DIR" ]] && BASE_DIR="Results/"
[[ ! -d "$BASE_DIR" ]] && mkdir -p "$BASE_DIR"
echo "$BASE_LINE" > "$BASE_DIR"/bash_call.dat


# are we using a param file?
if [[ -n "$PARAM_FILE" ]]; then
  PARAMS="$@"
  COUNT=$(cat "$PARAM_FILE" | wc -l)
  [[ $(echo "$PARAMS" | grep rand_seed) ]] && HAS_SEED=1 || HAS_SEED=0
  shift $#
else
  PARAMS=''
  COUNT=1
  HAS_SEED=0
fi

# parsing the arguments
KEYS=''
declare -A params
VARYING_PARAMS=1
while [[ $# -gt 1 ]]; do

  # is it a meta parameter?
  if [[ "$1" == '!' ]]; then
    VARYING_PARAMS=$((1-VARYING_PARAMS))
    shift
    continue
  fi

  echo "Params: $1 -> $2"

  # do we have the seed?
  if [[ "$1" = 'rand_seed' ]]; then
    HAS_SEED=1
  fi

  # store the parameter
  if [[ "$2" =~ , ]] && [[ $VARYING_PARAMS -eq 1 ]]; then
    # store in map
    seq_key="$1"
    seq_val=${2//,/ }
    echo "$seq_key = $seq_val"
    params["$seq_key"]="$seq_val"
    # store key
    KEYS="$KEYS $seq_key"
    # count elements (-1)
    seq_size=$(echo "$seq_val" | tr -cd ' ' | wc -c)
    # multiply full possibility count
    COUNT=$((COUNT * (seq_size + 1)))
  elif [[ "$2" =~ : ]] && [[ $VARYING_PARAMS -eq 1 ]]; then
    # store in map
    seq_key="$1"
    seq_val=$(seq -s ' ' $(echo "$2" | cut -d : -f 1-3 --output-delimiter ' '))
    echo "$seq_key = $seq_val"
    params["$seq_key"]="$seq_val"
    # store key
    KEYS="$KEYS $seq_key"
    # count elements (-1)
    seq_size=$(echo "$seq_val" | tr -cd ' ' | wc -c)
    # multiply full possibility count
    COUNT=$((COUNT * (seq_size + 1)))
  else
    # store directly in param list
    PARAMS="$PARAMS $1 $2"
  fi
  shift 2
done

# default seed
[[ $HAS_SEED -eq 0 ]] && PARAMS="rand_seed 0 $PARAMS"

function genParams() {
  local n=$(($1 - 1)) # 1-based argument, but 0-based for computations
  local str=''
  for key in $KEYS; do
    local vals=${params[$key]}
    local count=$(echo "$vals" | tr -cd ' ' | wc -c)
    count=$((count + 1))
    # find value index (0-based)
    local index=$((n % count))
    n=$(( (n - index) / count )) # integer division
    # pick the parameter value (1-based)
    index=$((1 + index))
    local value=$(echo "$vals" | cut -d ' ' -f $index)
    str="$str $key $value"
  done
  echo "$str"
}

function genName() {
  local parts="$1"
  local str=''
  for token in $parts; do
    if [[ "$token" =~ [0-9] ]]; then
      str="${str}${token}"
    else
      local key=''
      for noun in ${token//_/ }; do
        key="${key}${noun:0:1}"
      done
      if [[ ${#str} -gt 0 ]]; then
        str="${str}_${key}"
      else
        str="$key"
      fi
    fi
  done
  echo "$str"
}

# record when we start
start_time=$(date +%s)
[[ -z "$STARTUP_INDEX" ]] && STARTUP_INDEX=1
FULL_SEQ_COUNT=$COUNT
[[ -n "$SEQUENCE_COUNT" ]] && FULL_SEQ_COUNT=$((COUNT * SEQUENCE_COUNT)) || SEQUENCE_COUNT=1
# go over all parameter combinations
for j in $(seq "$STARTUP_INDEX" $FULL_SEQ_COUNT); do
  # indices
  i=$((j % COUNT)) # the parameter index
  [[ $i -eq 0 ]] && i=$COUNT
  SEQ_IDX=$(((j-1) / COUNT + 1))
  [[ $SEQUENCE_COUNT -gt 1 ]] && echo "Sequence #$SEQ_IDX"

  if [[ -n "$PARAM_FILE" ]]; then
    # use the ith line as parameter list
    VAR_PARAMS="$(sed -n $i'p' < $PARAM_FILE)"
    # check whether it starts with the directory name (@dirname params...)
    VAR_FIRST=$(cut -d' ' -f1 <<< "$VAR_PARAMS")
    if [[ "$VAR_FIRST" =~ ^@.+ ]]; then
      OUTPUT_DIR="$BASE_DIR/${VAR_FIRST:1}" # use @dir for dir name
      VAR_PARAMS="${VAR_PARAMS:${#VAR_FIRST}}" # remove dir name
    else
      OUTPUT_DIR="$BASE_DIR"/"p$i"
    fi
    [[ $SEQUENCE_COUNT -gt 1 ]] && OUTPUT_DIR="${OUTPUT_DIR}_i$SEQ_IDX"
  else
    # variable parameters
    VAR_PARAMS="$(genParams $i)"

    # output and parameter
    DIR_NAME=$(genName "$VAR_PARAMS")
    if [[ $SEQUENCE_COUNT -gt 1 ]]; then
      [[ -n "$DIR_NAME" ]] && DIR_NAME="${DIR_NAME}_i$SEQ_IDX" || DIR_NAME="i$SEQ_IDX"
    fi
    OUTPUT_DIR="$BASE_DIR/$DIR_NAME"
  fi

  # current effective parameters
  CUR_PARAMS="$PARAMS $VAR_PARAMS"

  # storing the parameters
  [[ ! -d "$OUTPUT_DIR" ]] && mkdir -p "$OUTPUT_DIR"
  echo "$CUR_PARAMS" > "$OUTPUT_DIR"/params.dat
  echo "- Parameters: $CUR_PARAMS"

  # generation for each texture
  for file in $FILE_LIST; do
    # the name
    RES_DIR=$(basename $file)
    RES_DIR="${RES_DIR/%.jpg/}" # strip file extension
    RES_DIR="${RES_DIR/%.png/}"

    # sequence file change
    if [[ $SEQ_IDX -gt 1 ]]; then
      PREV_DIR="${OUTPUT_DIR/%i$SEQ_IDX/i$((SEQ_IDX-1))}/$RES_DIR"
      # simlink for the name
      PWD=$(pwd)
      [[ ! -f "$PWD/$PREV_DIR/${RES_DIR}.png" ]] && ln -s "$PWD/$PREV_DIR/img.png" "$PWD/$PREV_DIR/${RES_DIR}.png"
      file="$PREV_DIR/${RES_DIR}.png" # use simlink to keep the correct name
    fi

    # work state (new, resume, skip)
    TEX_DIR="$OUTPUT_DIR/$RES_DIR"
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
    COMMAND="$SYNTH_COMMAND $file $OUTPUT_DIR num_threads 1 $CUR_PARAMS"
    if [[ "$MAX_NPROC" -le 1 ]]; then
      $COMMAND
    else
      waitForProcess
      [[ -n "$LET_ME_DIE" ]] && break
      mkdir -p "$TEX_DIR"
      $COMMAND 2> "$TEX_DIR"/proc_err.log > "$TEX_DIR"/proc.log &
    fi
  done
  [[ -n "$LET_ME_DIE" ]] && break
done

# wait for all the chidlren
echo "Waiting for the last synthesis"
wait
end_time=$(date +%s)
duration=$(($end_time - $start_time))
echo "Done. ($duration sec)"

# store time result
echo $duration > "$OUTPUT_DIR"/duration.dat
