#!/usr/bin/env bash

TARGET="$1"
shift

PARAMS=""
if [ $# -gt 1 ]; then
  PARAMS="$@"
fi
[[ -n "$TARGET" ]] && [[ "$TARGET" -eq 0 ]] && TARGET=""

# platform specific
os_name="$(uname)"
if [[ "$os_name" =~ Darwin ]]; then
  OPENMP_FLAG=""
  OPENMP_LIBS=""
else
  OPENMP_FLAG="$OPENMP_FLAG"
  OPENMP_LIBS="$OPENMP_LIBS"
fi

# general flags
OPTI_FLAGS="-O6 $OPENMP_FLAG"
BASE_FLAGS="$PARAMS -DNDEBUG -DUNIX_MODE -DMEXMODE $OPENMP_FLAG"
LIBS_FLAGS="-Wl,-e,mexFunction $OPENMP_LIBS"
LAST_FLAGS="-I./include/ -I/usr/include/"

function mex_call {
  echo mex "CXXOPTIM_FLAGS='$OPTI_FLAGS'" "CXXFLAGS='$BASE_FLAGS'" "CXXLIBS='\${CXXLIBS} $LIBS_FLAGS'" $LAST_FLAGS $@
  mex CXXOPTIM_FLAGS="$OPTI_FLAGS" CXXFLAGS="$BASE_FLAGS" CXXLIBS="\${CXXLIBS} $LIBS_FLAGS" $LAST_FLAGS $@
}
MEX_CALL="mex_call" # default

# os-dependent flags
case "$os_name" in
  *Linux*)
  OPTI_FLAGS="$OPTI_FLAGS -w -s -ffast-math -fomit-frame-pointer -fstrength-reduce -msse2 -funroll-loops -fPIC"
  BASE_FLAGS="$BASE_FLAGS -fPIC -ftls-model=global-dynamic"
  LIBS_FLAGS="-Wl,--export-dynamic $LIBS_FLAGS -shared"
  ;;
  *Darwin*)
  for d in /Applications/MATLAB*; do 
    MATLAB_LOCATION="$d"
    break
  done
  OPTI_FLAGS="$OPTI_FLAGS -w -s -ffast-math -fomit-frame-pointer -fstrength-reduce -msse2 -funroll-loops -fPIC"
  BASE_FLAGS="$BASE_FLAGS -fPIC -ftls-model=global-dynamic"
  LIBS_FLAGS="$LIBS_FLAGS -dynamiclib"
  LAST_FLAGS="-I/usr/local/Cellar/boost/1.49.0/include -L$MATLAB_LOCATION/bin/maci64 -L/usr/lib $LAST_FLAGS"
  ;;
  \?)
  MEX_CALL="echo"
  ;;
esac

# set output directory
[[ -z "$ODIR" ]] && ODIR=dist

# Optimized (slow to build)
if [[ -z "$TARGET" || "$TARGET" -eq 1 ]]; then 
  $MEX_CALL src/nnfMex.cpp -output "$ODIR/nnmex" -output "$ODIR/nnmex" &
fi
if [[ -z "$TARGET" || "$TARGET" -eq 2 ]]; then
  $MEX_CALL src/voteMex.cpp -output "$ODIR/votemex" -output "$ODIR/votemex" &
fi

if [[ -z "$TARGET" || "$TARGET" -eq 3 ]]; then
  $MEX_CALL src/segmentMex.cpp -output "$ODIR/segmentmex" -output "$ODIR/segmentmex" &
fi

if [[ -z "$TARGET" || "$TARGET" -eq 4 ]]; then
  $MEX_CALL src/ggdtMex.cpp -output "$ODIR/ggdtmex" -output "$ODIR/ggdtmex" &
fi

if [[ -z "$TARGET" || "$TARGET" -eq 5 ]]; then
  $MEX_CALL src/latticeMex.cpp -output "$ODIR/latticemex" -output "$ODIR/latticemex" &
fi

if [[ -z "$TARGET" || "$TARGET" -eq 6 ]]; then
  $MEX_CALL src/colorHistMex.cpp -output "$ODIR/colorhistmex" -output "$ODIR/colorhistmex" &
fi

if [[ "$TARGET" -eq 7 ]]; then # not officially compiled
  $MEX_CALL src/jigsawMex.cpp -output "$ODIR/jigsawmex" -output "$ODIR/jigsawmex" &
fi

wait
