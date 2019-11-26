#!/bin/bash
export TBBROOT="/users/sojavee/BayesRRcmd_genSparse/BayesRRcmd/tbb" #
tbb_bin="/users/sojavee/BayesRRcmd_genSparse/BayesRRcmd/src/tbb_cmake_build/tbb_cmake_build_subdir_debug" #
if [ -z "$CPATH" ]; then #
    export CPATH="${TBBROOT}/include" #
else #
    export CPATH="${TBBROOT}/include:$CPATH" #
fi #
if [ -z "$LIBRARY_PATH" ]; then #
    export LIBRARY_PATH="${tbb_bin}" #
else #
    export LIBRARY_PATH="${tbb_bin}:$LIBRARY_PATH" #
fi #
if [ -z "$LD_LIBRARY_PATH" ]; then #
    export LD_LIBRARY_PATH="${tbb_bin}" #
else #
    export LD_LIBRARY_PATH="${tbb_bin}:$LD_LIBRARY_PATH" #
fi #
 #
